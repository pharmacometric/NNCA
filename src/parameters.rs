use crate::{models::*, errors::NcaError, Result};
use nalgebra::{DMatrix, DVector};
use statrs::statistics::Statistics;

pub struct ParameterCalculator;

impl ParameterCalculator {
    /// Calculate terminal elimination rate constant (lambda_z)
    pub fn calculate_lambda_z(
        observations: &[Observation],
        selection: &LambdaZSelection,
    ) -> Result<(f64, f64, Vec<usize>)> {
        match selection {
            LambdaZSelection::Auto => Self::auto_lambda_z_selection(observations),
            LambdaZSelection::Manual(indices) => Self::manual_lambda_z_selection(observations, indices),
            LambdaZSelection::BestFit { min_points, r_squared_threshold } => {
                Self::best_fit_lambda_z_selection(observations, *min_points, *r_squared_threshold)
            }
        }
    }

    fn auto_lambda_z_selection(observations: &[Observation]) -> Result<(f64, f64, Vec<usize>)> {
        let n = observations.len();
        if n < 3 {
            return Err(NcaError::InsufficientData(
                "Need at least 3 points for lambda_z calculation".to_string()
            ));
        }

        let mut best_r_squared = 0.0;
        let mut best_lambda_z = 0.0;
        let mut best_indices = Vec::new();

        // Try different combinations of the last points
        for start_idx in 0..=(n.saturating_sub(3)) {
            let indices: Vec<usize> = (start_idx..n).collect();
            
            if let Ok((lambda_z, r_squared)) = Self::fit_lambda_z(observations, &indices) {
                if r_squared > best_r_squared && r_squared >= 0.8 {
                    best_r_squared = r_squared;
                    best_lambda_z = lambda_z;
                    best_indices = indices;
                }
            }
        }

        if best_indices.is_empty() {
            return Err(NcaError::CalculationError(
                "Could not find suitable points for lambda_z calculation".to_string()
            ));
        }

        Ok((best_lambda_z, best_r_squared, best_indices))
    }

    fn manual_lambda_z_selection(
        observations: &[Observation],
        indices: &[usize],
    ) -> Result<(f64, f64, Vec<usize>)> {
        let (lambda_z, r_squared) = Self::fit_lambda_z(observations, indices)?;
        Ok((lambda_z, r_squared, indices.to_vec()))
    }

    fn best_fit_lambda_z_selection(
        observations: &[Observation],
        min_points: usize,
        r_squared_threshold: f64,
    ) -> Result<(f64, f64, Vec<usize>)> {
        let n = observations.len();
        if n < min_points {
            return Err(NcaError::InsufficientData(
                format!("Need at least {} points for lambda_z calculation", min_points)
            ));
        }

        let mut best_r_squared = 0.0;
        let mut best_lambda_z = 0.0;
        let mut best_indices = Vec::new();

        // Try different combinations of points
        for start_idx in 0..=(n.saturating_sub(min_points)) {
            for end_idx in (start_idx + min_points - 1)..n {
                let indices: Vec<usize> = (start_idx..=end_idx).collect();
                
                if let Ok((lambda_z, r_squared)) = Self::fit_lambda_z(observations, &indices) {
                    if r_squared > best_r_squared && r_squared >= r_squared_threshold {
                        best_r_squared = r_squared;
                        best_lambda_z = lambda_z;
                        best_indices = indices;
                    }
                }
            }
        }

        if best_indices.is_empty() {
            return Err(NcaError::CalculationError(
                format!("Could not find suitable points with R² >= {}", r_squared_threshold)
            ));
        }

        Ok((best_lambda_z, best_r_squared, best_indices))
    }

    fn fit_lambda_z(observations: &[Observation], indices: &[usize]) -> Result<(f64, f64)> {
        let mut times = Vec::new();
        let mut ln_concentrations = Vec::new();

        for &idx in indices {
            if idx >= observations.len() {
                continue;
            }
            
            let obs = &observations[idx];
            if obs.concentration > 0.0 {
                times.push(obs.time);
                ln_concentrations.push(obs.concentration.ln());
            }
        }

        if times.len() < 2 {
            return Err(NcaError::InsufficientData(
                "Need at least 2 positive concentrations for lambda_z".to_string()
            ));
        }

        // Linear regression: ln(C) = ln(C0) - lambda_z * t
        let n = times.len() as f64;
        let sum_t = times.iter().sum::<f64>();
        let sum_ln_c = ln_concentrations.iter().sum::<f64>();
        let sum_t_ln_c = times.iter().zip(&ln_concentrations)
            .map(|(t, ln_c)| t * ln_c)
            .sum::<f64>();
        let sum_t2 = times.iter().map(|t| t * t).sum::<f64>();

        let slope = (n * sum_t_ln_c - sum_t * sum_ln_c) / (n * sum_t2 - sum_t * sum_t);
        let lambda_z = -slope; // Negative because we're fitting declining concentrations

        // Calculate R-squared
        let mean_ln_c = sum_ln_c / n;
        let ss_tot = ln_concentrations.iter()
            .map(|ln_c| (ln_c - mean_ln_c).powi(2))
            .sum::<f64>();
        
        let intercept = (sum_ln_c - slope * sum_t) / n;
        let ss_res = times.iter().zip(&ln_concentrations)
            .map(|(t, ln_c)| {
                let predicted = intercept + slope * t;
                (ln_c - predicted).powi(2)
            })
            .sum::<f64>();

        let r_squared = if ss_tot > 0.0 { 1.0 - (ss_res / ss_tot) } else { 0.0 };

        Ok((lambda_z, r_squared))
    }

    /// Calculate Cmax and Tmax
    pub fn calculate_cmax_tmax(observations: &[Observation]) -> Result<(f64, f64)> {
        let max_obs = observations
            .iter()
            .max_by(|a, b| a.concentration.partial_cmp(&b.concentration).unwrap())
            .ok_or_else(|| NcaError::InsufficientData("No observations available".to_string()))?;

        Ok((max_obs.concentration, max_obs.time))
    }

    /// Calculate half-life from lambda_z
    pub fn calculate_half_life(lambda_z: f64) -> Result<f64> {
        if lambda_z <= 0.0 {
            return Err(NcaError::CalculationError(
                "Lambda_z must be positive for half-life calculation".to_string()
            ));
        }
        
        Ok(0.693147 / lambda_z) // ln(2) / lambda_z
    }

    /// Calculate clearance for IV dosing
    pub fn calculate_clearance_iv(dose: f64, auc_inf: f64) -> Result<f64> {
        if auc_inf <= 0.0 {
            return Err(NcaError::CalculationError(
                "AUC_inf must be positive for clearance calculation".to_string()
            ));
        }
        
        Ok(dose / auc_inf)
    }

    /// Calculate apparent clearance for oral dosing
    pub fn calculate_clearance_oral(dose: f64, auc_inf: f64, bioavailability: Option<f64>) -> Result<f64> {
        let cl_f = Self::calculate_clearance_iv(dose, auc_inf)?;
        
        match bioavailability {
            Some(f) if f > 0.0 => Ok(cl_f / f),
            _ => Ok(cl_f), // Return CL/F if bioavailability unknown
        }
    }

    /// Calculate volume of distribution at steady state
    pub fn calculate_vss(clearance: f64, mrt: f64) -> Result<f64> {
        if clearance <= 0.0 || mrt <= 0.0 {
            return Err(NcaError::CalculationError(
                "Clearance and MRT must be positive for Vss calculation".to_string()
            ));
        }
        
        Ok(clearance * mrt)
    }

    /// Calculate terminal volume of distribution
    pub fn calculate_vz(clearance: f64, lambda_z: f64) -> Result<f64> {
        if clearance <= 0.0 || lambda_z <= 0.0 {
            return Err(NcaError::CalculationError(
                "Clearance and lambda_z must be positive for Vz calculation".to_string()
            ));
        }
        
        Ok(clearance / lambda_z)
    }

    /// Calculate mean residence time
    pub fn calculate_mrt(aumc_inf: f64, auc_inf: f64) -> Result<f64> {
        if auc_inf <= 0.0 {
            return Err(NcaError::CalculationError(
                "AUC_inf must be positive for MRT calculation".to_string()
            ));
        }
        
        Ok(aumc_inf / auc_inf)
    }

    /// Find time of last quantifiable concentration
    pub fn find_tlast_clast(observations: &[Observation]) -> Option<(f64, f64)> {
        observations
            .iter()
            .rev()
            .find(|obs| obs.concentration > 0.0 && !obs.bloq)
            .map(|obs| (obs.time, obs.concentration))
    }

    /// Calculate percentage of AUC extrapolated to infinity
    pub fn calculate_auc_percent_extrap(auc_last: f64, auc_inf: f64) -> Result<f64> {
        if auc_inf <= 0.0 {
            return Err(NcaError::CalculationError(
                "AUC_inf must be positive".to_string()
            ));
        }
        
        Ok(((auc_inf - auc_last) / auc_inf) * 100.0)
    }
}