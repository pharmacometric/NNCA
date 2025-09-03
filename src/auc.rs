use crate::{models::*, errors::NcaError, Result};
use std::collections::HashMap;

pub struct AucCalculator;

impl AucCalculator {
    /// Calculate AUC using multiple methods for comparison
    pub fn calculate_all_methods(
        observations: &[Observation],
        config: &AnalysisConfig,
    ) -> Result<HashMap<String, f64>> {
        let mut results = HashMap::new();
        
        // Filter valid observations (remove BLQ based on config)
        let filtered_obs = Self::filter_observations(observations, &config.lloq_handling);
        
        if filtered_obs.len() < 2 {
            return Err(NcaError::InsufficientData(
                "Need at least 2 data points for AUC calculation".to_string()
            ));
        }

        // Linear trapezoidal method
        results.insert(
            "linear_trapezoidal".to_string(),
            Self::linear_trapezoidal(&filtered_obs)?,
        );

        // Log trapezoidal method
        results.insert(
            "log_trapezoidal".to_string(),
            Self::log_trapezoidal(&filtered_obs)?,
        );

        // Linear-log trapezoidal method (Phoenix WinNonlin style)
        results.insert(
            "linear_log_trapezoidal".to_string(),
            Self::linear_log_trapezoidal(&filtered_obs)?,
        );

        // Linear up log down method
        results.insert(
            "linear_up_log_down".to_string(),
            Self::linear_up_log_down(&filtered_obs)?,
        );

        Ok(results)
    }

    fn filter_observations(observations: &[Observation], lloq_handling: &LloqHandling) -> Vec<Observation> {
        observations
            .iter()
            .filter_map(|obs| {
                if obs.bloq {
                    match lloq_handling {
                        LloqHandling::Drop => None,
                        LloqHandling::Zero => {
                            let mut modified_obs = obs.clone();
                            modified_obs.concentration = 0.0;
                            Some(modified_obs)
                        }
                        LloqHandling::HalfLloq => {
                            let mut modified_obs = obs.clone();
                            modified_obs.concentration = obs.lloq.unwrap_or(0.0) / 2.0;
                            Some(modified_obs)
                        }
                    }
                } else {
                    Some(obs.clone())
                }
            })
            .collect()
    }

    fn linear_trapezoidal(observations: &[Observation]) -> Result<f64> {
        let mut auc = 0.0;
        
        for i in 1..observations.len() {
            let t1 = observations[i - 1].time;
            let t2 = observations[i].time;
            let c1 = observations[i - 1].concentration;
            let c2 = observations[i].concentration;
            
            if t2 <= t1 {
                continue;
            }
            
            // Linear trapezoidal rule
            auc += (t2 - t1) * (c1 + c2) / 2.0;
        }
        
        Ok(auc)
    }

    fn log_trapezoidal(observations: &[Observation]) -> Result<f64> {
        let mut auc = 0.0;
        
        for i in 1..observations.len() {
            let t1 = observations[i - 1].time;
            let t2 = observations[i].time;
            let c1 = observations[i - 1].concentration;
            let c2 = observations[i].concentration;
            
            if t2 <= t1 || c1 <= 0.0 || c2 <= 0.0 {
                continue;
            }
            
            // Log trapezoidal rule
            let ln_c1 = c1.ln();
            let ln_c2 = c2.ln();
            
            if (ln_c1 - ln_c2).abs() < 1e-10 {
                // Concentrations are essentially equal, use linear
                auc += (t2 - t1) * (c1 + c2) / 2.0;
            } else {
                auc += (t2 - t1) * (c1 - c2) / (ln_c1 - ln_c2);
            }
        }
        
        Ok(auc)
    }

    fn linear_log_trapezoidal(observations: &[Observation]) -> Result<f64> {
        let mut auc = 0.0;
        
        for i in 1..observations.len() {
            let t1 = observations[i - 1].time;
            let t2 = observations[i].time;
            let c1 = observations[i - 1].concentration;
            let c2 = observations[i].concentration;
            
            if t2 <= t1 {
                continue;
            }
            
            // Use log trapezoidal if both concentrations > 0 and declining
            // Otherwise use linear trapezoidal
            if c1 > 0.0 && c2 > 0.0 && c2 < c1 {
                let ln_c1 = c1.ln();
                let ln_c2 = c2.ln();
                
                if (ln_c1 - ln_c2).abs() < 1e-10 {
                    auc += (t2 - t1) * (c1 + c2) / 2.0;
                } else {
                    auc += (t2 - t1) * (c1 - c2) / (ln_c1 - ln_c2);
                }
            } else {
                // Linear trapezoidal
                auc += (t2 - t1) * (c1 + c2) / 2.0;
            }
        }
        
        Ok(auc)
    }

    fn linear_up_log_down(observations: &[Observation]) -> Result<f64> {
        let mut auc = 0.0;
        
        for i in 1..observations.len() {
            let t1 = observations[i - 1].time;
            let t2 = observations[i].time;
            let c1 = observations[i - 1].concentration;
            let c2 = observations[i].concentration;
            
            if t2 <= t1 {
                continue;
            }
            
            // Use linear when concentration is increasing, log when decreasing
            if c2 >= c1 {
                // Linear trapezoidal for increasing concentrations
                auc += (t2 - t1) * (c1 + c2) / 2.0;
            } else if c1 > 0.0 && c2 > 0.0 {
                // Log trapezoidal for decreasing concentrations
                let ln_c1 = c1.ln();
                let ln_c2 = c2.ln();
                
                if (ln_c1 - ln_c2).abs() < 1e-10 {
                    auc += (t2 - t1) * (c1 + c2) / 2.0;
                } else {
                    auc += (t2 - t1) * (c1 - c2) / (ln_c1 - ln_c2);
                }
            } else {
                // Fall back to linear if log calculation isn't possible
                auc += (t2 - t1) * (c1 + c2) / 2.0;
            }
        }
        
        Ok(auc)
    }

    /// Calculate AUC to infinity using terminal elimination rate constant
    pub fn calculate_auc_inf(
        auc_last: f64,
        clast: f64,
        lambda_z: f64,
    ) -> Result<f64> {
        if lambda_z <= 0.0 {
            return Err(NcaError::CalculationError(
                "Lambda_z must be positive for AUC infinity calculation".to_string()
            ));
        }
        
        let auc_extrap = clast / lambda_z;
        Ok(auc_last + auc_extrap)
    }

    /// Calculate AUMC (Area Under Moment Curve)
    pub fn calculate_aumc(observations: &[Observation]) -> Result<f64> {
        let mut aumc = 0.0;
        
        for i in 1..observations.len() {
            let t1 = observations[i - 1].time;
            let t2 = observations[i].time;
            let c1 = observations[i - 1].concentration;
            let c2 = observations[i].concentration;
            
            if t2 <= t1 {
                continue;
            }
            
            // AUMC calculation using linear trapezoidal rule
            aumc += (t2 - t1) * (t1 * c1 + t2 * c2) / 2.0;
        }
        
        Ok(aumc)
    }

    /// Calculate AUMC to infinity
    pub fn calculate_aumc_inf(
        aumc_last: f64,
        tlast: f64,
        clast: f64,
        lambda_z: f64,
    ) -> Result<f64> {
        if lambda_z <= 0.0 {
            return Err(NcaError::CalculationError(
                "Lambda_z must be positive for AUMC infinity calculation".to_string()
            ));
        }
        
        let aumc_extrap = (tlast * clast / lambda_z) + (clast / (lambda_z * lambda_z));
        Ok(aumc_last + aumc_extrap)
    }
}