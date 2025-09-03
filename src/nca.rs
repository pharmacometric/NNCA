use crate::{models::*, parameters::ParameterCalculator, auc::AucCalculator, Result};
use std::collections::HashMap;

pub struct NcaAnalyzer;

impl NcaAnalyzer {
    /// Perform comprehensive NCA analysis for a single subject
    pub fn analyze_subject(
        subject: &Subject,
        config: &AnalysisConfig,
    ) -> Result<(NcaResults, Vec<String>)> {
        let observations = &subject.observations;
        let mut warnings = Vec::new();
        
        if observations.is_empty() {
            return Err(crate::errors::NcaError::InsufficientData(
                "No observations available for analysis".to_string()
            ));
        }

        // Sort observations by time
        let mut sorted_obs = observations.clone();
        sorted_obs.sort_by(|a, b| a.time.partial_cmp(&b.time).unwrap());

        // Check minimum quantifiable concentrations requirement
        let quantifiable_count = sorted_obs.iter()
            .filter(|obs| obs.concentration > 0.0 && !obs.bloq)
            .count();
        
        if quantifiable_count < 3 {
            return Err(crate::errors::NcaError::InsufficientData(
                format!("Subject {} has only {} quantifiable concentrations (minimum 3 required)", 
                    subject.id, quantifiable_count)
            ));
        }

        // Calculate primary parameters
        let individual_params = Self::calculate_individual_parameters(&sorted_obs, subject, config)?;
        
        // Calculate using all AUC methods for comparison
        let mut method_comparisons = HashMap::new();
        
        for auc_method in &config.auc_methods {
            let method_name = format!("{:?}", auc_method);
            let method_config = AnalysisConfig {
                auc_methods: vec![auc_method.clone()],
                ..config.clone()
            };
            
            if let Ok(params) = Self::calculate_individual_parameters(&sorted_obs, subject, &method_config) {
                method_comparisons.insert(method_name, params);
            }
        }

        let results = NcaResults {
            subject_id: subject.id.clone(),
            individual_parameters: individual_params,
            method_comparisons,
        };

        // Generate warnings for missing parameters
        let param_warnings = Self::check_parameter_completeness(&results);
        warnings.extend(param_warnings);

        Ok((results, warnings))
    }

    fn calculate_individual_parameters(
        observations: &[Observation],
        subject: &Subject,
        config: &AnalysisConfig,
    ) -> Result<IndividualParameters> {
        // Basic parameters
        let (cmax, tmax) = ParameterCalculator::calculate_cmax_tmax(observations)?;
        let (tlast, clast) = ParameterCalculator::find_tlast_clast(observations)
            .ok_or_else(|| crate::errors::NcaError::InsufficientData(
                "No quantifiable concentrations found".to_string()
            ))?;

        // AUC calculations
        let auc_methods = AucCalculator::calculate_all_methods(observations, config)?;
        let auc_last = auc_methods.get("linear_trapezoidal").copied()
            .or_else(|| auc_methods.values().next().copied())
            .unwrap_or(0.0);

        // Terminal elimination parameters
        let (lambda_z, lambda_z_r_squared, _) = ParameterCalculator::calculate_lambda_z(
            observations,
            &config.lambda_z_selection,
        ).unwrap_or((0.0, 0.0, Vec::new()));

        // Calculate AUC to infinity
        let (auc_inf, auc_inf_pred) = if lambda_z > 0.0 {
            let auc_inf = AucCalculator::calculate_auc_inf(auc_last, clast, lambda_z)?;
            (Some(auc_inf), Some(auc_inf))
        } else {
            (None, None)
        };

        // AUC extrapolation percentage
        let auc_percent_extrap = if let Some(auc_inf_val) = auc_inf {
            Some(ParameterCalculator::calculate_auc_percent_extrap(auc_last, auc_inf_val)?)
        } else {
            None
        };

        // AUMC calculations
        let aumc_last = AucCalculator::calculate_aumc(observations)?;
        let aumc_inf = if lambda_z > 0.0 {
            Some(AucCalculator::calculate_aumc_inf(aumc_last, tlast, clast, lambda_z)?)
        } else {
            None
        };

        // Half-life
        let half_life = if lambda_z > 0.0 {
            Some(ParameterCalculator::calculate_half_life(lambda_z)?)
        } else {
            None
        };

        // MRT
        let mrt = if let (Some(aumc_inf_val), Some(auc_inf_val)) = (aumc_inf, auc_inf) {
            Some(ParameterCalculator::calculate_mrt(aumc_inf_val, auc_inf_val)?)
        } else {
            None
        };

        // Clearance and volume calculations
        let total_dose = Self::calculate_total_dose(subject);
        let (clearance, volume_steady_state, volume_terminal) = 
            Self::calculate_clearance_and_volumes(total_dose, auc_inf, lambda_z, mrt)?;

        Ok(IndividualParameters {
            auc_last: Some(auc_last),
            auc_inf,
            auc_inf_pred,
            auc_percent_extrap,
            aumc_last: Some(aumc_last),
            aumc_inf,
            cmax: Some(cmax),
            tmax: Some(tmax),
            tlast: Some(tlast),
            clast: Some(clast),
            half_life,
            lambda_z: if lambda_z > 0.0 { Some(lambda_z) } else { None },
            lambda_z_r_squared: if lambda_z_r_squared > 0.0 { Some(lambda_z_r_squared) } else { None },
            clearance,
            volume_steady_state,
            volume_terminal,
            mrt,
            bioavailability: None, // Would need reference data
        })
    }

    fn calculate_total_dose(subject: &Subject) -> f64 {
        subject.dosing_events.iter().map(|dose| dose.dose).sum()
    }

    fn calculate_clearance_and_volumes(
        total_dose: f64,
        auc_inf: Option<f64>,
        lambda_z: f64,
        mrt: Option<f64>,
    ) -> Result<(Option<f64>, Option<f64>, Option<f64>)> {
        let clearance = if let Some(auc_inf_val) = auc_inf {
            if auc_inf_val > 0.0 {
                Some(ParameterCalculator::calculate_clearance_iv(total_dose, auc_inf_val)?)
            } else {
                None
            }
        } else {
            None
        };

        let volume_steady_state = if let (Some(cl), Some(mrt_val)) = (clearance, mrt) {
            Some(ParameterCalculator::calculate_vss(cl, mrt_val)?)
        } else {
            None
        };

        let volume_terminal = if let Some(cl) = clearance {
            if lambda_z > 0.0 {
                Some(ParameterCalculator::calculate_vz(cl, lambda_z)?)
            } else {
                None
            }
        } else {
            None
        };

        Ok((clearance, volume_steady_state, volume_terminal))
    }

    fn check_parameter_completeness(results: &NcaResults) -> Vec<String> {
        let mut warnings = Vec::new();
        let params = &results.individual_parameters;
        
        if params.auc_inf.is_none() {
            warnings.push("AUC_inf could not be calculated - insufficient terminal phase data".to_string());
        }
        
        if params.lambda_z.is_none() {
            warnings.push("Lambda_z could not be calculated - poor terminal phase fit".to_string());
        }
        
        if params.half_life.is_none() {
            warnings.push("Half-life could not be calculated - lambda_z unavailable".to_string());
        }
        
        if params.clearance.is_none() {
            warnings.push("Clearance could not be calculated - AUC_inf unavailable".to_string());
        }
        
        if params.mrt.is_none() {
            warnings.push("MRT could not be calculated - AUMC_inf or AUC_inf unavailable".to_string());
        }
        
        if let Some(extrap) = params.auc_percent_extrap {
            if extrap > 20.0 {
                warnings.push(format!("High AUC extrapolation ({}%) - results may be unreliable", extrap));
            }
        }
        
        if let Some(r_sq) = params.lambda_z_r_squared {
            if r_sq < 0.8 {
                warnings.push(format!("Poor terminal phase fit (R² = {:.3}) - lambda_z may be unreliable", r_sq));
            }
        }
        
        warnings
    }

    /// Validate analysis results for quality control
    pub fn validate_results(results: &NcaResults) -> Vec<String> {
        let mut warnings = Vec::new();
        let params = &results.individual_parameters;

        // Check AUC extrapolation
        if let Some(extrap) = params.auc_percent_extrap {
            if extrap > 20.0 {
                warnings.push(format!(
                    "High AUC extrapolation ({}%) for subject {}", 
                    extrap, results.subject_id
                ));
            }
        }

        // Check R-squared for lambda_z
        if let Some(r_sq) = params.lambda_z_r_squared {
            if r_sq < 0.8 {
                warnings.push(format!(
                    "Poor terminal phase fit (R² = {:.3}) for subject {}", 
                    r_sq, results.subject_id
                ));
            }
        }

        // Check for reasonable half-life values
        if let Some(t_half) = params.half_life {
            if t_half < 0.1 || t_half > 1000.0 {
                warnings.push(format!(
                    "Unusual half-life ({:.3} h) for subject {}", 
                    t_half, results.subject_id
                ));
            }
        }

        warnings
    }
}