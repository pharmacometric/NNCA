use crate::{models::*, Result};
use std::collections::HashMap;
use statrs::statistics::Statistics;

pub struct CovariateAnalyzer;

impl CovariateAnalyzer {
    /// Perform comprehensive covariate analysis
    pub fn analyze_covariates(
        results: &[NcaResults],
        subjects: &[Subject],
    ) -> Result<CovariateAnalysis> {
        let correlations = Self::calculate_covariate_correlations(results, subjects)?;
        let regression_analysis = Self::perform_regression_analysis(results, subjects)?;
        let dose_normalized_analysis = Self::perform_dose_normalization_analysis(results, subjects)?;

        Ok(CovariateAnalysis {
            correlations,
            regression_analysis,
            dose_normalized_analysis: Some(dose_normalized_analysis),
        })
    }

    fn calculate_covariate_correlations(
        results: &[NcaResults],
        subjects: &[Subject],
    ) -> Result<HashMap<String, CovariateCorrelation>> {
        let mut correlations = HashMap::new();
        
        let covariates = vec!["age", "weight", "height"];
        let parameters = vec!["auc_inf", "cmax", "clearance", "half_life", "volume_terminal"];

        for covariate in &covariates {
            let mut parameter_correlations = HashMap::new();
            let mut p_values = HashMap::new();

            for parameter in &parameters {
                let (covariate_values, parameter_values) = Self::extract_paired_values(
                    results, subjects, covariate, parameter
                );

                if covariate_values.len() >= 3 {
                    let correlation = Self::calculate_pearson_correlation(&covariate_values, &parameter_values);
                    let p_value = Self::correlation_p_value(correlation, covariate_values.len());
                    
                    parameter_correlations.insert(parameter.to_string(), correlation);
                    p_values.insert(parameter.to_string(), p_value);
                }
            }

            if !parameter_correlations.is_empty() {
                correlations.insert(
                    covariate.to_string(),
                    CovariateCorrelation {
                        covariate_name: covariate.to_string(),
                        parameter_correlations,
                        p_values,
                    },
                );
            }
        }

        Ok(correlations)
    }

    fn extract_paired_values(
        results: &[NcaResults],
        subjects: &[Subject],
        covariate: &str,
        parameter: &str,
    ) -> (Vec<f64>, Vec<f64>) {
        let mut covariate_values = Vec::new();
        let mut parameter_values = Vec::new();

        for (result, subject) in results.iter().zip(subjects.iter()) {
            let cov_value = match covariate {
                "age" => subject.demographics.age,
                "weight" => subject.demographics.weight,
                "height" => subject.demographics.height,
                _ => None,
            };

            let param_value = match parameter {
                "auc_inf" => result.individual_parameters.auc_inf,
                "cmax" => result.individual_parameters.cmax,
                "clearance" => result.individual_parameters.clearance,
                "half_life" => result.individual_parameters.half_life,
                "volume_terminal" => result.individual_parameters.volume_terminal,
                _ => None,
            };

            if let (Some(cov), Some(param)) = (cov_value, param_value) {
                covariate_values.push(cov);
                parameter_values.push(param);
            }
        }

        (covariate_values, parameter_values)
    }

    fn calculate_pearson_correlation(x: &[f64], y: &[f64]) -> f64 {
        if x.len() != y.len() || x.len() < 2 {
            return 0.0;
        }

        let mean_x = x.mean();
        let mean_y = y.mean();

        let numerator: f64 = x.iter().zip(y.iter())
            .map(|(xi, yi)| (xi - mean_x) * (yi - mean_y))
            .sum();

        let sum_sq_x: f64 = x.iter().map(|xi| (xi - mean_x).powi(2)).sum();
        let sum_sq_y: f64 = y.iter().map(|yi| (yi - mean_y).powi(2)).sum();

        let denominator = (sum_sq_x * sum_sq_y).sqrt();

        if denominator == 0.0 { 0.0 } else { numerator / denominator }
    }

    fn correlation_p_value(r: f64, n: usize) -> f64 {
        if n < 3 { return 1.0; }
        
        let df = n - 2;
        let t_stat = r * ((df as f64) / (1.0 - r * r)).sqrt();
        
        // Simplified p-value calculation
        2.0 * (1.0 - Self::t_cdf(t_stat.abs(), df as f64))
    }

    fn t_cdf(t: f64, df: f64) -> f64 {
        // Simplified t-distribution CDF approximation
        if df <= 0.0 { return 0.5; }
        
        // For large df, approximate with normal distribution
        if df > 30.0 {
            return Self::standard_normal_cdf(t);
        }
        
        // Rough approximation for small df
        0.5 + 0.5 * (t / (1.0 + t * t / df).sqrt()).tanh()
    }

    fn standard_normal_cdf(z: f64) -> f64 {
        0.5 * (1.0 + Self::erf(z / 2.0_f64.sqrt()))
    }

    fn erf(x: f64) -> f64 {
        // Error function approximation
        let a1 = 0.254829592;
        let a2 = -0.284496736;
        let a3 = 1.421413741;
        let a4 = -1.453152027;
        let a5 = 1.061405429;
        let p = 0.3275911;

        let sign = if x < 0.0 { -1.0 } else { 1.0 };
        let x = x.abs();

        let t = 1.0 / (1.0 + p * x);
        let y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * (-x * x).exp();

        sign * y
    }

    fn perform_regression_analysis(
        results: &[NcaResults],
        subjects: &[Subject],
    ) -> Result<HashMap<String, RegressionResults>> {
        let mut regression_results = HashMap::new();
        
        let covariates = vec!["age", "weight", "height"];
        let parameters = vec!["auc_inf", "cmax", "clearance"];

        for covariate in &covariates {
            for parameter in &parameters {
                let (x_values, y_values) = Self::extract_paired_values(
                    results, subjects, covariate, parameter
                );

                if x_values.len() >= 3 {
                    let mut regression = Self::simple_linear_regression(&x_values, &y_values);
                    regression.parameter = parameter.to_string();
                    regression.covariate = covariate.to_string();
                    regression_results.insert(
                        format!("{}_{}", parameter, covariate),
                        regression,
                    );
                }
            }
        }

        Ok(regression_results)
    }

    fn simple_linear_regression(x: &[f64], y: &[f64]) -> RegressionResults {
        if x.len() != y.len() || x.len() < 2 {
            return RegressionResults {
                parameter: "unknown".to_string(),
                covariate: "unknown".to_string(),
                slope: 0.0,
                intercept: 0.0,
                r_squared: 0.0,
                p_value: 1.0,
                confidence_interval: (0.0, 0.0),
            };
        }

        let n = x.len() as f64;
        let mean_x = x.mean();
        let mean_y = y.mean();

        let numerator: f64 = x.iter().zip(y.iter())
            .map(|(xi, yi)| (xi - mean_x) * (yi - mean_y))
            .sum();

        let denominator: f64 = x.iter()
            .map(|xi| (xi - mean_x).powi(2))
            .sum();

        let slope = if denominator != 0.0 { numerator / denominator } else { 0.0 };
        let intercept = mean_y - slope * mean_x;

        // Calculate R-squared
        let ss_tot: f64 = y.iter().map(|yi| (yi - mean_y).powi(2)).sum();
        let ss_res: f64 = x.iter().zip(y.iter())
            .map(|(xi, yi)| {
                let predicted = intercept + slope * xi;
                (yi - predicted).powi(2)
            })
            .sum();

        let r_squared = if ss_tot != 0.0 { 1.0 - (ss_res / ss_tot) } else { 0.0 };

        // Calculate standard error and p-value for slope
        let mse = if n > 2.0 { ss_res / (n - 2.0) } else { 0.0 };
        let se_slope = if denominator > 0.0 && mse > 0.0 {
            (mse / denominator).sqrt()
        } else {
            0.0
        };

        let t_stat = if se_slope > 0.0 { slope / se_slope } else { 0.0 };
        let p_value = if n > 2.0 {
            2.0 * (1.0 - Self::t_cdf(t_stat.abs(), n - 2.0))
        } else {
            1.0
        };

        // 95% confidence interval for slope
        let t_critical = 1.96; // Approximate for large samples
        let margin_error = t_critical * se_slope;
        let confidence_interval = (slope - margin_error, slope + margin_error);

        RegressionResults {
            parameter: "parameter".to_string(),
            covariate: "covariate".to_string(),
            slope,
            intercept,
            r_squared,
            p_value,
            confidence_interval,
        }
    }

    fn perform_dose_normalization_analysis(
        results: &[NcaResults],
        subjects: &[Subject],
    ) -> Result<DoseNormalizedAnalysis> {
        let mut dose_normalized_auc = HashMap::new();
        let mut dose_normalized_cmax = HashMap::new();
        let mut dose_linearity_assessment = HashMap::new();

        // Group subjects by treatment/formulation for dose linearity assessment
        let treatment_groups = Self::group_by_treatment(subjects);

        for (treatment, treatment_subjects) in treatment_groups {
            let treatment_results: Vec<&NcaResults> = results
                .iter()
                .filter(|r| treatment_subjects.iter().any(|s| s.id == r.subject_id))
                .collect();

            if treatment_results.len() < 3 {
                continue;
            }

            // Calculate dose-normalized parameters
            let (dn_auc_values, dn_cmax_values, doses) = Self::calculate_dose_normalized_values(
                &treatment_results, &treatment_subjects
            );

            if !dn_auc_values.is_empty() {
                let auc_stats = Self::calculate_parameter_stats(&dn_auc_values);
                dose_normalized_auc.insert(treatment.clone(), auc_stats);
            }

            if !dn_cmax_values.is_empty() {
                let cmax_stats = Self::calculate_parameter_stats(&dn_cmax_values);
                dose_normalized_cmax.insert(treatment.clone(), cmax_stats);
            }

            // Assess dose linearity
            if doses.len() >= 3 {
                let linearity = Self::assess_dose_linearity(&doses, &dn_auc_values);
                dose_linearity_assessment.insert(treatment, linearity);
            }
        }

        Ok(DoseNormalizedAnalysis {
            dose_normalized_auc,
            dose_normalized_cmax,
            dose_linearity_assessment,
        })
    }

    fn group_by_treatment(subjects: &[Subject]) -> HashMap<String, Vec<Subject>> {
        let mut groups = HashMap::new();

        for subject in subjects {
            let treatment = subject.demographics.treatment
                .clone()
                .unwrap_or_else(|| "Unknown".to_string());
            
            groups.entry(treatment).or_insert_with(Vec::new).push(subject.clone());
        }

        groups
    }

    fn calculate_dose_normalized_values(
        results: &[&NcaResults],
        subjects: &[Subject],
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let mut dn_auc_values = Vec::new();
        let mut dn_cmax_values = Vec::new();
        let mut doses = Vec::new();

        for result in results {
            if let Some(subject) = subjects.iter().find(|s| s.id == result.subject_id) {
                let total_dose: f64 = subject.dosing_events.iter().map(|d| d.dose).sum();
                
                if total_dose > 0.0 {
                    if let Some(auc) = result.individual_parameters.auc_inf {
                        dn_auc_values.push(auc / total_dose);
                        doses.push(total_dose);
                    }
                    
                    if let Some(cmax) = result.individual_parameters.cmax {
                        dn_cmax_values.push(cmax / total_dose);
                    }
                }
            }
        }

        (dn_auc_values, dn_cmax_values, doses)
    }

    fn assess_dose_linearity(doses: &[f64], dn_auc_values: &[f64]) -> LinearityAssessment {
        if doses.len() != dn_auc_values.len() || doses.len() < 3 {
            return LinearityAssessment {
                slope: 0.0,
                r_squared: 0.0,
                linearity_conclusion: "Insufficient data".to_string(),
            };
        }

        // Linear regression of dose-normalized AUC vs dose
        // If linear, slope should be close to 0
        let mean_dose = doses.mean();
        let mean_dn_auc = dn_auc_values.mean();

        let numerator: f64 = doses.iter().zip(dn_auc_values.iter())
            .map(|(d, auc)| (d - mean_dose) * (auc - mean_dn_auc))
            .sum();

        let denominator: f64 = doses.iter()
            .map(|d| (d - mean_dose).powi(2))
            .sum();

        let slope = if denominator != 0.0 { numerator / denominator } else { 0.0 };

        // Calculate R-squared
        let ss_tot: f64 = dn_auc_values.iter()
            .map(|auc| (auc - mean_dn_auc).powi(2))
            .sum();

        let ss_res: f64 = doses.iter().zip(dn_auc_values.iter())
            .map(|(d, auc)| {
                let predicted = mean_dn_auc + slope * (d - mean_dose);
                (auc - predicted).powi(2)
            })
            .sum();

        let r_squared = if ss_tot != 0.0 { 1.0 - (ss_res / ss_tot) } else { 0.0 };

        // Assess linearity
        let linearity_conclusion = if slope.abs() < 0.1 && r_squared < 0.3 {
            "Linear pharmacokinetics".to_string()
        } else if slope.abs() > 0.3 || r_squared > 0.7 {
            "Non-linear pharmacokinetics".to_string()
        } else {
            "Inconclusive".to_string()
        };

        LinearityAssessment {
            slope,
            r_squared,
            linearity_conclusion,
        }
    }

    fn calculate_parameter_stats(values: &[f64]) -> ParameterStats {
        let n = values.len();
        
        if n == 0 {
            return ParameterStats {
                n: 0,
                mean: 0.0, 
                arithmetic_mean: 0.0,
                arithmetic_std: 0.0,
                arithmetic_cv_percent: 0.0,
                std: 0.0, 
                cv_percent: 0.0, 
                median: 0.0,
                q25: 0.0,
                q75: 0.0,
                min: 0.0,
                max: 0.0,
                geometric_mean: None,
                geometric_cv_percent: None,
            };
        }

        let arithmetic_mean = values.mean();
        let arithmetic_std = values.std_dev();
        let arithmetic_cv_percent = if arithmetic_mean != 0.0 { (arithmetic_std / arithmetic_mean) * 100.0 } else { 0.0 };

        let mut sorted_values = values.to_vec();
        sorted_values.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let median = if n % 2 == 0 {
            (sorted_values[n / 2 - 1] + sorted_values[n / 2]) / 2.0
        } else {
            sorted_values[n / 2]
        };

        let q25_idx = ((n as f64 * 0.25) as usize).min(n - 1);
        let q75_idx = ((n as f64 * 0.75) as usize).min(n - 1);
        let q25 = sorted_values[q25_idx];
        let q75 = sorted_values[q75_idx];

        let min = sorted_values[0];
        let max = sorted_values[n - 1];

        // Geometric statistics
        let (geometric_mean, geometric_cv_percent) = if values.iter().all(|&v| v > 0.0) {
            let ln_values: Vec<f64> = values.iter().map(|v| v.ln()).collect();
            let ln_mean = ln_values.as_slice().mean();
            let ln_std = ln_values.as_slice().std_dev();
            let geo_mean = ln_mean.exp();
            let geo_cv = ((ln_std.exp().powi(2) - 1.0).sqrt()) * 100.0;
            (Some(geo_mean), Some(geo_cv))
        } else {
            (None, None)
        };

        ParameterStats {
            n,
            mean: arithmetic_mean, 
            arithmetic_mean,
            arithmetic_std,
            arithmetic_cv_percent,
            std: arithmetic_std, 
            cv_percent: arithmetic_cv_percent, 
            median,
            q25,
            q75,
            min,
            max,
            geometric_mean,
            geometric_cv_percent,
        }
    }
}