use crate::{models::*, nca::NcaAnalyzer, Result};
use crate::stratification::StratificationAnalyzer;
use crate::covariate::CovariateAnalyzer;
use rayon::prelude::*;
use statrs::statistics::Statistics;
use std::collections::HashMap;

pub struct PopulationAnalyzer;

impl PopulationAnalyzer {
    /// Perform population NCA analysis with parallel processing
    pub fn analyze_population(
        subjects: Vec<Subject>,
        config: &AnalysisConfig,
    ) -> Result<PopulationResults> {
        log::info!("Starting population analysis for {} subjects", subjects.len());

        // Parallel processing of individual subjects
        let mut individual_results = Vec::new();
        let mut failed_subjects = Vec::new();
        
        let analysis_results: Vec<_> = subjects
            .par_iter()
            .map(|subject| {
                match NcaAnalyzer::analyze_subject(subject, config) {
                    Ok((result, warnings)) => {
                        let validation_warnings = NcaAnalyzer::validate_results(&result);
                        let all_warnings = [warnings, validation_warnings].concat();
                        
                        if !all_warnings.is_empty() {
                            log::warn!("Warnings for subject {}: {:?}", subject.id, all_warnings);
                        }
                        Ok((result, all_warnings))
                    }
                    Err(e) => {
                        log::error!("Failed to analyze subject {}: {}", subject.id, e);
                        
                        // Count quantifiable concentrations for failed subject
                        let quantifiable_count = subject.observations.iter()
                            .filter(|obs| obs.concentration > 0.0 && !obs.bloq)
                            .count();
                        
                        let failed_analysis = FailedSubjectAnalysis {
                            subject_id: subject.id.clone(),
                            failure_reason: e.to_string(),
                            quantifiable_concentrations: quantifiable_count,
                            total_observations: subject.observations.len(),
                            failed_parameters: vec!["All parameters".to_string()],
                        };
                        
                        Err(failed_analysis)
                    }
                }
            })
            .collect();
        
        // Separate successful and failed analyses
        for result in analysis_results {
            match result {
                Ok((nca_result, _warnings)) => {
                    individual_results.push(nca_result);
                }
                Err(failed_analysis) => {
                    failed_subjects.push(failed_analysis);
                }
            }
        }

        log::info!("Successfully analyzed {} subjects", individual_results.len());
        if !failed_subjects.is_empty() {
            log::warn!("Failed to analyze {} subjects", failed_subjects.len());
        }

        // Calculate summary statistics
        let summary_statistics = Self::calculate_summary_statistics(&individual_results)?;

        // Method comparison across all subjects
        let method_comparison = Self::perform_method_comparison(&individual_results)?;

        // Stratified analysis
        let stratified_results = StratificationAnalyzer::analyze_stratified(&subjects, config)?;

        // Covariate analysis
        let covariate_analysis = if config.perform_covariate_analysis {
            CovariateAnalyzer::analyze_covariates(&individual_results, &subjects)?
        } else {
            CovariateAnalysis {
                correlations: HashMap::new(),
                regression_analysis: HashMap::new(),
                dose_normalized_analysis: None,
            }
        };

        Ok(PopulationResults {
            individual_results,
            failed_subjects,
            summary_statistics,
            method_comparison,
            stratified_results,
            covariate_analysis,
        })
    }

    fn calculate_summary_statistics(results: &[NcaResults]) -> Result<SummaryStatistics> {
        let mut parameter_stats = HashMap::new();

        // Define parameters to analyze
        let parameters: Vec<(&str, fn(&IndividualParameters) -> Option<f64>)> = vec![
            ("auc_last", |p| p.auc_last),
            ("auc_inf", |p| p.auc_inf),
            ("cmax", |p| p.cmax),
            ("tmax", |p| p.tmax),
            ("half_life", |p| p.half_life),
            ("clearance", |p| p.clearance),
            ("volume_terminal", |p| p.volume_terminal),
            ("mrt", |p| p.mrt),
        ];

        for (param_name, extractor) in parameters {
            let values: Vec<f64> = results
                .iter()
                .filter_map(|r| extractor(&r.individual_parameters))
                .collect();

            if !values.is_empty() {
                let stats = Self::calculate_parameter_stats(&values);
                parameter_stats.insert(param_name.to_string(), stats);
            }
        }

        Ok(SummaryStatistics { parameter_stats })
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

        let mean = values.mean();
        let std = values.std_dev();
        let cv_percent = if mean != 0.0 { (std / mean) * 100.0 } else { 0.0 };

        let mut sorted_values = values.to_vec();
        sorted_values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

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

        // Geometric statistics (for positive values only)
        let (geometric_mean, geometric_cv_percent) = if values.iter().all(|&v| v > 0.0) {
            let ln_values: Vec<f64> = values.iter().map(|v| v.ln()).collect();
            let ln_mean = (&ln_values).mean();
            let ln_std = (&ln_values).std_dev();
            let geo_mean = ln_mean.exp();
            let geo_cv = ((ln_std.exp().powi(2) - 1.0).sqrt()) * 100.0;
            (Some(geo_mean), Some(geo_cv))
        } else {
            (None, None)
        };

        ParameterStats {
            n,
            mean,
            arithmetic_mean: mean,
            arithmetic_std: std,
            arithmetic_cv_percent: cv_percent,
            std,
            cv_percent,
            median,
            q25,
            q75,
            min,
            max,
            geometric_mean,
            geometric_cv_percent,
        }
    }

    fn perform_method_comparison(results: &[NcaResults]) -> Result<MethodComparison> {
        let mut auc_methods = HashMap::new();
        let correlation_matrix = HashMap::new();
        let bias_analysis = HashMap::new();

        // Collect AUC values by method
        let mut method_values: HashMap<String, Vec<f64>> = HashMap::new();

        for result in results {
            for (method, params) in &result.method_comparisons {
                if let Some(auc) = params.auc_last {
                    method_values.entry(method.clone()).or_insert_with(Vec::new).push(auc);
                }
            }
        }

        // Calculate mean AUC by method
        for (method, values) in method_values.iter() {
            if !values.is_empty() {
                auc_methods.insert(method.clone(), values.mean());
            }
        }

        // For now, return simplified method comparison
        // Full correlation and bias analysis would require additional implementation
        Ok(MethodComparison {
            auc_methods,
            correlation_matrix,
            bias_analysis,
        })
    }

    #[allow(dead_code)]
    fn calculate_correlation(values1: &[f64], values2: &[f64]) -> f64 {
        if values1.len() != values2.len() || values1.len() < 2 {
            return 0.0;
        }

        let mean1 = values1.mean();
        let mean2 = values2.mean();

        let numerator: f64 = values1
            .iter()
            .zip(values2.iter())
            .map(|(x, y)| (x - mean1) * (y - mean2))
            .sum();

        let sum_sq1: f64 = values1.iter().map(|x| (x - mean1).powi(2)).sum();
        let sum_sq2: f64 = values2.iter().map(|y| (y - mean2).powi(2)).sum();

        let denominator = (sum_sq1 * sum_sq2).sqrt();

        if denominator == 0.0 {
            0.0
        } else {
            numerator / denominator
        }
    }

    #[allow(dead_code)]
    fn calculate_bias_analysis(values1: &[f64], values2: &[f64]) -> BiasAnalysis {
        if values1.len() != values2.len() || values1.is_empty() {
            return BiasAnalysis {
                mean_difference: 0.0,
                mean_percent_difference: 0.0,
                limits_of_agreement: (0.0, 0.0),
            };
        }

        let differences: Vec<f64> = values1
            .iter()
            .zip(values2.iter())
            .map(|(x, y)| x - y)
            .collect();

        let percent_differences: Vec<f64> = values1
            .iter()
            .zip(values2.iter())
            .map(|(x, y)| {
                let mean_val = (x + y) / 2.0;
                if mean_val != 0.0 {
                    ((x - y) / mean_val) * 100.0
                } else {
                    0.0
                }
            })
            .collect();

        let mean_diff = (&differences).mean();
        let std_diff = (&differences).std_dev();
        let mean_percent_diff = (&percent_differences).mean();

        // 95% limits of agreement (mean ± 1.96 * SD)
        let lower_limit = mean_diff - 1.96 * std_diff;
        let upper_limit = mean_diff + 1.96 * std_diff;

        BiasAnalysis {
            mean_difference: mean_diff,
            mean_percent_difference: mean_percent_diff,
            limits_of_agreement: (lower_limit, upper_limit),
        }
    }
}