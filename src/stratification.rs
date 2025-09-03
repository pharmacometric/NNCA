use crate::{models::*, population::PopulationAnalyzer, Result};
use std::collections::HashMap;
use rayon::prelude::*;
use statrs::statistics::Statistics;
use serde::{Serialize, Deserialize};

pub struct StratificationAnalyzer;

impl StratificationAnalyzer {
    /// Perform stratified analysis based on specified variables
    pub fn analyze_stratified(
        subjects: &[Subject],
        config: &AnalysisConfig,
    ) -> Result<HashMap<String, StratifiedResults>> {
        let stratification_config = match &config.stratification {
            Some(strat_config) => strat_config,
            None => return Ok(HashMap::new()),
        };

        let mut stratified_results = HashMap::new();

        // Single variable stratification
        for variable in &stratification_config.stratify_columns {
            let strata = Self::create_strata(subjects, variable);
            
            for (stratum_value, stratum_subjects) in strata {
                if stratum_subjects.len() < stratification_config.minimum_n_per_stratum {
                    log::warn!(
                        "Skipping stratum {}={} (n={}, minimum required: {})",
                        variable, stratum_value, stratum_subjects.len(), 
                        stratification_config.minimum_n_per_stratum
                    );
                    continue;
                }

                let stratum_key = format!("{}_{}", variable, stratum_value);
                let stratum_results = Self::analyze_stratum(
                    &stratum_subjects, 
                    config, 
                    variable, 
                    &stratum_value
                )?;
                
                stratified_results.insert(stratum_key, stratum_results);
            }
        }

        // Interaction analysis if requested
        if stratification_config.include_interactions && stratification_config.stratify_columns.len() >= 2 {
            let interaction_results = Self::analyze_interactions(subjects, config)?;
            stratified_results.extend(interaction_results);
        }

        Ok(stratified_results)
    }

    fn create_strata(subjects: &[Subject], variable: &str) -> HashMap<String, Vec<Subject>> {
        let mut strata = HashMap::new();

        for subject in subjects {
            let stratum_value = Self::get_stratum_value(subject, variable);
            
            if let Some(value) = stratum_value {
                strata.entry(value).or_insert_with(Vec::new).push(subject.clone());
            }
        }

        strata
    }

    fn get_stratum_value(subject: &Subject, variable: &str) -> Option<String> {
        match variable.to_uppercase().as_str() {
            "SEX" => subject.demographics.sex.clone(),
            "RACE" => subject.demographics.race.clone(),
            "TREATMENT" | "TRT" => subject.demographics.treatment.clone(),
            "PERIOD" => subject.demographics.period.map(|p| p.to_string()),
            "SEQUENCE" | "SEQ" => subject.demographics.sequence.clone(),
            "FORMULATION" | "FORM" => subject.demographics.formulation.clone(),
            "AGE_GROUP" => Self::categorize_age(subject.demographics.age),
            "WEIGHT_GROUP" => Self::categorize_weight(subject.demographics.weight),
            "DOSE_GROUP" => Self::categorize_dose(subject),
            _ => None,
        }
    }

    fn categorize_age(age: Option<f64>) -> Option<String> {
        age.map(|a| {
            if a < 18.0 { "Pediatric".to_string() }
            else if a < 65.0 { "Adult".to_string() }
            else { "Elderly".to_string() }
        })
    }

    fn categorize_weight(weight: Option<f64>) -> Option<String> {
        weight.map(|w| {
            if w < 60.0 { "Low".to_string() }
            else if w < 90.0 { "Normal".to_string() }
            else { "High".to_string() }
        })
    }

    fn categorize_dose(subject: &Subject) -> Option<String> {
        let total_dose: f64 = subject.dosing_events.iter().map(|d| d.dose).sum();
        
        if total_dose > 0.0 {
            if total_dose < 100.0 { Some("Low".to_string()) }
            else if total_dose < 500.0 { Some("Medium".to_string()) }
            else { Some("High".to_string()) }
        } else {
            None
        }
    }

    fn analyze_stratum(
        subjects: &[Subject],
        config: &AnalysisConfig,
        variable: &str,
        value: &str,
    ) -> Result<StratifiedResults> {
        log::info!("Analyzing stratum: {} = {} (n = {})", variable, value, subjects.len());

        // Perform population analysis for this stratum
        let population_results = PopulationAnalyzer::analyze_population(subjects.to_vec(), config)?;

        Ok(StratifiedResults {
            stratum_name: variable.to_string(),
            stratum_value: value.to_string(),
            n_subjects: subjects.len(),
            individual_results: population_results.individual_results,
            summary_statistics: population_results.summary_statistics,
            method_comparison: population_results.method_comparison,
        })
    }

    fn analyze_interactions(
        subjects: &[Subject],
        config: &AnalysisConfig,
    ) -> Result<HashMap<String, StratifiedResults>> {
        let mut interaction_results = HashMap::new();
        let variables = &config.stratification.as_ref().unwrap().stratify_columns;

        // Two-way interactions
        for i in 0..variables.len() {
            for j in (i + 1)..variables.len() {
                let var1 = &variables[i];
                let var2 = &variables[j];
                
                let interaction_strata = Self::create_interaction_strata(subjects, var1, var2);
                
                for (interaction_key, stratum_subjects) in interaction_strata {
                    if stratum_subjects.len() >= config.stratification.as_ref().unwrap().minimum_n_per_stratum {
                        let stratum_results = Self::analyze_stratum(
                            &stratum_subjects,
                            config,
                            &format!("{}_{}", var1, var2),
                            &interaction_key,
                        )?;
                        
                        interaction_results.insert(
                            format!("interaction_{}_{}_{}", var1, var2, interaction_key),
                            stratum_results,
                        );
                    }
                }
            }
        }

        Ok(interaction_results)
    }

    fn create_interaction_strata(
        subjects: &[Subject],
        var1: &str,
        var2: &str,
    ) -> HashMap<String, Vec<Subject>> {
        let mut strata = HashMap::new();

        for subject in subjects {
            let value1 = Self::get_stratum_value(subject, var1);
            let value2 = Self::get_stratum_value(subject, var2);
            
            if let (Some(v1), Some(v2)) = (value1, value2) {
                let interaction_key = format!("{}-{}", v1, v2);
                strata.entry(interaction_key).or_insert_with(Vec::new).push(subject.clone());
            }
        }

        strata
    }

    /// Perform statistical comparison between strata
    pub fn compare_strata(
        strata_results: &HashMap<String, StratifiedResults>,
        parameter: &str,
    ) -> Result<StrataComparison> {
        let mut comparisons = Vec::new();

        let strata_names: Vec<&String> = strata_results.keys().collect();
        
        for i in 0..strata_names.len() {
            for j in (i + 1)..strata_names.len() {
                let stratum1 = &strata_results[strata_names[i]];
                let stratum2 = &strata_results[strata_names[j]];
                
                let comparison = Self::perform_statistical_test(stratum1, stratum2, parameter)?;
                comparisons.push(comparison);
            }
        }

        Ok(StrataComparison {
            parameter: parameter.to_string(),
            pairwise_comparisons: comparisons,
        })
    }

    fn perform_statistical_test(
        stratum1: &StratifiedResults,
        stratum2: &StratifiedResults,
        parameter: &str,
    ) -> Result<PairwiseComparison> {
        let values1 = Self::extract_parameter_values(&stratum1.individual_results, parameter);
        let values2 = Self::extract_parameter_values(&stratum2.individual_results, parameter);

        if values1.is_empty() || values2.is_empty() {
            return Ok(PairwiseComparison {
                stratum1_name: stratum1.stratum_value.clone(),
                stratum2_name: stratum2.stratum_value.clone(),
                n1: values1.len(),
                n2: values2.len(),
                mean1: 0.0,
                mean2: 0.0,
                p_value: 1.0,
                test_statistic: 0.0,
                test_type: "insufficient_data".to_string(),
                significant: false,
                effect_size: 0.0,
            });
        }

        let mean1 = values1.as_slice().mean();
        let mean2 = values2.as_slice().mean();
        
        // Perform Welch's t-test (unequal variances)
        let (t_stat, p_value) = Self::welch_t_test(&values1, &values2);
        
        // Calculate effect size (Cohen's d)
        let pooled_std = Self::calculate_pooled_std(&values1, &values2);
        let effect_size = if pooled_std > 0.0 { (mean1 - mean2).abs() / pooled_std } else { 0.0 };

        Ok(PairwiseComparison {
            stratum1_name: stratum1.stratum_value.clone(),
            stratum2_name: stratum2.stratum_value.clone(),
            n1: values1.len(),
            n2: values2.len(),
            mean1,
            mean2,
            p_value,
            test_statistic: t_stat,
            test_type: "welch_t_test".to_string(),
            significant: p_value < 0.05,
            effect_size,
        })
    }

    fn extract_parameter_values(results: &[NcaResults], parameter: &str) -> Vec<f64> {
        results
            .iter()
            .filter_map(|r| {
                let params = &r.individual_parameters;
                match parameter {
                    "auc_last" => params.auc_last,
                    "auc_inf" => params.auc_inf,
                    "cmax" => params.cmax,
                    "tmax" => params.tmax,
                    "half_life" => params.half_life,
                    "clearance" => params.clearance,
                    "volume_terminal" => params.volume_terminal,
                    "mrt" => params.mrt,
                    _ => None,
                }
            })
            .collect()
    }

    fn welch_t_test(values1: &[f64], values2: &[f64]) -> (f64, f64) {
        if values1.len() < 2 || values2.len() < 2 {
            return (0.0, 1.0);
        }

        let mean1 = values1.mean();
        let mean2 = values2.mean();
        let var1 = values1.variance();
        let var2 = values2.variance();
        let n1 = values1.len() as f64;
        let n2 = values2.len() as f64;

        let se = ((var1 / n1) + (var2 / n2)).sqrt();
        let t_stat = if se > 0.0 { (mean1 - mean2) / se } else { 0.0 };

        // Welch-Satterthwaite degrees of freedom
        let df = if var1 > 0.0 && var2 > 0.0 {
            let numerator = ((var1 / n1) + (var2 / n2)).powi(2);
            let denominator = (var1 / n1).powi(2) / (n1 - 1.0) + (var2 / n2).powi(2) / (n2 - 1.0);
            numerator / denominator
        } else {
            n1 + n2 - 2.0
        };

        // Approximate p-value using t-distribution
        let p_value = Self::t_distribution_p_value(t_stat.abs(), df);

        (t_stat, p_value)
    }

    fn calculate_pooled_std(values1: &[f64], values2: &[f64]) -> f64 {
        if values1.len() < 2 || values2.len() < 2 {
            return 0.0;
        }

        let n1 = values1.len() as f64;
        let n2 = values2.len() as f64;
        let var1 = values1.variance();
        let var2 = values2.variance();

        let pooled_variance = ((n1 - 1.0) * var1 + (n2 - 1.0) * var2) / (n1 + n2 - 2.0);
        pooled_variance.sqrt()
    }

    fn t_distribution_p_value(t_abs: f64, df: f64) -> f64 {
        // Simplified approximation for p-value calculation
        // In production, use a proper statistical library
        if df <= 0.0 { return 1.0; }
        
        // Very rough approximation - replace with proper implementation
        let z_approx = t_abs * (1.0 - 1.0 / (4.0 * df));
        2.0 * (1.0 - Self::standard_normal_cdf(z_approx))
    }

    fn standard_normal_cdf(z: f64) -> f64 {
        // Approximation of standard normal CDF
        0.5 * (1.0 + Self::erf(z / 2.0_f64.sqrt()))
    }

    fn erf(x: f64) -> f64 {
        // Approximation of error function
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
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StrataComparison {
    pub parameter: String,
    pub pairwise_comparisons: Vec<PairwiseComparison>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PairwiseComparison {
    pub stratum1_name: String,
    pub stratum2_name: String,
    pub n1: usize,
    pub n2: usize,
    pub mean1: f64,
    pub mean2: f64,
    pub p_value: f64,
    pub test_statistic: f64,
    pub test_type: String,
    pub significant: bool,
    pub effect_size: f64,
}