use crate::{models::*, Result};
use serde_json;
use std::fs::{self, File};
use std::io::Write;
use std::path::Path;
use std::collections::HashMap;

pub struct OutputManager;

impl OutputManager {
    pub fn save_results<P: AsRef<Path>>(
        results: &PopulationResults,
        config: &AnalysisConfig,
        output_path: P,
    ) -> Result<()> {
        let output_dir = output_path.as_ref();
        fs::create_dir_all(output_dir)?;

        // Save individual results
        Self::save_individual_results(&results.individual_results, output_dir)?;
        
        // Save failed subjects log
        Self::save_failed_subjects_log(&results.failed_subjects, output_dir)?;
        
        // Save summary statistics
        Self::save_summary_statistics(&results.summary_statistics, output_dir)?;
        
        // Save method comparison
        Self::save_method_comparison(&results.method_comparison, output_dir)?;
        
        // Save stratified results
        Self::save_stratified_results(&results.stratified_results, output_dir)?;
        
        // Save covariate analysis
        Self::save_covariate_analysis(&results.covariate_analysis, output_dir)?;
        
        // Save complete results as JSON
        Self::save_json_results(results, output_dir)?;
        
        // Save CSV summary
        Self::save_csv_summary(results, output_dir)?;
        
        // Generate analysis report
        Self::generate_analysis_report(results, config, output_dir)?;

        log::info!("Results saved to: {}", output_dir.display());
        Ok(())
    }

    fn save_individual_results(
        results: &[NcaResults],
        output_dir: &Path,
    ) -> Result<()> {
        let file_path = output_dir.join("individual_results.csv");
        let mut file = File::create(file_path)?;
        
        // Write header
        writeln!(file, "SUBJECT_ID,AUC_LAST,AUC_INF,AUC_INF_PRED,AUC_EXTRAP_PERCENT,AUMC_LAST,AUMC_INF,CMAX,TMAX,TLAST,CLAST,HALF_LIFE,LAMBDA_Z,LAMBDA_Z_R2,CLEARANCE,VSS,VZ,MRT")?;
        
        // Write data
        for result in results {
            let p = &result.individual_parameters;
            writeln!(
                file,
                "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                result.subject_id,
                p.auc_last.map_or("NA".to_string(), |v| v.to_string()),
                p.auc_inf.map_or("NA".to_string(), |v| v.to_string()),
                p.auc_inf_pred.map_or("NA".to_string(), |v| v.to_string()),
                p.auc_percent_extrap.map_or("NA".to_string(), |v| v.to_string()),
                p.aumc_last.map_or("NA".to_string(), |v| v.to_string()),
                p.aumc_inf.map_or("NA".to_string(), |v| v.to_string()),
                p.cmax.map_or("NA".to_string(), |v| v.to_string()),
                p.tmax.map_or("NA".to_string(), |v| v.to_string()),
                p.tlast.map_or("NA".to_string(), |v| v.to_string()),
                p.clast.map_or("NA".to_string(), |v| v.to_string()),
                p.half_life.map_or("NA".to_string(), |v| v.to_string()),
                p.lambda_z.map_or("NA".to_string(), |v| v.to_string()),
                p.lambda_z_r_squared.map_or("NA".to_string(), |v| v.to_string()),
                p.clearance.map_or("NA".to_string(), |v| v.to_string()),
                p.volume_steady_state.map_or("NA".to_string(), |v| v.to_string()),
                p.volume_terminal.map_or("NA".to_string(), |v| v.to_string()),
                p.mrt.map_or("NA".to_string(), |v| v.to_string()),
            )?;
        }
        
        Ok(())
    }

    fn save_summary_statistics(
        summary: &SummaryStatistics,
        output_dir: &Path,
    ) -> Result<()> {
        let file_path = output_dir.join("summary_statistics.csv");
        let mut file = File::create(file_path)?;
        
        writeln!(file, "PARAMETER,N,MEAN,STD,CV_PERCENT,MEDIAN,Q25,Q75,MIN,MAX,GEO_MEAN,GEO_CV_PERCENT")?;
        
        for (param, stats) in &summary.parameter_stats {
            writeln!(
                file,
                "{},{},{:.6},{:.6},{:.2},{:.6},{:.6},{:.6},{:.6},{:.6},{},{}",
                param,
                stats.n,
                stats.arithmetic_mean,
                stats.arithmetic_std,
                stats.arithmetic_cv_percent,
                stats.median,
                stats.q25,
                stats.q75,
                stats.min,
                stats.max,
                stats.geometric_mean.map_or("NA".to_string(), |v| format!("{:.6}", v)),
                stats.geometric_cv_percent.map_or("NA".to_string(), |v| format!("{:.2}", v)),
            )?;
        }
        
        Ok(())
    }

    fn save_failed_subjects_log(
        failed_subjects: &[FailedSubjectAnalysis],
        output_dir: &Path,
    ) -> Result<()> {
        if failed_subjects.is_empty() {
            return Ok(());
        }

        let file_path = output_dir.join("failed_subjects.log");
        let mut file = File::create(file_path)?;
        
        writeln!(file, "FAILED SUBJECT ANALYSIS LOG")?;
        writeln!(file, "==========================")?;
        writeln!(file)?;
        writeln!(file, "Total failed subjects: {}", failed_subjects.len())?;
        writeln!(file)?;
        
        for failed in failed_subjects {
            writeln!(file, "Subject ID: {}", failed.subject_id)?;
            writeln!(file, "Failure Reason: {}", failed.failure_reason)?;
            writeln!(file, "Quantifiable Concentrations: {}", failed.quantifiable_concentrations)?;
            writeln!(file, "Total Observations: {}", failed.total_observations)?;
            writeln!(file, "Failed Parameters: {}", failed.failed_parameters.join(", "))?;
            writeln!(file, "---")?;
        }
        
        Ok(())
    }

    fn save_method_comparison(
        comparison: &MethodComparison,
        output_dir: &Path,
    ) -> Result<()> {
        // Save method means
        let file_path = output_dir.join("method_comparison.csv");
        let mut file = File::create(file_path)?;
        
        writeln!(file, "METHOD,MEAN_AUC")?;
        for (method, mean_auc) in &comparison.auc_methods {
            writeln!(file, "{},{:.6}", method, mean_auc)?;
        }
        
        // Save correlation matrix
        let corr_path = output_dir.join("method_correlations.csv");
        let mut corr_file = File::create(corr_path)?;
        
        let methods: Vec<&String> = comparison.correlation_matrix.keys().collect();
        write!(corr_file, "METHOD")?;
        for method in &methods {
            write!(corr_file, ",{}", method)?;
        }
        writeln!(corr_file)?;
        
        for method1 in &methods {
            write!(corr_file, "{}", method1)?;
            for method2 in &methods {
                let corr = comparison.correlation_matrix
                    .get(*method1)
                    .and_then(|m| m.get(*method2))
                    .unwrap_or(&0.0);
                write!(corr_file, ",{:.4}", corr)?;
            }
            writeln!(corr_file)?;
        }
        
        Ok(())
    }

    fn save_stratified_results(
        stratified_results: &HashMap<String, StratifiedResults>,
        output_dir: &Path,
    ) -> Result<()> {
        if stratified_results.is_empty() {
            return Ok(());
        }

        let file_path = output_dir.join("stratified_analysis.csv");
        let mut file = File::create(file_path)?;
        
        writeln!(file, "STRATUM,STRATUM_VALUE,N,PARAMETER,MEAN,STD,CV_PERCENT,MEDIAN,GEO_MEAN,GEO_CV_PERCENT")?;
        
        for (stratum_key, stratum_results) in stratified_results {
            for (param, stats) in &stratum_results.summary_statistics.parameter_stats {
                writeln!(
                    file,
                    "{},{},{},{},{:.6},{:.6},{:.2},{:.6},{},{}",
                    stratum_results.stratum_name,
                    stratum_results.stratum_value,
                    stratum_results.n_subjects,
                    param,
                    stats.mean,
                    stats.std,
                    stats.cv_percent,
                    stats.median,
                    stats.geometric_mean.map_or("NA".to_string(), |v| format!("{:.6}", v)),
                    stats.geometric_cv_percent.map_or("NA".to_string(), |v| format!("{:.2}", v)),
                )?;
            }
        }
        
        // Save detailed stratified results
        for (stratum_key, stratum_results) in stratified_results {
            let stratum_file_path = output_dir.join(format!("stratum_{}.csv", stratum_key));
            let mut stratum_file = File::create(stratum_file_path)?;
            
            writeln!(stratum_file, "SUBJECT_ID,AUC_LAST,AUC_INF,CMAX,TMAX,HALF_LIFE,CLEARANCE,VSS,VZ,MRT")?;
            
            for result in &stratum_results.individual_results {
                let p = &result.individual_parameters;
                writeln!(
                    stratum_file,
                    "{},{},{},{},{},{},{},{},{},{}",
                    result.subject_id,
                    p.auc_last.map_or("NA".to_string(), |v| v.to_string()),
                    p.auc_inf.map_or("NA".to_string(), |v| v.to_string()),
                    p.cmax.map_or("NA".to_string(), |v| v.to_string()),
                    p.tmax.map_or("NA".to_string(), |v| v.to_string()),
                    p.half_life.map_or("NA".to_string(), |v| v.to_string()),
                    p.clearance.map_or("NA".to_string(), |v| v.to_string()),
                    p.volume_steady_state.map_or("NA".to_string(), |v| v.to_string()),
                    p.volume_terminal.map_or("NA".to_string(), |v| v.to_string()),
                    p.mrt.map_or("NA".to_string(), |v| v.to_string()),
                )?;
            }
        }
        
        Ok(())
    }

    fn save_covariate_analysis(
        covariate_analysis: &CovariateAnalysis,
        output_dir: &Path,
    ) -> Result<()> {
        // Save correlations
        let corr_path = output_dir.join("covariate_correlations.csv");
        let mut corr_file = File::create(corr_path)?;
        
        writeln!(corr_file, "COVARIATE,PARAMETER,CORRELATION,P_VALUE,SIGNIFICANCE")?;
        
        for (covariate, correlation_data) in &covariate_analysis.correlations {
            for (parameter, &corr_value) in &correlation_data.parameter_correlations {
                let p_value = correlation_data.p_values.get(parameter).copied().unwrap_or(1.0);
                let significant = if p_value < 0.05 { "Yes" } else { "No" };
                
                writeln!(
                    corr_file,
                    "{},{},{:.4},{:.4},{}",
                    covariate, parameter, corr_value, p_value, significant
                )?;
            }
        }
        
        // Save regression analysis
        let reg_path = output_dir.join("regression_analysis.csv");
        let mut reg_file = File::create(reg_path)?;
        
        writeln!(reg_file, "PARAMETER,COVARIATE,SLOPE,INTERCEPT,R_SQUARED,P_VALUE,CI_LOWER,CI_UPPER")?;
        
        for (key, regression) in &covariate_analysis.regression_analysis {
            writeln!(
                reg_file,
                "{},{},{:.6},{:.6},{:.4},{:.4},{:.6},{:.6}",
                regression.parameter,
                regression.covariate,
                regression.slope,
                regression.intercept,
                regression.r_squared,
                regression.p_value,
                regression.confidence_interval.0,
                regression.confidence_interval.1,
            )?;
        }
        
        // Save dose normalization analysis
        if let Some(dose_analysis) = &covariate_analysis.dose_normalized_analysis {
            let dose_path = output_dir.join("dose_normalized_analysis.csv");
            let mut dose_file = File::create(dose_path)?;
            
            writeln!(dose_file, "TREATMENT,PARAMETER,N,MEAN,STD,CV_PERCENT,LINEARITY_ASSESSMENT")?;
            
            for (treatment, stats) in &dose_analysis.dose_normalized_auc {
                let linearity = dose_analysis.dose_linearity_assessment
                    .get(treatment)
                    .map(|l| l.linearity_conclusion.clone())
                    .unwrap_or_else(|| "Unknown".to_string());
                
                writeln!(
                    dose_file,
                    "{},AUC_DN,{},{:.6},{:.6},{:.2},{}",
                    treatment, stats.n, stats.mean, stats.std, stats.cv_percent, linearity
                )?;
            }
            
            for (treatment, stats) in &dose_analysis.dose_normalized_cmax {
                writeln!(
                    dose_file,
                    "{},CMAX_DN,{},{:.6},{:.6},{:.2},NA",
                    treatment, stats.n, stats.mean, stats.std, stats.cv_percent
                )?;
            }
        }
        
        Ok(())
    }

    fn save_json_results(
        results: &PopulationResults,
        output_dir: &Path,
    ) -> Result<()> {
        let file_path = output_dir.join("complete_results.json");
        let json_string = serde_json::to_string_pretty(results)?;
        fs::write(file_path, json_string)?;
        Ok(())
    }

    fn save_csv_summary(
        results: &PopulationResults,
        output_dir: &Path,
    ) -> Result<()> {
        let file_path = output_dir.join("population_summary.csv");
        let mut file = File::create(file_path)?;
        
        writeln!(file, "ANALYSIS_SUMMARY")?;
        writeln!(file, "Total Subjects,{}", results.individual_results.len())?;
        writeln!(file, "Successful Analyses,{}", results.individual_results.len())?;
        writeln!(file)?;
        
        writeln!(file, "PARAMETER,N,MEAN,MEDIAN,CV%,GEO_MEAN,GEO_CV%")?;
        for (param, stats) in &results.summary_statistics.parameter_stats {
            writeln!(
                file,
                "{},{},{:.3},{:.3},{:.1},{},{}",
                param,
                stats.n,
                stats.arithmetic_mean,
                stats.median,
                stats.arithmetic_cv_percent,
                stats.geometric_mean.map_or("NA".to_string(), |v| format!("{:.3}", v)),
                stats.geometric_cv_percent.map_or("NA".to_string(), |v| format!("{:.1}", v)),
            )?;
        }
        
        Ok(())
    }

    fn generate_analysis_report(
        results: &PopulationResults,
        config: &AnalysisConfig,
        output_dir: &Path,
    ) -> Result<()> {
        let file_path = output_dir.join("analysis_report.txt");
        let mut file = File::create(file_path)?;
        
        writeln!(file, "PHARMACOKINETICS NON-COMPARTMENTAL ANALYSIS REPORT")?;
        writeln!(file, "==================================================")?;
        writeln!(file)?;
        
        writeln!(file, "Analysis Configuration:")?;
        writeln!(file, "- Time units: {}", config.time_units)?;
        writeln!(file, "- Concentration units: {}", config.concentration_units)?;
        writeln!(file, "- LLOQ handling: {:?}", config.lloq_handling)?;
        writeln!(file, "- Lambda_z selection: {:?}", config.lambda_z_selection)?;
        writeln!(file)?;
        
        writeln!(file, "Population Summary:")?;
        writeln!(file, "- Total subjects analyzed: {}", results.individual_results.len())?;
        if !results.failed_subjects.is_empty() {
            writeln!(file, "- Failed subjects: {}", results.failed_subjects.len())?;
        }
        writeln!(file)?;
        
        writeln!(file, "Key Parameters (Geometric Mean ± Geometric CV%):")?;
        for (param, stats) in &results.summary_statistics.parameter_stats {
            if let (Some(geo_mean), Some(geo_cv)) = (stats.geometric_mean, stats.geometric_cv_percent) {
                writeln!(file, "- {} (Arithmetic): {:.3} ± {:.1}%", param, stats.arithmetic_mean, stats.arithmetic_cv_percent)?;
                writeln!(file, "- {} (Geometric): {:.3} ± {:.1}%", param, geo_mean, geo_cv)?;
            }
        }
        
        writeln!(file)?;
        writeln!(file, "Method Comparison:")?;
        for (method, mean_auc) in &results.method_comparison.auc_methods {
            writeln!(file, "- {}: {:.3}", method, mean_auc)?;
        }
        
        Ok(())
    }
}