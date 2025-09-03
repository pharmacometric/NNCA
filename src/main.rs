use clap::{Arg, Command};
use env_logger;
use nca_analysis::{
    models::*,
    parser::NonmemParser,
    population::PopulationAnalyzer,
    output::OutputManager,
    example_data::ExampleDataGenerator,
    Result,
};
use std::path::PathBuf;

fn main() -> Result<()> {
    env_logger::init();

    let matches = Command::new("NCA Analysis Tool")
        .version("1.0")
        .author("Pharmacokinetics Analysis Suite")
        .about("Comprehensive non-compartmental pharmacokinetics analysis")
        .arg(
            Arg::new("input")
                .short('i')
                .long("input")
                .value_name("FILE")
                .help("Input NONMEM dataset file")
                .required_unless_present("generate-example"),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("DIR")
                .help("Output directory for results")
                .default_value("./nca_results"),
        )
        .arg(
            Arg::new("generate-example")
                .long("generate-example")
                .help("Generate example dataset with 20 subjects")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("subjects")
                .short('n')
                .long("subjects")
                .value_name("NUMBER")
                .help("Number of subjects for example dataset")
                .default_value("20"),
        )
        .arg(
            Arg::new("lloq-handling")
                .long("lloq-handling")
                .value_name("METHOD")
                .help("LLOQ handling method: zero, drop, half-lloq")
                .default_value("half-lloq"),
        )
        .arg(
            Arg::new("lambda-z-method")
                .long("lambda-z-method")
                .value_name("METHOD")
                .help("Lambda_z selection method: auto, best-fit")
                .default_value("auto"),
        )
        .arg(
            Arg::new("time-units")
                .long("time-units")
                .value_name("UNITS")
                .help("Time units")
                .default_value("h"),
        )
        .arg(
            Arg::new("conc-units")
                .long("conc-units")
                .value_name("UNITS")
                .help("Concentration units")
                .default_value("ng/mL"),
        )
        .arg(
            Arg::new("dose-normalization")
                .long("dose-normalization")
                .help("Enable dose normalization")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("covariate-analysis")
                .long("covariate-analysis")
                .help("Perform covariate analysis")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("stratify-by")
                .long("stratify-by")
                .value_name("COLUMN")
                .help("Stratify results by the specified column")
                .action(clap::ArgAction::Append),
        )
        .get_matches();

    let output_dir = PathBuf::from(matches.get_one::<String>("output").unwrap());

    // Generate example dataset if requested
    if matches.get_flag("generate-example") {
        let n_subjects: usize = matches.get_one::<String>("subjects")
            .unwrap()
            .parse()
            .expect("Invalid number of subjects");
        
        let example_file = output_dir.join("example_dataset.csv");
        std::fs::create_dir_all(&output_dir)?;
        
        ExampleDataGenerator::generate_dataset(&example_file, n_subjects)?;
        println!("Generated example dataset: {}", example_file.display());
        
        if !matches.contains_id("input") {
            // If only generating example, perform analysis on it
            return run_analysis(&example_file, &output_dir, &matches);
        }
    }

    // Run analysis on input file
    if let Some(input_file) = matches.get_one::<String>("input") {
        let input_path = PathBuf::from(input_file);
        run_analysis(&input_path, &output_dir, &matches)
    } else {
        println!("No input file specified. Use --generate-example to create sample data.");
        Ok(())
    }
}

fn run_analysis(
    input_path: &PathBuf,
    output_dir: &PathBuf,
    matches: &clap::ArgMatches,
) -> Result<()> {
    println!("Starting NCA analysis...");
    println!("Input file: {}", input_path.display());
    println!("Output directory: {}", output_dir.display());

    // Parse configuration
    let config = create_analysis_config(matches, output_dir)?;

    // Parse dataset
    println!("Parsing dataset...");
    let subjects = NonmemParser::parse_dataset(input_path)?;
    println!("Loaded {} subjects", subjects.len());

    // Perform population analysis
    println!("Performing NCA analysis...");
    let start_time = std::time::Instant::now();
    
    let results = PopulationAnalyzer::analyze_population(subjects, &config)?;
    
    let duration = start_time.elapsed();
    println!("Analysis completed in {:.2} seconds", duration.as_secs_f64());

    // Save results
    println!("Saving results...");
    OutputManager::save_results(&results, &config, output_dir)?;

    // Print summary
    print_analysis_summary(&results);

    Ok(())
}

fn create_analysis_config(
    matches: &clap::ArgMatches,
    output_dir: &PathBuf,
) -> Result<AnalysisConfig> {
    let lloq_handling = match matches.get_one::<String>("lloq-handling").unwrap().as_str() {
        "zero" => LloqHandling::Zero,
        "drop" => LloqHandling::Drop,
        "half-lloq" => LloqHandling::HalfLloq,
        _ => LloqHandling::HalfLloq,
    };

    let lambda_z_selection = match matches.get_one::<String>("lambda-z-method").unwrap().as_str() {
        "auto" => LambdaZSelection::Auto,
        "best-fit" => LambdaZSelection::BestFit { 
            min_points: 3, 
            r_squared_threshold: 0.8 
        },
        _ => LambdaZSelection::Auto,
    };

    // Get stratification columns if specified
    let stratification = if let Some(columns) = matches.get_many::<String>("stratify-by") {
        let column_names: Vec<String> = columns.cloned().collect();
        Some(StratificationConfig {
            stratify_columns: column_names,
            include_interactions: false,
            minimum_n_per_stratum: 3,
            perform_statistical_tests: true,
        })
    } else {
        None
    };

    Ok(AnalysisConfig {
        auc_methods: vec![
            AucMethod::LinearTrapezoidal,
            AucMethod::LogTrapezoidal,
            AucMethod::LinearLogTrapezoidal,
            AucMethod::LinearUpLogDown,
        ],
        lambda_z_selection,
        interpolation_method: InterpolationMethod::Linear,
        output_path: output_dir.to_string_lossy().to_string(),
        lloq_handling,
        time_units: matches.get_one::<String>("time-units").unwrap().clone(),
        concentration_units: matches.get_one::<String>("conc-units").unwrap().clone(),
        dose_normalization: matches.get_flag("dose-normalization"),
        perform_covariate_analysis: matches.get_flag("covariate-analysis"),
        stratification,
    })
}

fn print_analysis_summary(results: &PopulationResults) {
    println!("\n=== ANALYSIS SUMMARY ===");
    println!("Subjects analyzed: {}", results.individual_results.len());
    if !results.failed_subjects.is_empty() {
        println!("Failed subjects: {}", results.failed_subjects.len());
        println!("  (See failed_subjects.log for details)");
    }
    
    println!("\nKey Parameters:");
    for (param, stats) in &results.summary_statistics.parameter_stats {
        println!("  {} (Arithmetic): {:.3} ± {:.1}%", param, stats.arithmetic_mean, stats.arithmetic_cv_percent);
        if let (Some(geo_mean), Some(geo_cv)) = (stats.geometric_mean, stats.geometric_cv_percent) {
            println!("  {} (Geometric): {:.3} ± {:.1}%", param, geo_mean, geo_cv);
        }
    }
    
    println!("\nAUC Method Comparison:");
    for (method, mean_auc) in &results.method_comparison.auc_methods {
        println!("  {}: {:.3}", method, mean_auc);
    }
    
    // Print stratified results summary
    if !results.stratified_results.is_empty() {
        println!("\nStratified Analysis Summary:");
        for (stratum_key, stratum_results) in &results.stratified_results {
            println!("  {} = {}: n = {}", 
                stratum_results.stratum_name, 
                stratum_results.stratum_value, 
                stratum_results.n_subjects
            );
        }
    }
    
    // Print covariate analysis summary
    if !results.covariate_analysis.correlations.is_empty() {
        println!("\nSignificant Covariate Correlations (p < 0.05):");
        for (covariate, correlation_data) in &results.covariate_analysis.correlations {
            for (parameter, &corr_value) in &correlation_data.parameter_correlations {
                let p_value = correlation_data.p_values.get(parameter).copied().unwrap_or(1.0);
                if p_value < 0.05 {
                    println!("  {} vs {}: r = {:.3}, p = {:.3}", 
                        covariate, parameter, corr_value, p_value);
                }
            }
        }
    }
    
    println!("\nResults saved to output directory.");
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn test_example_data_generation() {
        let temp_dir = TempDir::new().unwrap();
        let example_file = temp_dir.path().join("test_dataset.csv");
        
        ExampleDataGenerator::generate_dataset(&example_file, 5).unwrap();
        assert!(example_file.exists());
    }

    #[test]
    fn test_dataset_parsing() {
        let temp_dir = TempDir::new().unwrap();
        let example_file = temp_dir.path().join("test_dataset.csv");
        
        ExampleDataGenerator::generate_dataset(&example_file, 3).unwrap();
        let subjects = NonmemParser::parse_dataset(&example_file).unwrap();
        
        assert_eq!(subjects.len(), 3);
        assert!(!subjects[0].observations.is_empty());
        assert!(!subjects[0].dosing_events.is_empty());
    }
}