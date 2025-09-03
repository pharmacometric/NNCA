use nca_analysis::{
    models::*,
    parser::NonmemParser,
    population::PopulationAnalyzer,
    output::OutputManager,
    example_data::ExampleDataGenerator,
};
use tempfile::TempDir;
use std::path::PathBuf;

#[test]
fn test_complete_nca_workflow() {
    // Create temporary directory
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();
    
    // Generate example dataset
    let dataset_path = temp_path.join("test_dataset.csv");
    ExampleDataGenerator::generate_dataset(&dataset_path, 5).unwrap();
    
    // Parse dataset
    let subjects = NonmemParser::parse_dataset(&dataset_path).unwrap();
    assert_eq!(subjects.len(), 5);
    
    // Create analysis configuration
    let config = AnalysisConfig {
        auc_methods: vec![
            AucMethod::LinearTrapezoidal,
            AucMethod::LogTrapezoidal,
            AucMethod::LinearLogTrapezoidal,
        ],
        lambda_z_selection: LambdaZSelection::Auto,
        interpolation_method: InterpolationMethod::Linear,
        output_path: temp_path.to_string_lossy().to_string(),
        lloq_handling: LloqHandling::HalfLloq,
        time_units: "h".to_string(),
        concentration_units: "ng/mL".to_string(),
    };
    
    // Perform analysis
    let results = PopulationAnalyzer::analyze_population(subjects, &config).unwrap();
    
    // Verify results
    assert_eq!(results.individual_results.len(), 5);
    assert!(!results.summary_statistics.parameter_stats.is_empty());
    
    // Save results
    let output_path = temp_path.join("test_output");
    OutputManager::save_results(&results, &config, &output_path).unwrap();
    
    // Verify output files exist
    assert!(output_path.join("individual_results.csv").exists());
    assert!(output_path.join("summary_statistics.csv").exists());
    assert!(output_path.join("complete_results.json").exists());
}

#[test]
fn test_auc_calculation_methods() {
    use nca_analysis::auc::AucCalculator;
    
    // Create test observations
    let observations = vec![
        Observation {
            time: 0.0,
            concentration: 100.0,
            lloq: Some(0.1),
            bloq: false,
            evid: 0,
            dv: 100.0,
        },
        Observation {
            time: 1.0,
            concentration: 75.0,
            lloq: Some(0.1),
            bloq: false,
            evid: 0,
            dv: 75.0,
        },
        Observation {
            time: 2.0,
            concentration: 50.0,
            lloq: Some(0.1),
            bloq: false,
            evid: 0,
            dv: 50.0,
        },
        Observation {
            time: 4.0,
            concentration: 25.0,
            lloq: Some(0.1),
            bloq: false,
            evid: 0,
            dv: 25.0,
        },
    ];
    
    let config = AnalysisConfig {
        auc_methods: vec![AucMethod::LinearTrapezoidal],
        lambda_z_selection: LambdaZSelection::Auto,
        interpolation_method: InterpolationMethod::Linear,
        output_path: "/tmp".to_string(),
        lloq_handling: LloqHandling::HalfLloq,
        time_units: "h".to_string(),
        concentration_units: "ng/mL".to_string(),
    };
    
    let auc_results = AucCalculator::calculate_all_methods(&observations, &config).unwrap();
    
    // Verify that we get AUC results
    assert!(auc_results.contains_key("linear_trapezoidal"));
    assert!(auc_results["linear_trapezoidal"] > 0.0);
}

#[test]
fn test_parameter_calculation() {
    use nca_analysis::parameters::ParameterCalculator;
    
    // Test Cmax/Tmax calculation
    let observations = vec![
        Observation {
            time: 0.0,
            concentration: 0.0,
            lloq: Some(0.1),
            bloq: false,
            evid: 0,
            dv: 0.0,
        },
        Observation {
            time: 1.0,
            concentration: 100.0,
            lloq: Some(0.1),
            bloq: false,
            evid: 0,
            dv: 100.0,
        },
        Observation {
            time: 2.0,
            concentration: 75.0,
            lloq: Some(0.1),
            bloq: false,
            evid: 0,
            dv: 75.0,
        },
    ];
    
    let (cmax, tmax) = ParameterCalculator::calculate_cmax_tmax(&observations).unwrap();
    assert_eq!(cmax, 100.0);
    assert_eq!(tmax, 1.0);
    
    // Test half-life calculation
    let lambda_z = 0.1;
    let half_life = ParameterCalculator::calculate_half_life(lambda_z).unwrap();
    assert!((half_life - 6.93147).abs() < 0.001);
}