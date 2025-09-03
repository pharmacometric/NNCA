use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Subject {
    pub id: String,
    pub observations: Vec<Observation>,
    pub dosing_events: Vec<DosingEvent>,
    pub demographics: Demographics,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Observation {
    pub time: f64,
    pub concentration: f64,
    pub lloq: Option<f64>,
    pub bloq: bool,
    pub evid: i32,
    pub dv: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DosingEvent {
    pub time: f64,
    pub dose: f64,
    pub route: DosingRoute,
    pub infusion_duration: Option<f64>,
    pub evid: i32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum DosingRoute {
    #[serde(rename = "IV")]
    IntravenousBolus,
    #[serde(rename = "INFUSION")]
    IntravenousInfusion,
    #[serde(rename = "ORAL")]
    Oral,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Demographics {
    pub age: Option<f64>,
    pub weight: Option<f64>,
    pub height: Option<f64>,
    pub sex: Option<String>,
    pub race: Option<String>,
    pub treatment: Option<String>,
    pub study_day: Option<i32>,
    pub period: Option<i32>,
    pub sequence: Option<String>,
    pub formulation: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NcaResults {
    pub subject_id: String,
    pub individual_parameters: IndividualParameters,
    pub method_comparisons: HashMap<String, IndividualParameters>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IndividualParameters {
    pub auc_last: Option<f64>,
    pub auc_inf: Option<f64>,
    pub auc_inf_pred: Option<f64>,
    pub auc_percent_extrap: Option<f64>,
    pub aumc_last: Option<f64>,
    pub aumc_inf: Option<f64>,
    pub cmax: Option<f64>,
    pub tmax: Option<f64>,
    pub tlast: Option<f64>,
    pub clast: Option<f64>,
    pub half_life: Option<f64>,
    pub lambda_z: Option<f64>,
    pub lambda_z_r_squared: Option<f64>,
    pub clearance: Option<f64>,
    pub volume_steady_state: Option<f64>,
    pub volume_terminal: Option<f64>,
    pub mrt: Option<f64>,
    pub bioavailability: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PopulationResults {
    pub individual_results: Vec<NcaResults>,
    pub failed_subjects: Vec<FailedSubjectAnalysis>,
    pub summary_statistics: SummaryStatistics,
    pub method_comparison: MethodComparison,
    pub stratified_results: HashMap<String, StratifiedResults>,
    pub covariate_analysis: CovariateAnalysis,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FailedSubjectAnalysis {
    pub subject_id: String,
    pub failure_reason: String,
    pub quantifiable_concentrations: usize,
    pub total_observations: usize,
    pub failed_parameters: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StratifiedResults {
    pub stratum_name: String,
    pub stratum_value: String,
    pub n_subjects: usize,
    pub individual_results: Vec<NcaResults>,
    pub summary_statistics: SummaryStatistics,
    pub method_comparison: MethodComparison,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CovariateAnalysis {
    pub correlations: HashMap<String, CovariateCorrelation>,
    pub regression_analysis: HashMap<String, RegressionResults>,
    pub dose_normalized_analysis: Option<DoseNormalizedAnalysis>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CovariateCorrelation {
    pub covariate_name: String,
    pub parameter_correlations: HashMap<String, f64>,
    pub p_values: HashMap<String, f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RegressionResults {
    pub parameter: String,
    pub covariate: String,
    pub slope: f64,
    pub intercept: f64,
    pub r_squared: f64,
    pub p_value: f64,
    pub confidence_interval: (f64, f64),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DoseNormalizedAnalysis {
    pub dose_normalized_auc: HashMap<String, ParameterStats>,
    pub dose_normalized_cmax: HashMap<String, ParameterStats>,
    pub dose_linearity_assessment: HashMap<String, LinearityAssessment>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LinearityAssessment {
    pub slope: f64,
    pub r_squared: f64,
    pub linearity_conclusion: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct StratificationConfig {
    pub stratify_columns: Vec<String>,
    pub include_interactions: bool,
    pub minimum_n_per_stratum: usize,
    pub perform_statistical_tests: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SummaryStatistics {
    pub parameter_stats: HashMap<String, ParameterStats>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParameterStats {
    pub n: usize,
    pub mean: f64,
    pub arithmetic_mean: f64,
    pub arithmetic_std: f64,
    pub arithmetic_cv_percent: f64,
    pub std: f64,
    pub cv_percent: f64,
    pub median: f64,
    pub q25: f64,
    pub q75: f64,
    pub min: f64,
    pub max: f64,
    pub geometric_mean: Option<f64>,
    pub geometric_cv_percent: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MethodComparison {
    pub auc_methods: HashMap<String, f64>,
    pub correlation_matrix: HashMap<String, HashMap<String, f64>>,
    pub bias_analysis: HashMap<String, BiasAnalysis>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BiasAnalysis {
    pub mean_difference: f64,
    pub mean_percent_difference: f64,
    pub limits_of_agreement: (f64, f64),
}

#[derive(Debug, Clone, PartialEq)]
pub struct AnalysisConfig {
    pub auc_methods: Vec<AucMethod>,
    pub lambda_z_selection: LambdaZSelection,
    pub interpolation_method: InterpolationMethod,
    pub output_path: String,
    pub lloq_handling: LloqHandling,
    pub time_units: String,
    pub concentration_units: String,
    pub stratification: Option<StratificationConfig>,
    pub perform_covariate_analysis: bool,
    pub dose_normalization: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub enum AucMethod {
    LinearTrapezoidal,
    LogTrapezoidal,
    LinearLogTrapezoidal,
    LinearUpLogDown,
}

#[derive(Debug, Clone, PartialEq)]
pub enum LambdaZSelection {
    Auto,
    Manual(Vec<usize>),
    BestFit { min_points: usize, r_squared_threshold: f64 },
}

#[derive(Debug, Clone, PartialEq)]
pub enum InterpolationMethod {
    Linear,
    LogLinear,
}

#[derive(Debug, Clone, PartialEq)]
pub enum LloqHandling {
    Zero,
    Drop,
    HalfLloq,
}