#!/bin/bash

# Enhanced example script with stratification analysis

echo "=== Enhanced Rust Pharmacokinetics NCA Analysis Tool ==="
echo "=== With Stratification and Covariate Analysis ==="
echo

# Build the project
echo "Building the project..."
cargo build --release

if [ $? -ne 0 ]; then
    echo "Build failed. Please check for compilation errors."
    exit 1
fi

# Create output directory
mkdir -p ./nca_stratified_output

echo "=== Example 1: Basic Analysis with Sex Stratification ==="
cargo run --release -- \
    --generate-example \
    --subjects 50 \
    --output ./nca_stratified_output/sex_analysis \
    --stratify SEX \
    --min-stratum-size 5 \
    --covariate-analysis \
    --dose-normalization

echo
echo "=== Example 2: Treatment and Sex Stratification with Interactions ==="
cargo run --release -- \
    --generate-example \
    --subjects 100 \
    --output ./nca_stratified_output/treatment_sex_analysis \
    --stratify TREATMENT,SEX \
    --include-interactions \
    --min-stratum-size 8 \
    --covariate-analysis \
    --dose-normalization

echo
echo "=== Example 3: Comprehensive Multi-variable Stratification ==="
cargo run --release -- \
    --generate-example \
    --subjects 200 \
    --output ./nca_stratified_output/comprehensive_analysis \
    --stratify TREATMENT,SEX,AGE_GROUP,FORMULATION \
    --include-interactions \
    --min-stratum-size 10 \
    --covariate-analysis \
    --dose-normalization \
    --lloq-handling half-lloq \
    --lambda-z-method best-fit

echo
echo "Analysis complete! Check the output directories for results:"
echo
echo "Generated files include:"
echo "- individual_results.csv (individual subject parameters)"
echo "- summary_statistics.csv (population statistics)"
echo "- stratified_analysis.csv (stratified summary)"
echo "- stratum_*.csv (detailed results by stratum)"
echo "- covariate_correlations.csv (covariate-parameter correlations)"
echo "- regression_analysis.csv (regression analysis results)"
echo "- dose_normalized_analysis.csv (dose linearity assessment)"
echo "- method_comparison.csv (AUC method comparison)"
echo "- complete_results.json (full results with all analyses)"
echo "- analysis_report.txt (comprehensive report)"

echo
echo "Example usage for your own dataset:"
echo "cargo run --release -- \\"
echo "  --input your_dataset.csv \\"
echo "  --output ./your_results \\"
echo "  --stratify SEX,TREATMENT \\"
echo "  --covariate-analysis \\"
echo "  --dose-normalization"