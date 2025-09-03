#!/bin/bash

# Example script to run the NCA analysis tool

echo "=== Rust Pharmacokinetics NCA Analysis Tool ==="
echo

# Build the project
echo "Building the project..."
cargo build --release

if [ $? -ne 0 ]; then
    echo "Build failed. Please check for compilation errors."
    exit 1
fi

# Create output directory
mkdir -p ./nca_output

# Generate example dataset and run analysis
echo "Generating example dataset with 20 subjects and running analysis..."
cargo run --release -- \
    --generate-example \
    --subjects 20 \
    --output ./nca_output \
    --lloq-handling half-lloq \
    --lambda-z-method auto \
    --time-units h \
    --conc-units ng/mL

echo
echo "Analysis complete! Check the ./nca_output directory for results."
echo
echo "Generated files:"
echo "- example_dataset.csv (input data)"
echo "- individual_results.csv (individual subject parameters)"
echo "- summary_statistics.csv (population statistics)"
echo "- method_comparison.csv (AUC method comparison)"
echo "- complete_results.json (full results)"
echo "- analysis_report.txt (comprehensive report)"

# Optional: Run analysis on custom dataset
# echo "To analyze your own dataset:"
# echo "cargo run --release -- --input your_dataset.csv --output ./your_results"