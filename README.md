# Comprehensive Rust Pharmacokinetics NCA Analysis Tool

A high-performance, multi-threaded Rust implementation for non-compartmental pharmacokinetics analysis supporting NONMEM-style datasets.

## Features

### Core Analysis Capabilities
- **Individual and Population NCA Analysis**: Comprehensive parameter calculation for single subjects and population statistics
- **Multiple Dosing Routes**: Support for IV bolus, IV infusion, and oral administration
- **Robust AUC Calculation**: Four different AUC calculation methods with comparison capabilities
- **Advanced Parameter Estimation**: Automatic lambda_z selection with multiple algorithms
- **Quality Control**: Built-in validation and warning system for analysis results

### Supported Parameters
- AUC (Area Under Curve) - last, infinity, predicted
- AUMC (Area Under Moment Curve)
- Cmax/Tmax (Maximum concentration and time)
- Lambda_z (Terminal elimination rate constant)
- Half-life (Terminal elimination half-life)
- Clearance (Total body clearance)
- Volume of distribution (steady-state and terminal)
- MRT (Mean residence time)
- Bioavailability assessment

### AUC Calculation Methods
1. **Linear Trapezoidal**: Traditional linear interpolation
2. **Log Trapezoidal**: Logarithmic interpolation for declining phases
3. **Linear-Log Trapezoidal**: Phoenix WinNonlin compatible method
4. **Linear Up Log Down**: Linear for increasing, log for decreasing concentrations

### Performance Features
- **Parallel Processing**: Multi-threaded analysis using Rayon
- **Memory Efficient**: Optimized data structures and algorithms
- **Fast Execution**: Designed for rapid analysis of large datasets
- **Robust Error Handling**: Comprehensive error management and validation

## Installation

```bash
# Clone the repository
git clone <repository-url>
cd nca-analysis

# Build the project
cargo build --release

# Run tests
cargo test
```

## Usage

### Generate Example Dataset
```bash
# Generate example dataset with 20 subjects
cargo run -- --generate-example --subjects 20 --output ./results

# Generate and analyze in one step
cargo run -- --generate-example --subjects 50 --output ./results
.\target\release\nca-analysis --generate-example --subjects 50 --output ./results2
```

### Analyze Existing Dataset
```bash
# Basic analysis
cargo run -- --input dataset.csv --output ./results

# Advanced configuration
cargo run -- \
  --input dataset.csv \
  --output ./results \
  --lloq-handling half-lloq \
  --lambda-z-method best-fit \
  --time-units h \
  --conc-units ng/mL
```

### Command Line Options

- `--input, -i`: Input NONMEM dataset file
- `--output, -o`: Output directory for results (default: ./nca_results)
- `--generate-example`: Generate example dataset
- `--subjects, -n`: Number of subjects for example dataset (default: 20)
- `--lloq-handling`: LLOQ handling method (zero, drop, half-lloq)
- `--lambda-z-method`: Lambda_z selection method (auto, best-fit)
- `--time-units`: Time units for output (default: h)
- `--conc-units`: Concentration units for output (default: ng/mL)
- `--stratify`: Comma-separated stratification variables (e.g., SEX,TREATMENT)
- `--min-stratum-size`: Minimum subjects per stratum (default: 3)
- `--covariate-analysis`: Enable covariate analysis
- `--dose-normalization`: Enable dose normalization analysis
- `--include-interactions`: Include interaction analysis

## Input Dataset Format

The program expects NONMEM-style CSV datasets with the following columns:

### Required Columns
- `ID`: Subject identifier
- `TIME`: Time since first dose
- `DV`: Dependent variable (concentration)
- `AMT`: Dose amount
- `EVID`: Event ID (0=observation, 1=dose)
- `CMT`: Compartment number
- `RATE`: Infusion rate (-1=bolus, -2=oral, >0=infusion rate)

### Optional Columns
- `BLQ`: Below limit of quantification flag
- `LLOQ`: Lower limit of quantification
- `AGE`: Subject age
- `WT`: Body weight
- `HT`: Height
- `SEX`: Sex (M/F)
- `RACE`: Race/ethnicity

## Output Files

The analysis generates multiple output files:

1. **individual_results.csv**: Individual subject parameters
2. **summary_statistics.csv**: Population summary statistics
3. **method_comparison.csv**: AUC method comparison
4. **method_correlations.csv**: Correlation matrix between methods
5. **complete_results.json**: Complete results in JSON format
6. **population_summary.csv**: High-level population summary
7. **analysis_report.txt**: Comprehensive analysis report
8. **stratified_analysis.csv**: Summary of stratified analysis
9. **stratum_*.csv**: Detailed results for each stratum
10. **covariate_correlations.csv**: Covariate-parameter correlations
11. **regression_analysis.csv**: Regression analysis results
12. **dose_normalized_analysis.csv**: Dose linearity assessment

## Example Dataset

The program can generate realistic example datasets with:
- Random demographics (age, weight, height, sex, race)
- Various dosing regimens (IV bolus, infusion, oral)
- Realistic PK profiles with inter-individual variability
- Residual error and LLOQ handling
- 20 subjects by default (configurable)

## Performance Benchmarks

- **Single Subject**: ~1-5ms per subject
- **Population (100 subjects)**: ~100-500ms total
- **Large Dataset (1000 subjects)**: ~1-5 seconds total
- **Memory Usage**: <100MB for datasets up to 10,000 subjects

## Algorithm Implementation

### Lambda_z Selection
- **Auto Method**: Tests multiple point combinations, selects best R²
- **Best Fit Method**: Systematic search with minimum points and R² threshold
- **Manual Method**: User-specified time points

### AUC Calculation Robustness
- Multiple interpolation methods with automatic fallback
- Comprehensive LLOQ handling strategies
- Extrapolation percentage validation
- Cross-method validation and comparison

### Quality Control
- R² thresholds for lambda_z acceptance
- AUC extrapolation percentage limits
- Parameter range validation
- Method agreement assessment

## Dependencies

- **serde**: Serialization/deserialization
- **csv**: CSV file parsing
- **rayon**: Parallel processing
- **nalgebra**: Linear algebra operations
- **statrs**: Statistical functions
- **clap**: Command line interface
- **chrono**: Date/time handling
- **anyhow/thiserror**: Error handling

## License

MIT License - see LICENSE file for details.

## Contributing

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Submit a pull request

## Citation

If you use this tool in research, please cite:
```
Rust NCA Analysis Tool v1.0
Non-compartmental pharmacokinetics analysis implementation
```