use thiserror::Error;

#[derive(Error, Debug)]
pub enum NcaError {
    #[error("Data parsing error: {0}")]
    ParseError(String),
    
    #[error("Calculation error: {0}")]
    CalculationError(String),
    
    #[error("File I/O error: {0}")]
    IoError(#[from] std::io::Error),
    
    #[error("CSV error: {0}")]
    CsvError(#[from] csv::Error),
    
    #[error("Serialization error: {0}")]
    SerdeError(#[from] serde_json::Error),
    
    #[error("Invalid dosing regimen: {0}")]
    InvalidDosing(String),
    
    #[error("Insufficient data points for calculation: {0}")]
    InsufficientData(String),
    
    #[error("Mathematical error: {0}")]
    MathError(String),
}