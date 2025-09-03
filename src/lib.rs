//! Non-compartmental pharmacokinetics analysis library
//! 
//! This library provides comprehensive tools for performing individual and population
//! pharmacokinetics analysis using non-compartmental methods.

pub mod models;
pub mod parser;
pub mod nca;
pub mod auc;
pub mod parameters;
pub mod population;
pub mod output;
pub mod example_data;
pub mod errors;
pub mod stratification;
pub mod covariate;

pub use models::*;
pub use nca::*;
pub use errors::*;

/// Re-export commonly used types
pub type Result<T> = std::result::Result<T, NcaError>;