use crate::{models::*, errors::NcaError, Result};
use csv::ReaderBuilder;
use std::collections::HashMap;
use std::fs::File;
use std::path::Path;

pub struct NonmemParser;

impl NonmemParser {
    pub fn parse_dataset<P: AsRef<Path>>(file_path: P) -> Result<Vec<Subject>> {
        let file = File::open(file_path)?;
        let mut reader = ReaderBuilder::new()
            .has_headers(true)
            .from_reader(file);

        let mut subjects_map: HashMap<String, Subject> = HashMap::new();

        for result in reader.records() {
            let record = result?;
            let row = Self::parse_record(&record)?;
            
            let subject_id = row.get("ID")
                .ok_or_else(|| NcaError::ParseError("Missing ID column".to_string()))?
                .to_string();

            let subject = subjects_map.entry(subject_id.clone()).or_insert_with(|| Subject {
                id: subject_id.clone(),
                observations: Vec::new(),
                dosing_events: Vec::new(),
                demographics: Demographics::default(),
            });

            Self::process_row(&row, subject)?;
        }

        Ok(subjects_map.into_values().collect())
    }

    fn parse_record(record: &csv::StringRecord) -> Result<HashMap<String, String>> {
        let headers = vec![
            "ID", "TIME", "DV", "AMT", "EVID", "CMT", "RATE", "SS", "II", "ADDL",
            "MDV", "BLQ", "LLOQ", "AGE", "WT", "HT", "SEX", "RACE", "TRT", "TREAT", 
            "TREATMENT", "STDAY", "PERIOD", "SEQ", "SEQUENCE", "FORM", "FORMULATION"
        ];

        let mut row = HashMap::new();
        for (i, value) in record.iter().enumerate() {
            if i < headers.len() {
                row.insert(headers[i].to_string(), value.to_string());
            }
        }

        Ok(row)
    }

    fn process_row(row: &HashMap<String, String>, subject: &mut Subject) -> Result<()> {
        let time = Self::parse_float(row, "TIME")?;
        let evid = Self::parse_int(row, "EVID").unwrap_or(0);

        match evid {
            0 => {
                // Observation record
                let concentration = Self::parse_float(row, "DV")?;
                let lloq = Self::parse_float_optional(row, "LLOQ");
                let bloq = Self::parse_bool(row, "BLQ").unwrap_or(false);
                
                subject.observations.push(Observation {
                    time,
                    concentration,
                    lloq,
                    bloq,
                    evid,
                    dv: concentration,
                });
            }
            1 => {
                // Dosing record
                let dose = Self::parse_float(row, "AMT")?;
                let rate = Self::parse_float_optional(row, "RATE");
                
                let (route, infusion_duration) = Self::determine_dosing_route(rate, dose);
                
                subject.dosing_events.push(DosingEvent {
                    time,
                    dose,
                    route,
                    infusion_duration,
                    evid,
                });
            }
            _ => {
                // Other event types (reset, additional dose, etc.)
            }
        }

        // Update demographics if available
        Self::update_demographics(row, &mut subject.demographics)?;

        Ok(())
    }

    fn determine_dosing_route(rate: Option<f64>, dose: f64) -> (DosingRoute, Option<f64>) {
        match rate {
            Some(r) if r > 0.0 => {
                let duration = dose / r;
                (DosingRoute::IntravenousInfusion, Some(duration))
            }
            Some(-1.0) => (DosingRoute::IntravenousBolus, None),
            Some(-2.0) => (DosingRoute::Oral, None),
            _ => (DosingRoute::IntravenousBolus, None),
        }
    }

    fn update_demographics(row: &HashMap<String, String>, demographics: &mut Demographics) -> Result<()> {
        if let Some(age_str) = row.get("AGE") {
            if let Ok(age) = age_str.parse::<f64>() {
                demographics.age = Some(age);
            }
        }

        if let Some(wt_str) = row.get("WT") {
            if let Ok(weight) = wt_str.parse::<f64>() {
                demographics.weight = Some(weight);
            }
        }

        if let Some(ht_str) = row.get("HT") {
            if let Ok(height) = ht_str.parse::<f64>() {
                demographics.height = Some(height);
            }
        }

        if let Some(sex) = row.get("SEX") {
            if !sex.is_empty() {
                demographics.sex = Some(sex.clone());
            }
        }

        if let Some(race) = row.get("RACE") {
            if !race.is_empty() {
                demographics.race = Some(race.clone());
            }
        }

        // Treatment/study design variables
        for trt_col in &["TRT", "TREAT", "TREATMENT"] {
            if let Some(treatment) = row.get(*trt_col) {
                if !treatment.is_empty() {
                    demographics.treatment = Some(treatment.clone());
                    break;
                }
            }
        }

        if let Some(stday_str) = row.get("STDAY") {
            if let Ok(stday) = stday_str.parse::<i32>() {
                demographics.study_day = Some(stday);
            }
        }

        if let Some(period_str) = row.get("PERIOD") {
            if let Ok(period) = period_str.parse::<i32>() {
                demographics.period = Some(period);
            }
        }

        for seq_col in &["SEQ", "SEQUENCE"] {
            if let Some(sequence) = row.get(*seq_col) {
                if !sequence.is_empty() {
                    demographics.sequence = Some(sequence.clone());
                    break;
                }
            }
        }

        for form_col in &["FORM", "FORMULATION"] {
            if let Some(formulation) = row.get(*form_col) {
                if !formulation.is_empty() {
                    demographics.formulation = Some(formulation.clone());
                    break;
                }
            }
        }

        Ok(())
    }

    fn parse_float(row: &HashMap<String, String>, key: &str) -> Result<f64> {
        row.get(key)
            .ok_or_else(|| NcaError::ParseError(format!("Missing column: {}", key)))?
            .parse::<f64>()
            .map_err(|_| NcaError::ParseError(format!("Invalid float value for {}", key)))
    }

    fn parse_float_optional(row: &HashMap<String, String>, key: &str) -> Option<f64> {
        row.get(key)?.parse::<f64>().ok()
    }

    fn parse_int(row: &HashMap<String, String>, key: &str) -> Result<i32> {
        row.get(key)
            .ok_or_else(|| NcaError::ParseError(format!("Missing column: {}", key)))?
            .parse::<i32>()
            .map_err(|_| NcaError::ParseError(format!("Invalid integer value for {}", key)))
    }

    fn parse_bool(row: &HashMap<String, String>, key: &str) -> Option<bool> {
        match row.get(key)?.to_lowercase().as_str() {
            "1" | "true" | "yes" => Some(true),
            "0" | "false" | "no" => Some(false),
            _ => None,
        }
    }
}

impl Default for Demographics {
    fn default() -> Self {
        Self {
            age: None,
            weight: None,
            height: None,
            sex: None,
            race: None,
            treatment: None,
            study_day: None,
            period: None,
            sequence: None,
            formulation: None,
        }
    }
}