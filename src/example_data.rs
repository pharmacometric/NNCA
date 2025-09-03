use crate::{models::*, Result};
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;
use std::fs::File;
use std::io::Write;
use std::path::Path;

pub struct ExampleDataGenerator;

impl ExampleDataGenerator {
    pub fn generate_dataset<P: AsRef<Path>>(
        output_path: P,
        n_subjects: usize,
    ) -> Result<()> {
        let mut rng = StdRng::seed_from_u64(42); // Reproducible results
        let mut file = File::create(output_path)?;
        
        // Write header
        writeln!(file, "ID,TIME,DV,AMT,EVID,CMT,RATE,SS,II,ADDL,MDV,BLQ,LLOQ,AGE,WT,HT,SEX,RACE,TRT,STDAY,PERIOD,SEQ,FORM")?;
        
        for subject_id in 1..=n_subjects {
            let subject_data = Self::generate_subject_data(&mut rng, subject_id)?;
            Self::write_subject_data(&mut file, &subject_data)?;
        }
        
        log::info!("Generated example dataset with {} subjects", n_subjects);
        Ok(())
    }

    fn generate_subject_data(rng: &mut StdRng, subject_id: usize) -> Result<Subject> {
        // Demographics
        let age = rng.gen_range(18.0..80.0);
        let weight = rng.gen_range(50.0..120.0);
        let height = rng.gen_range(150.0..200.0);
        let sex = if rng.gen_bool(0.5) { "M" } else { "F" };
        let races = ["White", "Black", "Asian", "Hispanic"];
        let race = races[rng.gen_range(0..races.len())];
        let treatments = ["Treatment_A", "Treatment_B", "Placebo"];
        let treatment = treatments[rng.gen_range(0..treatments.len())];
        let formulations = ["Tablet", "Capsule", "Solution"];
        let formulation = formulations[rng.gen_range(0..formulations.len())];

        let demographics = Demographics {
            age: Some(age),
            weight: Some(weight),
            height: Some(height),
            sex: Some(sex.to_string()),
            race: Some(race.to_string()),
            treatment: Some(treatment.to_string()),
            study_day: Some(1),
            period: Some(rng.gen_range(1..=3)),
            sequence: Some(format!("SEQ{}", rng.gen_range(1..=4))),
            formulation: Some(formulation.to_string()),
        };

        // Generate dosing event
        let dose = rng.gen_range(10.0..500.0);
        let dosing_routes = [
            (DosingRoute::IntravenousBolus, None),
            (DosingRoute::IntravenousInfusion, Some(rng.gen_range(0.5..4.0))),
            (DosingRoute::Oral, None),
        ];
        let (route, infusion_duration) = dosing_routes[rng.gen_range(0..dosing_routes.len())].clone();

        let dosing_event = DosingEvent {
            time: 0.0,
            dose,
            route: route.clone(),
            infusion_duration,
            evid: 1,
        };

        // Generate concentration-time profile
        let observations = Self::generate_concentration_profile(rng, &route, dose, weight)?;

        Ok(Subject {
            id: subject_id.to_string(),
            observations,
            dosing_events: vec![dosing_event],
            demographics,
        })
    }

    fn generate_concentration_profile(
        rng: &mut StdRng,
        route: &DosingRoute,
        dose: f64,
        weight: f64,
    ) -> Result<Vec<Observation>> {
        let mut observations = Vec::new();
        
        // Typical PK parameters (population values with variability)
        let cl = Self::log_normal_random(rng, 10.0, 0.3) * (weight / 70.0).powf(0.75); // Allometric scaling
        let vd = Self::log_normal_random(rng, 50.0, 0.4) * (weight / 70.0);
        let ka = match route {
            DosingRoute::Oral => Self::log_normal_random(rng, 1.0, 0.5),
            _ => 0.0,
        };
        let f = match route {
            DosingRoute::Oral => Self::log_normal_random(rng, 0.8, 0.2).min(1.0),
            _ => 1.0,
        };

        // Time points
        let time_points = match route {
            DosingRoute::Oral => vec![0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 24.0, 36.0, 48.0],
            _ => vec![0.0, 0.083, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0, 48.0, 72.0],
        };

        for time in time_points {
            let concentration = Self::calculate_concentration(
                time, dose, cl, vd, ka, f, route
            );
            
            // Add residual error
            let cv_error = 0.15; // 15% CV
            let error_factor = Self::log_normal_random(rng, 1.0, cv_error);
            let final_concentration = (concentration * error_factor).max(0.0);
            
            // LLOQ handling
            let lloq = 0.1;
            let bloq = final_concentration < lloq;
            
            observations.push(Observation {
                time,
                concentration: if bloq { lloq / 2.0 } else { final_concentration },
                lloq: Some(lloq),
                bloq,
                evid: 0,
                dv: if bloq { lloq / 2.0 } else { final_concentration },
            });
        }
        
        Ok(observations)
    }

    fn calculate_concentration(
        time: f64,
        dose: f64,
        cl: f64,
        vd: f64,
        ka: f64,
        f: f64,
        route: &DosingRoute,
    ) -> f64 {
        let k = cl / vd;
        
        match route {
            DosingRoute::IntravenousBolus => {
                // One-compartment IV bolus: C = (Dose/Vd) * exp(-k*t)
                (dose / vd) * (-k * time).exp()
            }
            DosingRoute::IntravenousInfusion => {
                // Assume 1-hour infusion for simplicity
                let t_inf = 1.0;
                if time <= t_inf {
                    // During infusion: C = (Dose/t_inf/CL) * (1 - exp(-k*t))
                    (dose / t_inf / cl) * (1.0 - (-k * time).exp())
                } else {
                    // After infusion: C = (Dose/t_inf/CL) * (1 - exp(-k*t_inf)) * exp(-k*(t-t_inf))
                    let c_end_inf = (dose / t_inf / cl) * (1.0 - (-k * t_inf).exp());
                    c_end_inf * (-k * (time - t_inf)).exp()
                }
            }
            DosingRoute::Oral => {
                // One-compartment oral: C = (F*Dose*ka/(Vd*(ka-k))) * (exp(-k*t) - exp(-ka*t))
                if ka == k {
                    (f * dose / vd) * time * (-k * time).exp()
                } else {
                    (f * dose * ka) / (vd * (ka - k)) * ((-k * time).exp() - (-ka * time).exp())
                }
            }
        }
    }

    fn log_normal_random(rng: &mut StdRng, median: f64, cv: f64) -> f64 {
        let sigma = (1.0 + cv * cv).ln().sqrt();
        let mu = median.ln() - 0.5 * sigma * sigma;
        let normal_sample: f64 = rng.gen(); // This should use a proper normal distribution
        (mu + sigma * Self::inverse_normal_cdf(normal_sample)).exp()
    }

    fn inverse_normal_cdf(p: f64) -> f64 {
        // Approximation of inverse normal CDF (Box-Muller transform would be better)
        if p <= 0.0 { return f64::NEG_INFINITY; }
        if p >= 1.0 { return f64::INFINITY; }
        
        // Simple approximation - in production use a proper statistical library
        let t = (-2.0 * (1.0 - p).ln()).sqrt();
        t * if p > 0.5 { 1.0 } else { -1.0 }
    }

    fn write_subject_data(file: &mut File, subject: &Subject) -> Result<()> {
        // Write dosing record
        for dose_event in &subject.dosing_events {
            let (rate, cmt) = match (&dose_event.route, dose_event.infusion_duration) {
                (DosingRoute::IntravenousBolus, _) => (-1.0, 1),
                (DosingRoute::IntravenousInfusion, Some(duration)) => (dose_event.dose / duration, 1),
                (DosingRoute::Oral, _) => (-2.0, 1),
                _ => (0.0, 1),
            };

            writeln!(
                file,
                "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                subject.id,
                dose_event.time,
                0.0, // DV
                dose_event.dose, // AMT
                1, // EVID
                cmt, // CMT
                rate, // RATE
                0, // SS
                0, // II
                0, // ADDL
                0, // MDV
                0, // BLQ
                0.1, // LLOQ
                subject.demographics.age.unwrap_or(30.0),
                subject.demographics.weight.unwrap_or(70.0),
                subject.demographics.height.unwrap_or(170.0),
                subject.demographics.sex.as_ref().unwrap_or(&"M".to_string()),
                subject.demographics.race.as_ref().unwrap_or(&"White".to_string()),
                subject.demographics.treatment.as_ref().unwrap_or(&"Treatment_A".to_string()),
                subject.demographics.study_day.unwrap_or(1),
                subject.demographics.period.unwrap_or(1),
                subject.demographics.sequence.as_ref().unwrap_or(&"SEQ1".to_string()),
                subject.demographics.formulation.as_ref().unwrap_or(&"Tablet".to_string()),
            )?;
        }

        // Write observation records
        for obs in &subject.observations {
            writeln!(
                file,
                "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                subject.id,
                obs.time,
                obs.concentration, // DV
                0, // AMT
                0, // EVID
                1, // CMT
                0, // RATE
                0, // SS
                0, // II
                0, // ADDL
                0, // MDV
                if obs.bloq { 1 } else { 0 }, // BLQ
                obs.lloq.unwrap_or(0.1), // LLOQ
                subject.demographics.age.unwrap_or(30.0),
                subject.demographics.weight.unwrap_or(70.0),
                subject.demographics.height.unwrap_or(170.0),
                subject.demographics.sex.as_ref().unwrap_or(&"M".to_string()),
                subject.demographics.race.as_ref().unwrap_or(&"White".to_string()),
                subject.demographics.treatment.as_ref().unwrap_or(&"Treatment_A".to_string()),
                subject.demographics.study_day.unwrap_or(1),
                subject.demographics.period.unwrap_or(1),
                subject.demographics.sequence.as_ref().unwrap_or(&"SEQ1".to_string()),
                subject.demographics.formulation.as_ref().unwrap_or(&"Tablet".to_string()),
            )?;
        }
        
        Ok(())
    }
}