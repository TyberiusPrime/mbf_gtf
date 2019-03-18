#![feature(nll)]
extern crate pyo3;
//#[macro_use]
extern crate serde_derive;
extern crate serde;

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use pyo3::{PyErr, PyResult};
use std::collections::{HashMap, HashSet};

use pyo3::exceptions::ValueError;
use std::fs::File;
use serde::{Deserialize};

use serde::{de, Deserializer};
use std::iter::FromIterator;


//use pyo3::types::PyTuple;

mod categorical;
use categorical::Categorical;

enum GTFColumn {
    Categorical(Categorical),
    Vu64(Vec<u64>),
    Vf32(Vec<f32>),
    Vi8(Vec<i8>),
}

impl IntoPyObject for GTFColumn {
    fn into_object(self, py: Python) -> PyObject {
        match self {
            GTFColumn::Categorical(s) => s.into_object(py),
            GTFColumn::Vu64(s) => s.into_object(py),
            GTFColumn::Vf32(s) => s.into_object(py),
            GTFColumn::Vi8(s) => s.into_object(py),
        }
    }
}

#[derive(Deserialize, Debug)]
struct EnsemblGTFRecord<'a> {
    seqname: &'a str,
    source: &'a str,
    feature: &'a str,
    start: u64,
    end: u64,
    #[serde(deserialize_with = "parse_nan_float")]
    score: f32,
    #[serde(deserialize_with = "parse_strand")]
    strand: i8,
    #[serde(deserialize_with = "parse_frame")]
    frame: i8,
    attributes: &'a str,
}

pub fn parse_nan_float<'de, D>(
        deserializer: D,
    ) -> Result<f32, D::Error>
    where
        D: Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        if s == "." {
            Ok(std::f32::NAN)
        }
        else {
            let f = s.parse::<f32>();
            f.map_err(de::Error::custom)
        }
    }

pub fn parse_strand<'de, D>(
        deserializer: D,
    ) -> Result<i8, D::Error>
    where
        D: Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        if s == "+" {
            Ok(1 as i8)
        }
        else if s == "-" {
            Ok(-1 as i8)
        }
        else {
            Ok(0 as i8)
        }
    }

pub fn parse_frame<'de, D>(
        deserializer: D,
    ) -> Result<i8, D::Error>
    where
        D: Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        let i = s.parse::<i8>();
        match i {
            Ok(v) => Ok(v),
            Err(_e) => Ok(-1)
        }
        
    }

fn parse_ensembl_gtf(filename: &str, include_ssf:bool,  accepted_features: Option<HashSet<String>>) -> Result<HashMap<String, GTFColumn>, csv::Error>
{

    let mut out: HashMap<String, GTFColumn> = HashMap::new();

    let mut seqname = Categorical::new();
    let mut source = Categorical::new();
    let mut feature = Categorical::new();
    let mut start: Vec<u64> = Vec::new();
    let mut end: Vec<u64> = Vec::new();
    let mut score: Vec<f32> = Vec::new();
    let mut strand: Vec<i8> = Vec::new();
    let mut frame: Vec<i8> = Vec::new();

    let f = File::open(filename)?;
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .comment(Some(b'#'))
        .has_headers(false)
        .from_reader(f);
    for result in rdr.records() {
        let row = result?;
        let parsed_row: EnsemblGTFRecord = row.deserialize(None)?;
        let keep = match &accepted_features {
            Some(af) => af.contains(parsed_row.feature),
            None => true
        };
        if keep {
            seqname.push(parsed_row.seqname);
            feature.push(parsed_row.feature);
            start.push(parsed_row.start);
            end.push(parsed_row.end);
            strand.push(parsed_row.strand);
            if include_ssf{
                source.push(parsed_row.source);
                score.push(parsed_row.score);
                frame.push(parsed_row.frame);
            }
        }
    }

    out.insert("seqname".to_string(), GTFColumn::Categorical(seqname));
    out.insert("feature".to_string(), GTFColumn::Categorical(feature));
    out.insert("start".to_string(), GTFColumn::Vu64(start));
    out.insert("stop".to_string(), GTFColumn::Vu64(end));
    out.insert("strand".to_string(), GTFColumn::Vi8(strand));
    if include_ssf {
        out.insert("source".to_string(), GTFColumn::Categorical(source));
        out.insert("score".to_string(), GTFColumn::Vf32(score));
        out.insert("frame".to_string(), GTFColumn::Vi8(frame));
    }

    Ok(out)

}

#[pyfunction]
//return a categorical
fn return_a_categorical(filename: &str, include_ssf:bool,  accepted_features: Vec<String>) -> PyResult<HashMap<String, GTFColumn>> {

    let hm_accepted_features: HashSet<String> = HashSet::from_iter(accepted_features.iter().cloned());
    let hmo_accepted_features;
    if hm_accepted_features.len() > 0{
        hmo_accepted_features = Some(hm_accepted_features);
    }
    else{
        hmo_accepted_features = None;
    }
    let parse_result = parse_ensembl_gtf(filename, include_ssf, hmo_accepted_features);
    match parse_result {
        Ok(r) => return Ok(r),
        Err(e) => return Err(PyErr::new::<ValueError, _>(e.to_string())),
    }
}

/// This module is a python module implemented in Rust.
#[pymodule]
fn mbf_gtf(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(return_a_categorical))?;

    Ok(())
}
