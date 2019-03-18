#![feature(nll)]
extern crate pyo3;
#[macro_use]
extern crate serde_derive;

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use pyo3::{PyErr, PyResult};
use std::collections::HashMap;

use pyo3::exceptions::ValueError;
use std::fs::File;

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
    end: u64,
    score: u64,
    strand: &'a str,
    frame: &'a str,
    attributes: &'a str,
}

#[pyfunction]
//return a categorical
fn return_a_categorical() -> PyResult<HashMap<String, GTFColumn>> {
    let filename = "/project/genes.gtf";

    let mut result: HashMap<String, GTFColumn> = HashMap::new();

    let mut seqname = Categorical::new();
    let mut source = Categorical::new();
    let mut feature = Categorical::new();
    let mut start: Vec<u64> = Vec::new();
    let mut end: Vec<u64> = Vec::new();
    let mut score: Vec<f32> = Vec::new();
    let mut strand: Vec<i8> = Vec::new();
    let mut frame: Vec<i8> = Vec::new();

    result.insert("seqname".to_string(), GTFColumn::Categorical(seqname));
    result.insert("source".to_string(), GTFColumn::Categorical(source));
    result.insert("feature".to_string(), GTFColumn::Categorical(feature));
    result.insert("start".to_string(), GTFColumn::Vu64(start));
    result.insert("stop".to_string(), GTFColumn::Vu64(end));
    result.insert("score".to_string(), GTFColumn::Vf32(score));
    result.insert("strand".to_string(), GTFColumn::Vi8(strand));
    result.insert("frame".to_string(), GTFColumn::Vi8(frame));

    let mut f = File::open(filename).unwrap();
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .comment(Some(b'#'))
        .from_reader(f);
    println!("parsing");
    for result in rdr.records() {
        // The iterator yields Result<StringRecord, Error>, so we check the
        // error here.
        match result
        {
            Ok(row) => {
                let parsed_row: Result<EnsemblGTFRecord, csv::Error> = row.deserialize(None);
                match parsed_row {
                    Ok(rcd) => {
                            println!("{:?}", rcd);
                    },
                    Err(e) => return Err(PyErr::new::<ValueError, _>("CSV parsing error2"))
                }
            }
            Err(e) => return Err(PyErr::new::<ValueError, _>("CSV parsing error")),
        }
    }

    Ok(result)
}

/// This module is a python module implemented in Rust.
#[pymodule]
fn mbf_gtf(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(return_a_categorical))?;

    Ok(())
}
