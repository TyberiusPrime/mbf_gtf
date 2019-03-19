#![feature(nll)]
extern crate pyo3;
//#[macro_use]
extern crate serde;
extern crate serde_derive;

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use pyo3::{PyErr, PyResult};
use std::collections::{HashMap, HashSet};
use std::io::BufRead;

use pyo3::exceptions::ValueError;
use serde::Deserialize;
use std::fs::File;

use serde::{de, Deserializer};
use std::iter::FromIterator;
use std::error;

//use pyo3::types::PyTuple;

mod categorical;
use categorical::Categorical;

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

pub fn parse_nan_float<'de, D>(deserializer: D) -> Result<f32, D::Error>
where
    D: Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;
    if s == "." {
        Ok(std::f32::NAN)
    } else {
        let f = s.parse::<f32>();
        f.map_err(de::Error::custom)
    }
}

pub fn parse_strand<'de, D>(deserializer: D) -> Result<i8, D::Error>
where
    D: Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;
    if s == "+" {
        Ok(1 as i8)
    } else if s == "-" {
        Ok(-1 as i8)
    } else {
        Ok(0 as i8)
    }
}

pub fn parse_frame<'de, D>(deserializer: D) -> Result<i8, D::Error>
where
    D: Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;
    let i = s.parse::<i8>();
    match i {
        Ok(v) => Ok(v),
        Err(_e) => Ok(-1),
    }
}

struct GTFEntrys {
    seqname: Categorical,
    start: Vec<u64>,
    end: Vec<u64>,
    strand: Vec<i8>,
    cat_attributes: HashMap<String, Categorical>,
    vec_attributes: HashMap<String, Vec<String>>,
    count: u32,
}

impl GTFEntrys {
    pub fn new() -> GTFEntrys {
        GTFEntrys {
            seqname: Categorical::new(),
            start: Vec::new(),
            end: Vec::new(),
            strand: Vec::new(),
            cat_attributes: HashMap::new(),
            vec_attributes: HashMap::new(),
            count: 0,
        }
    }
}

impl IntoPyObject for GTFEntrys {
    fn into_object(self, py: Python) -> PyObject {
        let mut hm = HashMap::new();
        hm.insert("seqname", self.seqname.into_object(py));
        hm.insert("start", self.start.into_object(py));
        hm.insert("end", self.end.into_object(py));
        hm.insert("strand", self.strand.into_object(py));
        hm.insert("cat_attributes", self.cat_attributes.into_object(py));
        hm.insert("vec_attributes", self.vec_attributes.into_object(py));
        hm.into_object(py)
    }
}

fn new_resized_string_vector(count: u32, value: String) -> Vec<String> {
    let mut res = Vec::new();
    res.resize(count as usize, "".to_string());
    res.push(value);
    res
}

fn inner_parse_ensembl_gtf_manual(filename: &str, accepted_features: HashSet<String>) 
        -> Result<HashMap<String, GTFEntrys>, Box<error::Error>> 
{
    let f = File::open(filename)?;
    let f = std::io::BufReader::new(f);
    let mut out: HashMap<String, GTFEntrys> = HashMap::new();
    for line in f.lines() {
        let line = line?;
        if line.starts_with('#') || line.len() == 0 {
            continue;
        }
        let mut parts = line.splitn(9, "\t");
        let seqname = parts.next().ok_or("Failed to find seqname")?;
        parts.next(); //consume source
        let feature = parts.next().ok_or("Failed to find feature")?;
        if !out.contains_key(feature) {
            if (accepted_features.len() > 0) && (!accepted_features.contains(feature)) {
                continue;
            }
            let hm: GTFEntrys = GTFEntrys::new();
            out.insert(feature.to_string(), hm);
        }
        let start: u64 = parts.next().ok_or("Failed to find start")?.parse()?;
        let end: u64 = parts.next().ok_or("Failed to find start")?.parse()?;
        parts.next(); //consume score
        let strand = parts.next().ok_or("Failed to find start")?;
        let strand:i8 = if strand == "+" { 1 } else if strand == "-" {-1} else {0};
        let mut target = out.get_mut(feature).unwrap();
        target.seqname.push(seqname);
        target.start.push(start);
        target.end.push(end);
        target.strand.push(strand);
        let mut tag_count = 0;
        parts.next(); //consume frame
        let attributes = parts.next().ok_or("Failed to find attributes")?;
        let it = attributes
                .split_terminator(';')
                .map(|x| x.trim_start())
                .filter(|x| x.len() > 0);
        for attr_value in it {
            let mut kv = attr_value.splitn(2, ' ');
            let mut key = kv.next().unwrap();
            if key == "tag" {
                if tag_count == 0 {
                    key = "tag0"
                } else if tag_count == 1 {
                    key = "tag1"
                } else if tag_count == 2 {
                    key = "tag2"
                    } else if tag_count == 3 {
                    key = "tag3"
                    } else if tag_count == 4 {
                    key = "tag4"
                    } else if tag_count == 5 {
                    key = "tag5"
                } else {
                    continue; // silently swallow further tags
                }
                tag_count += 1;
            }
            /*if (feature == "exon")
                & ((key == "gene_biotype")
                    | (key == "gene_version")
                    | (key == "transcript_version")
                    | (key == "gene_name")
                    | (key == "gene_source")
                    | (key == "transcript_name")
                    | (key == "transcript_source")
                    | (key == "transcript_biotype")
                    | (key == "transcript_support_level"))
            {
                continue;
            } else if (feature == "transcript")
                & ((key == "gene_version")
                    | (key == "gene_name")
                    | (key == "gene_source")
                    | (key == "gene_biotype")
                    | (key == "transcript_source"))
            {
                continue;
            }
            */

            let value = kv.next().unwrap().trim_matches('"');
            if false & key.ends_with("_id") {
                target
                    .vec_attributes
                    .get_mut(key)
                    .map(|at| {
                        at.push(value.to_string());
                    })
                    .unwrap_or_else(|| {
                        target.vec_attributes.insert(
                            key.to_string(),
                            new_resized_string_vector(target.count, value.to_string()),
                        );
                    });
            } else {
                target
                    .cat_attributes
                    .get_mut(key)
                    .map(|at| {
                        at.push(value);
                    })
                    .unwrap_or_else(|| {
                        target.cat_attributes.insert(
                            key.to_string(),
                            Categorical::new_empty_push(target.count, value),
                        );
                    });
            }
        }
        target.count += 1;
        for (_key, value) in target.cat_attributes.iter_mut() {
            if (value.len() as u32) < target.count {
                value.push("");
            }
        }
        for (_key, value) in target.vec_attributes.iter_mut() {
            if (value.len() as u32) < target.count {
                value.push("".to_string());
            }
        }
    }


    Ok(out)
}

#[pyfunction]
//return a categorical
fn parse_ensembl_gtf(
    filename: &str,
    accepted_features: Vec<String>,
) -> PyResult<HashMap<String, GTFEntrys>> {
    let hm_accepted_features: HashSet<String> =
        HashSet::from_iter(accepted_features.iter().cloned());
    let parse_result = inner_parse_ensembl_gtf_manual(filename, hm_accepted_features);
    match parse_result {
        Ok(r) => return Ok(r),
        Err(e) => return Err(PyErr::new::<ValueError, _>(e.to_string())),
    }
}

/// This module is a python module implemented in Rust.
#[pymodule]
fn mbf_gtf(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(parse_ensembl_gtf))?;

    Ok(())
}
