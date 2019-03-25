#![feature(nll)]
extern crate pyo3;
extern crate hashbrown;
extern crate flate2;

use std::collections::HashMap as RustHashMap;
use std::fs::File;
use std::io::{BufRead, Read, BufReader};
use std::error;
use std::iter::FromIterator;

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use pyo3::{PyErr, PyResult};
use pyo3::exceptions::ValueError;
use pyo3::types::{PyDict, PyTuple};

use hashbrown::{HashMap, HashSet}; //hashbrown offers a very modest speedup of about 0.7 seconds (from 10.28 to 9.5) 
use flate2::read::GzDecoder;

mod categorical;
mod numpy;
use categorical::Categorical;
use numpy::{
    numpy_from_vec_u64,
    numpy_from_vec_i8
};


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
    fn into_object(mut self, py: Python) -> PyObject {
        let mut hm = RustHashMap::new();
        hm.insert("seqname", self.seqname.into_object(py));
        hm.insert("start", numpy_from_vec_u64(self.start).unwrap());
                  //self.start.into_object(py));
        hm.insert("end", 
                  numpy_from_vec_u64(self.end).unwrap());
                  //self.end.into_object(py));
        hm.insert("strand", //self.strand.into_object(py));
                  numpy_from_vec_i8(self.strand).unwrap());
        let cat_attributes: RustHashMap<String, Categorical> = self.cat_attributes.drain().collect();
        let vec_attributes: RustHashMap<String, Vec<String>> = self.vec_attributes.drain().collect();
        hm.insert("cat_attributes", cat_attributes.into_object(py));
        hm.insert("vec_attributes", vec_attributes.into_object(py));
        hm.into_object(py)
    }
}

/// a helper that creates a vector, fills it with empty strings up to count
/// then adds value
/// similar to Categorical.new_empty_push
fn vector_new_empty_push(count: u32, value: String) -> Vec<String> {
    let mut res = Vec::new();
    res.resize(count as usize, "".to_string());
    res.push(value);
    res
}

fn inner_parse_ensembl_gtf(
    filename: &str,
    accepted_features: HashSet<String>,
) -> Result<RustHashMap<String, GTFEntrys>, Box<error::Error>> {
    // this is good but it still iterates through parts of the input
    // three times!
    let f = File::open(filename)?;
    let f: Box<Read> = if filename.ends_with(".gz") {
        Box::new(GzDecoder::new(f))
    } else {Box::new(f)};

    let f = BufReader::new(f);
    let mut out: HashMap<String, GTFEntrys> = HashMap::new();
    for line in f.lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let mut parts = line.splitn(9, '\t');
        let seqname = parts.next().ok_or("Failed to find seqname")?;
        parts.next(); //consume source
        let feature = parts.next().ok_or("Failed to find feature")?;
        if !out.contains_key(feature) {
            if (!accepted_features.is_empty()) && (!accepted_features.contains(feature)) {
                continue;
            }
            let hm: GTFEntrys = GTFEntrys::new();
            out.insert(feature.to_string(), hm);
        }
        let start: u64 = parts.next().ok_or("Failed to find start")?.parse()?;
        let start = start -1;
        let end: u64 = parts.next().ok_or("Failed to find start")?.parse()?;
        parts.next(); //consume score
        let strand = parts.next().ok_or("Failed to find start")?;
        let strand: i8 = if strand == "+" {
            1
        } else if strand == "-" {
            -1
        } else {
            0
        };
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
            .map(str::trim_start)
            .filter(|x| !x.is_empty());
        for attr_value in it {
            let mut kv = attr_value.splitn(2, ' ');
            let mut key: &str = kv.next().unwrap();
            if key == "tag" {
                if feature != "transcript"{ // only transcripts have tags!
                    continue;
                }
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
            if (key.starts_with("gene") & (key != "gene_id") & (feature != "gene"))
                | (key.starts_with("transcript")
                    & (key != "transcript_id")
                    & (feature != "transcript"))
            {
                continue;
            }
            let value: &str = kv.next().unwrap().trim_matches('"');
            if key.ends_with("_id") { // vec vs categorical seems to be almost performance neutral
                //just htis push here (and I guess the fill-er-up below
                //takes about 3 seconds.
                target
                    .vec_attributes
                    .get_mut(key)
                    .map(|at| {
                        at.push(value.to_string());
                    })
                    .unwrap_or_else(|| {
                        target.vec_attributes.insert(
                            key.to_string(),
                            vector_new_empty_push(target.count, value.to_string()),
                        );
                    });
            } else { // these tributes take about 1.5s to store (nd fill-er-up)
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

    let res: RustHashMap<String, GTFEntrys> = out.drain().collect();
    Ok(res)
}

/// parse_ensembl_gtf
///
/// parse a Ensembl GTF file to a dict of DataFrames
///
/// # arguments
/// `filename - A filename (uncompressed gtf)
/// `accepted_features` - a list of features to fetch, or an empty list for all 
#[pyfunction]
fn parse_ensembl_gtf(
    filename: &str,
    accepted_features: Vec<String>,
) -> PyResult<PyObject> {
    let hm_accepted_features: HashSet<String> =
        HashSet::from_iter(accepted_features.iter().cloned());
    let parse_result = inner_parse_ensembl_gtf(filename, hm_accepted_features);
    let parse_result = match parse_result {
        Ok(r) => r,
        Err(e) => return Err(PyErr::new::<ValueError, _>(e.to_string())),
    };
    let gil = Python::acquire_gil();
    let py = gil.python();
    let locals = PyDict::new(py);
    py.run("
    
def all_to_pandas(dict_of_frames):
    import pandas as pd
    result = {}
    for k, frame in dict_of_frames.items():
        d = {
            'seqname': pd.Categorical.from_codes(frame['seqname'][0], frame['seqname'][1]),
            'start': frame['start'],
            'end': frame['end'],
            'strand': frame['strand'],
        }
        for c in frame['cat_attributes']:
            x = pd.Categorical.from_codes(
                    frame['cat_attributes'][c][0], frame['cat_attributes'][c][1]
                )
            d[c] = x
        for c in frame['vec_attributes']:
            d[c] = frame['vec_attributes'][c]
        result[k] = pd.DataFrame(d)
    return result
    ", None, Some(locals))?;
    let all_to_pandas = locals.get_item("all_to_pandas").unwrap();
    let args = PyTuple::new(py, &[parse_result.into_object(py)]);
    let result = all_to_pandas.call1(args);
    match result {
        Ok(r) => Ok(r.to_object(py)),
        Err(e) => Err(e)
    }
}



/// This module is a python module implemented in Rust.
#[pymodule]
fn mbf_gtf(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(parse_ensembl_gtf))?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    Ok(())
}
