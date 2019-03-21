#![feature(nll)]
extern crate pyo3;
extern crate hashbrown;


use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use pyo3::{PyErr, PyResult};
use pyo3::exceptions::ValueError;
use std::collections::HashMap as RustHashMap;
use std::fs::File;
use std::io::BufRead;
use std::error;
use std::iter::FromIterator;

use hashbrown::{HashMap, HashSet}; //hashbrown offers a very modest speedup of about 0.7 seconds (from 10.28 to 9.5) 

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

//a helper that creates a vector, fills it with empty strings up to count
//then adds value
//similar to Categorical.new_empty_push
fn vector_new_empty_push(count: u32, value: String) -> Vec<String> {
    let mut res = Vec::new();
    res.resize(count as usize, "".to_string());
    res.push(value);
    res
}

fn inner_parse_ensembl_gtf_manual(
    filename: &str,
    accepted_features: HashSet<String>,
) -> Result<RustHashMap<String, GTFEntrys>, Box<error::Error>> {
    /*
    let input: String = std::fs::read_to_string(filename)?;
    let mut out: HashMap<String, GTFEntrys> = HashMap::new();

    let mut mode = 0;
    let mut count = 0;
    let mut offset:usize = 0;
    let mut last_start:usize = 0;
    let mut seqname: &str = "";
    let mut feature: &str = "";
    let mut start: u64 = 0;
    let mut stop: u64 = 0;
    let mut strand: i8 = 0;
    let mut attribute_name: &str = "";
    let mut attribute_value: &str = "";
    let mut attributes: HashMap<&str, &str> = HashMap::new();
    let bts = input.into_bytes();
    let last_b: u8 = 0;
    for b in bts.iter(){
        if (mode == 0) && (*b == ('\t' as u8)) { //read seqname
            seqname = std::str::from_utf8(&bts[last_start..offset])?;
            last_start = offset+1;
            mode = 1;
        }
        else if (mode == 1) && (*b == ('\t' as u8)) { //red source
            //source = std::str::from_utf8(&bts[last_start..offset])?; 
            last_start = offset+1;
            mode = 2;
        }
        else if (mode == 2) && (*b == ('\t' as u8)) { //read feature
            feature = std::str::from_utf8(&bts[last_start..offset])?;
            last_start = offset+1;
            mode = 3;
        }
        else if (mode == 3) && (*b == ('\t' as u8)) { //read start
            start = std::str::from_utf8(&bts[last_start..offset])?.parse()?;
            last_start = offset+1;
            mode = 4;
        }
        else if (mode == 4) && (*b == ('\t' as u8)) { //read stop
            stop = std::str::from_utf8(&bts[last_start..offset])?.parse()?;
            last_start = offset+1;
            mode = 5;
        }
        else if (mode == 5) && (*b == ('\t' as u8)) { //read score
            last_start = offset+1;
            mode = 6;
        }
        else if (mode == 6) && (*b == ('\t' as u8)) { //read strand
            if bts[last_start] == '+' as u8 {
                strand = 1;
            }
            else if bts[last_start] == '-' as u8 {
                strand = -1;}
            else {
                strand =  0
            }
            last_start = offset+1;
            mode = 7;
        }
        else if (mode == 7) && (*b == ('\t' as u8)) { //red frame
            last_start = offset+1;
            mode = 8;
        }
        else if (mode == 8) && (*b == ' ' as u8) {
            attribute_name = std::str::from_utf8(&bts[last_start..offset])?;
            last_start = offset+1;
            mode = 9;
        }
        else if (mode == 9) && (*b == ';' as u8) {
            attribute_value = std::str::from_utf8(&bts[last_start+1..offset-2])?;
            attributes.insert(attribute_name, attribute_value);
            last_start = offset+1;
            mode = 8;
        }
        else if (mode == 8) && (*b == ('\n' as u8)) {
            count += 1;
            mode = 0;
            last_start = offset+1;

            if !out.contains_key(feature) {
                if (accepted_features.len() > 0) && (!accepted_features.contains(feature)) {
                    continue;
                }
                let hm: GTFEntrys = GTFEntrys::new();
                out.insert(feature.to_string(), hm);
            }
            let target = out.get_mut(feature).unwrap();
            target.seqname.push(seqname);
            target.start.push(start);
            target.end.push(stop);
            target.strand.push(strand);
            }
        offset += 1;
    }
    if mode != 0 {
        println!("last line had no \n at the end");
    }
    println!("{}", count);



    Ok(out)
        */
    // this is good but it still iterates through parts of the input
    // three times!
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
            .map(|x| x.trim_start())
            .filter(|x| x.len() > 0);
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

#[pyfunction]
//return a categorical
fn parse_ensembl_gtf(
    filename: &str,
    accepted_features: Vec<String>,
) -> PyResult<RustHashMap<String, GTFEntrys>> {
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
