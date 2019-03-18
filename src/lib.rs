extern crate pyo3;

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use std::collections::HashMap;
//use pyo3::types::PyTuple;

mod categorical;
use categorical::Categorical;

#[pyfunction]
/// Formats the sum of two numbers as string
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

#[pyfunction]
//return a vector test
fn return_a_vec() -> PyResult<Vec<u32>> {
    let result = vec![23,12,1000];
    Ok(result)
}

#[pyfunction]
//return a dict test
fn return_a_dict() -> PyResult<HashMap<String, u32>> {
    let mut result: HashMap<String, u32>   = HashMap::new();
    result.insert("I can't drive".to_string(), 55);
    Ok(result)
}

#[pyfunction]
//return a categorical
fn return_a_categorical() -> PyResult<Categorical> {
    let mut result = Categorical::new();
    result.push("hello");
    result.push("how");
    result.push("are");
    result.push("you");
    result.push("hello");
    //result.insert("I can't drive".to_string(), 55);
    Ok(result)
}

/// This module is a python module implemented in Rust.
#[pymodule]
fn mbf_gtf(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(sum_as_string))?;
    m.add_wrapped(wrap_pyfunction!(return_a_vec))?;
    m.add_wrapped(wrap_pyfunction!(return_a_dict))?;
    m.add_wrapped(wrap_pyfunction!(return_a_categorical))?;

    Ok(())
}
