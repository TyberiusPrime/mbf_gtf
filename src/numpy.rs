use pyo3::prelude::*;
use pyo3::{PyResult};

use pyo3::types::PyDict;

pub fn numpy_from_vec_u32(input: Vec<u32>) -> PyResult<PyObject> {
    let len = input.len();

    //import numpy
    let gil = Python::acquire_gil();
    let py = gil.python();
    let locals = PyDict::new(py);
    locals.set_item("numpy", py.import("numpy")?)?;

    //create the array
    let code = format!("numpy.zeros(({},), numpy.uint32)", len);
    let rr: PyObject = py.eval(&code, None, Some(&locals))?.extract()?;
    locals.set_item("arr", rr)?;

    //figure out where the pointer is and turn it into  pointer
    let code = "arr.ctypes.data";
    let ptr_py: PyObject = py.eval(code, None, Some(&locals))?.extract()?;
    let ptr_us: usize = ptr_py.extract(py)?;
    let ptr: *mut u32 = unsafe { std::mem::transmute(ptr_us) };

    //write the data
    let da_array = unsafe { std::slice::from_raw_parts_mut(ptr, len) };
    for i in 0..len {
        da_array[i] = input[i];
    }
    //return the python object
    let result: PyObject = py.eval("arr", None, Some(&locals))?.extract()?;
    Ok(result)
}

pub fn numpy_from_vec_u64(input: Vec<u64>) -> PyResult<PyObject> {
    let len = input.len();

    //import numpy
    let gil = Python::acquire_gil();
    let py = gil.python();
    let locals = PyDict::new(py);
    locals.set_item("numpy", py.import("numpy")?)?;

    //create the array
    let code = format!("numpy.zeros(({},), numpy.uint64)", len);
    let rr: PyObject = py.eval(&code, None, Some(&locals))?.extract()?;
    locals.set_item("arr", rr)?;

    //figure out where the pointer is and turn it into  pointer
    let code = "arr.ctypes.data";
    let ptr_py: PyObject = py.eval(code, None, Some(&locals))?.extract()?;
    let ptr_us: usize = ptr_py.extract(py)?;
    let ptr: *mut u64 = unsafe { std::mem::transmute(ptr_us) };

    //write the data
    let da_array = unsafe { std::slice::from_raw_parts_mut(ptr, len) };
    for i in 0..len {
        da_array[i] = input[i];
    }
    //return the python object
    let result: PyObject = py.eval("arr", None, Some(&locals))?.extract()?;
    Ok(result)
}

pub fn numpy_from_vec_i8(input: Vec<i8>) -> PyResult<PyObject> {
    let len = input.len();

    //import numpy
    let gil = Python::acquire_gil();
    let py = gil.python();
    let locals = PyDict::new(py);
    locals.set_item("numpy", py.import("numpy")?)?;

    //create the array
    let code = format!("numpy.zeros(({},), numpy.int8)", len);
    let rr: PyObject = py.eval(&code, None, Some(&locals))?.extract()?;
    locals.set_item("arr", rr)?;

    //figure out where the pointer is and turn it into  pointer
    let code = "arr.ctypes.data";
    let ptr_py: PyObject = py.eval(code, None, Some(&locals))?.extract()?;
    let ptr_us: usize = ptr_py.extract(py)?;
    let ptr: *mut i8 = unsafe { std::mem::transmute(ptr_us) };

    //write the data
    let da_array = unsafe { std::slice::from_raw_parts_mut(ptr, len) };
    for i in 0..len {
        da_array[i] = input[i];
    }
    //return the python object
    let result: PyObject = py.eval("arr", None, Some(&locals))?.extract()?;
    Ok(result)
}

