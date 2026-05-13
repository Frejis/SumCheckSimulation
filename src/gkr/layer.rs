use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, Polynomial};
use ark_poly::univariate::DensePolynomial;
use crate::structures::circuit_structures::Gate;

#[derive(Clone, Debug)]
pub struct LayerReductionMessage<F: Field> {
    z1: F,
    z2: F,
    qt: DensePolynomial<F>,
}

impl<F: Field> LayerReductionMessage<F> {
    pub fn z1(&self) -> F {
        self.z1
    }

    pub fn z2(&self) -> F {
        self.z2
    }

    pub fn qt(&self) -> &DensePolynomial<F> {
        &self.qt
    }
}

impl<F: Field> LayerReductionMessage<F> {
    pub fn new(z1: F, z2: F, qt: DensePolynomial<F>) -> Self {
        Self {
            z1,
            z2,
            qt,
        }
    }
}


pub struct Layer {
    /// Gate wiring for this layer (size = 2^{k_i}).
    pub gates: Vec<Gate>
}


#[derive(Clone)]
pub struct EvaluatedLayer<F: Field> {
    pub values: Vec<F>
}

impl<F: Field> EvaluatedLayer<F> {
    pub fn new(values: Vec<F>) -> Self {
        Self { values }
    }

    pub fn value_extension(&self, s_i: usize) -> DenseMultilinearExtension<F> {
        DenseMultilinearExtension::from_evaluations_vec(s_i, self.values.clone())
    }
}

impl<F: Field> EvaluatedLayer<F> {
    pub fn empty() -> Self {
        Self { 
            values: Vec::new(), 
        }
    }
}

pub struct InputLayer<F: Field> {
    pub values: Vec<F>,
}

pub struct OutputLayer<F: Field> {
    pub values: Vec<F>,
}