use ark_ff::Field;
use crate::structures::circuit_structures::Gate;

#[derive(Clone, Debug)]
pub struct LayerReductionMessage<F: Field> {
    pub z1: F,            // W(b*)
    pub z2: F,            // W(c*)
    pub q_coeffs: Vec<F>, // coefficients of q(t) = W(b* + t(c* - b*))
}

pub struct Layer<F: Field> {
    /// Gate wiring for this layer (size = 2^{k_i}).
    pub gates: Vec<Gate>,
    /// Gate values at this layer.
    pub values: Vec<F>,
}

impl<F: Field> LayerReductionMessage<F> {
    pub fn new(z1: F, z2: F, q_coeffs: Vec<F>) -> Self {
        assert!(!q_coeffs.is_empty(), "q_coeffs cannot be empty");
        LayerReductionMessage { z1, z2, q_coeffs }
    }
    
}

