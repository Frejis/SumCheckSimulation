use ark_ff::Field;

#[derive(Clone, Debug)]
pub struct LayerReductionMessage<F: Field> {
    pub z1: F,            // W(b*)
    pub z2: F,            // W(c*)
    pub q_coeffs: Vec<F>, // coefficients of q(t) = W(b* + t(c* - b*))
}

impl<F: Field> LayerReductionMessage<F> {
    pub fn new(z1: F, z2: F, q_coeffs: Vec<F>) -> Self {
        assert!(!q_coeffs.is_empty(), "q_coeffs cannot be empty");
        LayerReductionMessage { z1, z2, q_coeffs }
    }
    
}

