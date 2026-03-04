use ark_ff::Field;
use crate::structures::circuit_structures::Gate;

#[derive(Clone, Debug)]
pub struct LayerReductionMessage<F: Field> {
    z1: F,
    z2: F,
    qt: Vec<F>,
}

impl<F: Field> LayerReductionMessage<F> {
    pub fn z1(&self) -> F {
        self.z1
    }

    pub fn z2(&self) -> F {
        self.z2
    }

    pub fn qt(&self) -> &Vec<F> {
        &self.qt
    }
}

impl<F: Field> LayerReductionMessage<F> {
    pub fn new(z1: F, z2: F, qt: Vec<F>) -> Self {
        Self {
            z1,
            z2,
            qt,
        }
    }
}


pub struct Layer<F: Field> {
    /// Gate wiring for this layer (size = 2^{k_i}).
    pub gates: Vec<Gate>,
    /// Gate values at this layer.
    pub values: Vec<F>,
}
