use ark_ff::Field;
use crate::gkr::layer::Layer;

/// Gate type: add or multiply child outputs.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum GateType {
    Add,
    Mul,
}

/// A gate at one layer, referencing children at the next layer.
#[derive(Clone)]
pub struct Gate {
    pub(crate) left: usize,
    pub(crate) right: usize,
    pub(crate) typ: GateType,
}

/// A layered arithmetic circuit.
pub struct GKRCircuit<F: Field> {
    pub layers: Vec<Layer<F>>,
}