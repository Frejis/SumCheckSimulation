use ark_ff::Field;
use serde::{Deserialize, Serialize};
use crate::gkr::layer::Layer;

/// Gate type: add or multiply child outputs.
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum GateType {
    Add,
    Mul,
}

/// A gate at one layer, referencing children at the next layer.
#[derive(Clone, Serialize, Deserialize)]
pub struct Gate {
    pub(crate) left: usize,
    pub(crate) right: usize,
    pub(crate) typ: GateType,
}

/// A layered arithmetic circuit.
pub struct GKRCircuit<F: Field> {
    pub layers: Vec<Layer<F>>,
}