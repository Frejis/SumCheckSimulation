use ark_ff::Field;
use serde::{Deserialize, Serialize};
use crate::gkr::layer::{EvaluatedLayer, InputLayer, Layer};

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
    pub(crate) predicate: GateType,
}

impl Gate {
    pub fn new(left: usize, right: usize, predicate: GateType) -> Self {
        Self { left, right, predicate }
    }
}

/// A layered arithmetic circuit.
#[derive(Clone)]
pub struct GKRCircuit<F: Field> {
    pub layers: Vec<Layer>,
    pub field: F, // Dummy just here to make type checker happy.
}

impl<F: Field> GKRCircuit<F> {
    pub fn new(layers: Vec<Layer>) -> Self {
        Self { layers, field: F::zero() }
    }
}

impl<F: Field> GKRCircuit<F> {
    pub fn evaluate_circuit(&self, input: &InputLayer<F>) -> EvaluatedGKRCircuit<F> {
        let layers = self.layers.len();
        // Create n empty evaluated layers so we can fill them from bottom to top.
        let mut evaluated_layers: Vec<EvaluatedLayer<F>> = vec![EvaluatedLayer::empty(); layers];
        for i in 0..layers {
            // Working from bottom to up.
            let i = layers - i - 1;

            // Create an empty vector with 0 for each gate as the default assignment of the output of a gate
            let gate_amount = self.layers[i].gates.len();
            let mut evaluated_gates: Vec<F> = vec![F::zero(); gate_amount];

            for (idx, gate) in self.layers[i].gates.iter().enumerate() {
                let (left_child, right_child) = self.get_input_values(i, &gate, &evaluated_layers, &input);
                let value = match gate.predicate {
                    GateType::Add => left_child + right_child,
                    GateType::Mul => left_child * right_child,
                };

                evaluated_gates[idx] = value;
            }

            evaluated_layers[i] = EvaluatedLayer::new(evaluated_gates)
        }

        EvaluatedGKRCircuit::new(evaluated_layers)
    }

    fn get_input_values(&self, i: usize, gate: &Gate, computed_layers: &Vec<EvaluatedLayer<F>>, input_layer: &InputLayer<F>) -> (F, F) {
        return if i == self.layers.len() - 1 {
            // This means we are at the bottom layer at the circuit and needs the inputs to come
            // from the input layer.
            (input_layer.values[gate.left], input_layer.values[gate.right])
        } else {
            // We need to get the values from a previously computed layer.
            // Remember layer i + 1 is the layer below.
            (computed_layers[i + 1].values[gate.left], computed_layers[i + 1].values[gate.right])
        }
    }
}

#[derive(Clone)]
pub struct EvaluatedGKRCircuit<F: Field> {
    pub layers: Vec<EvaluatedLayer<F>>,
}

impl<F: Field> EvaluatedGKRCircuit<F> {
    pub fn new(layers: Vec<EvaluatedLayer<F>>) -> Self {
        Self { layers }
    }

    pub fn empty() -> Self {
        Self {
            layers: Vec::new(),
        }
    }
}