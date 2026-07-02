//! Layered arithmetic circuits and their evaluation.

use std::marker::PhantomData;

use ark_ff::Field;
use ark_poly::DenseMultilinearExtension;
use ark_std::rand::Rng;
use ark_std::test_rng;

use crate::util::log2_pow2;

/// Gate type: add or multiply the two child outputs.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum GateType {
    Add,
    Mul,
}

/// A gate at one layer, referencing children at the layer below.
#[derive(Clone)]
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

/// One layer of gate wiring (size must be a power of two).
#[derive(Clone)]
pub struct Layer {
    pub gates: Vec<Gate>,
}

impl Layer {
    pub fn new(gates: Vec<Gate>) -> Self {
        Self { gates }
    }
}

/// A layered arithmetic circuit. Layer 0 is the output layer; the last layer
/// reads from the input layer.
#[derive(Clone)]
pub struct GKRCircuit<F: Field> {
    pub layers: Vec<Layer>,
    _field: PhantomData<F>,
}

impl<F: Field> GKRCircuit<F> {
    pub fn new(layers: Vec<Layer>) -> Self {
        Self { layers, _field: PhantomData }
    }

    /// Generates a random layered circuit with the given sizes.
    /// Each element of `layer_sizes` must be a power of two; the last entry is
    /// the size of the input layer.
    pub fn random<R: Rng>(layer_sizes: &[usize], rng: &mut R) -> Self {
        assert!(layer_sizes.len() >= 2);
        let mut layers: Vec<Layer> = Vec::new();
        for i in 0..layer_sizes.len() - 1 {
            let gates = (0..layer_sizes[i])
                .map(|_| {
                    let predicate = if rng.gen_bool(0.5) { GateType::Add } else { GateType::Mul };
                    let left = rng.gen_range(0..layer_sizes[i + 1]);
                    let right = rng.gen_range(0..layer_sizes[i + 1]);
                    Gate::new(left, right, predicate)
                })
                .collect();
            layers.push(Layer::new(gates));
        }
        GKRCircuit::new(layers)
    }

    /// The small hand-built circuit from the report's figure (only
    /// multiplication gates).
    pub fn figure_circuit() -> Self {
        let output_layer = vec![Gate::new(0, 1, GateType::Mul), Gate::new(2, 3, GateType::Mul)];
        let middle_layer = vec![
            Gate::new(0, 0, GateType::Mul),
            Gate::new(1, 1, GateType::Mul),
            Gate::new(1, 2, GateType::Mul),
            Gate::new(3, 3, GateType::Mul),
        ];
        GKRCircuit::new(vec![Layer::new(output_layer), Layer::new(middle_layer)])
    }

    /// Evaluates every gate of the circuit bottom-up for the given input layer.
    pub fn evaluate_circuit(&self, input: &InputLayer<F>) -> EvaluatedGKRCircuit<F> {
        let layers = self.layers.len();
        // Create n empty evaluated layers so they can be filled bottom-up.
        let mut evaluated_layers: Vec<EvaluatedLayer<F>> = vec![EvaluatedLayer::empty(); layers];
        for i in (0..layers).rev() {
            let mut evaluated_gates: Vec<F> = vec![F::zero(); self.layers[i].gates.len()];

            for (idx, gate) in self.layers[i].gates.iter().enumerate() {
                let (left, right) = self.get_input_values(i, gate, &evaluated_layers, input);
                evaluated_gates[idx] = match gate.predicate {
                    GateType::Add => left + right,
                    GateType::Mul => left * right,
                };
            }

            evaluated_layers[i] = EvaluatedLayer::new(evaluated_gates);
        }

        EvaluatedGKRCircuit::new(evaluated_layers)
    }

    fn get_input_values(
        &self,
        layer: usize,
        gate: &Gate,
        evaluated_layers: &[EvaluatedLayer<F>],
        input_layer: &InputLayer<F>,
    ) -> (F, F) {
        if layer == self.layers.len() - 1 {
            // The bottom layer of the circuit reads from the input layer.
            (input_layer.values[gate.left], input_layer.values[gate.right])
        } else {
            // Layer `layer + 1` is the layer below.
            (
                evaluated_layers[layer + 1].values[gate.left],
                evaluated_layers[layer + 1].values[gate.right],
            )
        }
    }
}

/// The gate values of one evaluated layer.
#[derive(Clone)]
pub struct EvaluatedLayer<F: Field> {
    pub values: Vec<F>,
}

impl<F: Field> EvaluatedLayer<F> {
    pub fn new(values: Vec<F>) -> Self {
        Self { values }
    }

    pub fn empty() -> Self {
        Self { values: Vec::new() }
    }

    /// The multilinear extension `~W_i` of this layer's values.
    pub fn value_extension(&self) -> DenseMultilinearExtension<F> {
        let variables = log2_pow2(self.values.len());
        DenseMultilinearExtension::from_evaluations_vec(variables, self.values.clone())
    }
}

/// The gate values of every layer of an evaluated circuit.
#[derive(Clone)]
pub struct EvaluatedGKRCircuit<F: Field> {
    pub layers: Vec<EvaluatedLayer<F>>,
}

impl<F: Field> EvaluatedGKRCircuit<F> {
    pub fn new(layers: Vec<EvaluatedLayer<F>>) -> Self {
        Self { layers }
    }
}

/// The circuit's input values.
#[derive(Clone)]
pub struct InputLayer<F: Field> {
    pub values: Vec<F>,
}

impl<F: Field> InputLayer<F> {
    pub fn new(values: Vec<F>) -> Self {
        Self { values }
    }

    /// Random input layer; testing/simulation only (uses `test_rng`).
    pub fn random(input_size: usize) -> Self {
        let mut rng = test_rng();
        Self::new((0..input_size).map(|_| F::rand(&mut rng)).collect())
    }

    /// The multilinear extension of the input values.
    pub fn value_extension(&self) -> DenseMultilinearExtension<F> {
        let variables = log2_pow2(self.values.len());
        DenseMultilinearExtension::from_evaluations_vec(variables, self.values.clone())
    }
}
