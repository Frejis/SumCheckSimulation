//! Wiring predicates: the sparse MLEs `~add_i` and `~mult_i` of each layer.

use ark_ff::Field;
use ark_poly::SparseMultilinearExtension;

use crate::gkr::{GKRCircuit, GateType, InputLayer, Layer};
use crate::util::log2_pow2;

/// The add/mult wiring predicates of one circuit layer, as sparse MLEs over
/// the packed index `(g, b, c)`.
#[derive(Clone)]
pub struct LayerPredicates<F: Field> {
    pub add: SparseMultilinearExtension<F>,
    pub mult: SparseMultilinearExtension<F>,
}

/// Precomputes the add/mult wiring predicates for every layer of the circuit.
pub fn compute_predicates<F: Field>(
    circuit: &GKRCircuit<F>,
    input: &InputLayer<F>,
) -> Vec<LayerPredicates<F>> {
    circuit
        .layers
        .iter()
        .enumerate()
        .map(|(i, layer)| {
            // s_i is the number of bits addressing a gate in layer i.
            let s_i = log2_pow2(layer.gates.len());
            let next_s_i = log2_pow2(next_layer_size(circuit, input, i));
            layer_predicates(layer, s_i, next_s_i)
        })
        .collect()
}

fn layer_predicates<F: Field>(layer: &Layer, s_i: usize, next_s_i: usize) -> LayerPredicates<F> {
    let mut add_terms = Vec::<(usize, F)>::new();
    let mut mult_terms = Vec::<(usize, F)>::new();

    // The predicate index packs (g, b, c) with g in the lowest s_i bits, then
    // b, then c — so g is fixed first when fixing variables.
    for (gate_idx, gate) in layer.gates.iter().enumerate() {
        let index = gate_idx | (gate.left << s_i) | (gate.right << (s_i + next_s_i));
        match gate.predicate {
            GateType::Add => add_terms.push((index, F::one())),
            GateType::Mul => mult_terms.push((index, F::one())),
        }
    }

    let total_vars = s_i + 2 * next_s_i;
    LayerPredicates {
        add: SparseMultilinearExtension::from_evaluations(total_vars, &add_terms),
        mult: SparseMultilinearExtension::from_evaluations(total_vars, &mult_terms),
    }
}

fn next_layer_size<F: Field>(circuit: &GKRCircuit<F>, input: &InputLayer<F>, layer: usize) -> usize {
    if layer == circuit.layers.len() - 1 {
        input.values.len()
    } else {
        circuit.layers[layer + 1].gates.len()
    }
}
