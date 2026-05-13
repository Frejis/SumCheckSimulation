use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, SparseMultilinearExtension};
use rand::Rng;
use crate::gkr::layer::Layer;
use crate::structures::circuit_structures::{GKRCircuit, Gate, GateType};

impl<F: Field> GKRCircuit<F> {
    /// Generate a random layered circuit with the given sizes.
    /// Each element of `layer_sizes` must be a power of 2.
    pub fn random<R: Rng>(layer_sizes: &[usize], rng: &mut R) -> Self {
        assert!(layer_sizes.len() >= 2);
        todo!()
    }
}