//! The per-layer instance that a sum-check execution runs over.

use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, SparseMultilinearExtension};
use ark_std::test_rng;

/// The inputs to one sum-check execution for a GKR layer `i`:
/// the add/mult wiring predicates (indexed by `(g, b, c)`) and the MLEs
/// `vi`/`vj` of the next layer's gate values (the `W_{i+1}(b)` and
/// `W_{i+1}(c)` factors of the GKR sum).
#[derive(Clone)]
pub struct GKRRound<F: Field> {
    mult_predicate: SparseMultilinearExtension<F>,
    add_predicate: SparseMultilinearExtension<F>,
    pub vi: DenseMultilinearExtension<F>,
    pub vj: DenseMultilinearExtension<F>,
    gate_labels: usize,
}

impl<F: Field> GKRRound<F> {
    pub fn new(
        mult_predicate: &SparseMultilinearExtension<F>,
        add_predicate: &SparseMultilinearExtension<F>,
        vi: &DenseMultilinearExtension<F>,
        vj: &DenseMultilinearExtension<F>,
    ) -> Self {
        Self {
            mult_predicate: mult_predicate.clone(),
            add_predicate: add_predicate.clone(),
            gate_labels: vi.num_vars,
            vi: vi.clone(),
            vj: vj.clone(),
        }
    }

    /// Random instance with `dim` variables per gate label.
    ///
    /// Only for testing/simulation: relies on the deterministic
    /// `ark_std::test_rng`. Adapted from
    /// [arkworks](https://github.com/arkworks-rs/sumcheck/blob/master/src/gkr_round_sumcheck/test.rs),
    /// modified to return both an addition and a multiplication predicate so
    /// circuits with both gate types are covered.
    pub fn new_rand(dim: usize) -> Self {
        let rng = &mut test_rng();
        let mle = DenseMultilinearExtension::rand(dim, rng);
        let mult_predicate = SparseMultilinearExtension::rand_with_config(dim * 3, 1 << dim, rng);
        let add_predicate = SparseMultilinearExtension::rand_with_config(dim * 3, 1 << dim, rng);
        Self {
            mult_predicate,
            add_predicate,
            vi: mle.clone(),
            vj: mle,
            gate_labels: dim,
        }
    }

    pub fn mult_predicate(&self) -> &SparseMultilinearExtension<F> {
        &self.mult_predicate
    }

    pub fn add_predicate(&self) -> &SparseMultilinearExtension<F> {
        &self.add_predicate
    }

    pub fn vi(&self) -> &DenseMultilinearExtension<F> {
        &self.vi
    }

    pub fn vj(&self) -> &DenseMultilinearExtension<F> {
        &self.vj
    }

    /// Number of variables in a gate label of this layer.
    pub fn gate_labels(&self) -> usize {
        self.gate_labels
    }
}
