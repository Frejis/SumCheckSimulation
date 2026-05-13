use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, SparseMultilinearExtension};
use crate::structures::circuit_structures::GateType;
use crate::util::random_gkr_round_gates;

#[derive(Clone)]
pub struct GKRRound<F: Field> {
    mult_predicate: SparseMultilinearExtension<F>,
    add_predicate: SparseMultilinearExtension<F>,
    pub vi: DenseMultilinearExtension<F>,
    pub vj: DenseMultilinearExtension<F>,
    gate_labes: usize,
}

impl<F: Field> GKRRound<F> {
    pub fn set_mult_predicate(&mut self, mult_predicate: SparseMultilinearExtension<F>) {
        self.mult_predicate = mult_predicate;
    }

    pub fn set_add_predicate(&mut self, add_predicate: SparseMultilinearExtension<F>) {
        self.add_predicate = add_predicate;
    }

    pub fn set_vi(&mut self, vi: DenseMultilinearExtension<F>) {
        self.vi = vi;
    }

    pub fn set_vj(&mut self, vj: DenseMultilinearExtension<F>) {
        self.vj = vj;
    }

    pub fn set_gate_labes(&mut self, gate_labes: usize) {
        self.gate_labes = gate_labes;
    }
}

impl<F: Field> GKRRound<F> {
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

    pub fn gate_labes(&self) -> usize {
        self.gate_labes
    }
}

impl<F: Field> GKRRound<F> {
    pub fn new(
        mult_predicate: &SparseMultilinearExtension<F>,
        add_predicate: &SparseMultilinearExtension<F>,
        vi: &DenseMultilinearExtension<F>,
        vj: &DenseMultilinearExtension<F>,
    ) -> GKRRound<F> {
        GKRRound {
            mult_predicate: mult_predicate.clone(),
            add_predicate: add_predicate.clone(),
            gate_labes: vi.num_vars,
            vi: vi.clone(),
            vj: vj.clone(),
        }
    }

    /// This function should only be used for testing purposes.
    pub fn new_rand() -> GKRRound<F> {
        let (mult_pred, add_pred, vi, vj) = random_gkr_round_gates(7);
        GKRRound {
            mult_predicate: mult_pred,
            add_predicate: add_pred,
            vi,
            vj,
            gate_labes: 7,
        }
    }

    /// This function should only be used for testing purposes.
    pub fn new_rand_var_size(var_size: usize) -> GKRRound<F> {
        
        let (mult_pred, add_pred, vi, vj) = random_gkr_round_gates(var_size);
        GKRRound {
            mult_predicate: mult_pred,
            add_predicate: add_pred,
            vi,
            vj,
            gate_labes: var_size,
        }
    }

}