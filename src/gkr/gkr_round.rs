use rand::Rng;
use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, SparseMultilinearExtension};
use ark_std::test_rng;
use crate::structures::circuit_structures::GateType;
use crate::util::random_gkr_round_gates;

#[derive(Clone)]
pub struct GKRRound<F: Field> {
    mult: SparseMultilinearExtension<F>,
    pub(crate) vi: DenseMultilinearExtension<F>,
    pub(crate) vj: DenseMultilinearExtension<F>,
    gate_labes: usize,
    pub(crate) gate_type: GateType
}

impl<F: Field> GKRRound<F> {
    pub fn set_mult(&mut self, mult: SparseMultilinearExtension<F>) {
        self.mult = mult;
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
    pub fn mult(&self) -> &SparseMultilinearExtension<F> {
        &self.mult
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
        mult: &SparseMultilinearExtension<F>,
        vi: &DenseMultilinearExtension<F>,
        vj: &DenseMultilinearExtension<F>,
        gate_type: &GateType,
    ) -> GKRRound<F> {
        GKRRound {
            mult: mult.clone(),
            gate_labes: vi.num_vars,
            vi: vi.clone(),
            vj: vj.clone(),
            gate_type: gate_type.clone(),
        }
    }

    /// This function should only be used for testing purposes.
    pub fn new_rand() -> GKRRound<F> {

        let typ = if test_rng().r#gen::<bool>() {
            GateType::Mul // Should be an `Add` but rn fast_prover does not work with add.
        } else {
            GateType::Mul
        };
        let (mult, vi, vj) = random_gkr_round_gates(7);
        GKRRound {
            mult,
            vi,
            vj,
            gate_labes: 7,
            gate_type: typ,
        }
    }

    /// This function should only be used for testing purposes.
    pub fn new_rand_var_size(var_size: usize) -> GKRRound<F> {

        let typ = GateType::Mul; // TODO update this when the `Add` functionality gets implemented.
        let (mult, vi, vj) = random_gkr_round_gates(var_size);
        GKRRound {
            mult,
            vi,
            vj,
            gate_labes: var_size,
            gate_type: typ,
        }
    }

}