/*
This file will contain structures relevant to setting up the proof system.
*/
use rand::{Rng, RngCore};
use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, SparseMultilinearExtension};
use ark_std::test_rng;
use crate::circuit_structures::GateType;
use crate::util::random_gkr_round_gates;

pub trait Prover<F: Field> {
    // Computes the sum so we can have an alleged claim of the functions.
    fn compute_sum(&mut self) -> F;

    // Creates a function that has one variable (meaning it fixes all other variables)
    fn get_verifier_function(&self) -> DenseMultilinearExtension<F>;

    fn fix_variable(&mut self, random_field_element: F);

}

pub trait Verifier<F: Field> {
    // Has to check the degree of the function to ensure no one cheats.
    fn verify_degree(&self, fx: &DenseMultilinearExtension<F>) -> bool;

    // Returns a random field element from the verifier
    fn get_random_field_element(&mut self) -> F;

    // Takes as input a multilinear extension and checks that for each field their sum is the claim.
    fn check_claimed_value(&self, fx: &DenseMultilinearExtension<F>) -> bool;

    /// Should ideally take a function by the prover and do all necessary checks
    /// If any fails then it panics, and if everything is good then it returns a random field element.
    fn handle_round(&mut self, fx: &DenseMultilinearExtension<F>) -> F;

    fn set_claim(&mut self, claim: F);

    fn final_check(&self);
}

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
            GateType::Add
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


}

pub trait GKRRoundProver<F: Field> {
    /*
    mult: SparseMultilinearExtension<F>, // mle of mult for round k with (r, i , j)
    vi: DenseMultilinearExtension<F>, // mle of v_(k-1)(i)
    vj: DenseMultilinearExtension<F>, // mle of v_(k-1)(i)
    r: [F], // The gate "r" that is fixed.
    */
    fn set_mult(func: SparseMultilinearExtension<F>);
    fn get_mult() -> SparseMultilinearExtension<F>;
    fn set_vi(func: DenseMultilinearExtension<F>);
    fn get_vi() -> DenseMultilinearExtension<F>;
    fn set_vj(func: DenseMultilinearExtension<F>);
    fn get_vj() -> DenseMultilinearExtension<F>;

    fn set_gate(gate: [F]);
    fn get_gate() -> [F];

}