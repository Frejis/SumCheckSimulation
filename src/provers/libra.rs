use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, Polynomial, SparseMultilinearExtension};
use blake2::digest::generic_array::arr;

use crate::{structures::data_structures::SumCheckProver, util::index_to_field_element};
use crate::gkr::gkr_round::GKRRound;
use crate::provers::prover_phases::ProverPhase;

struct Libra<F: Field> {
    f1: SparseMultilinearExtension<F>,
    f2: DenseMultilinearExtension<F>,
    f3: DenseMultilinearExtension<F>,
    a_hg: DenseMultilinearExtension<F>,
    g: Vec<F>,
    phase: ProverPhase,
}

impl<F: Field> Libra<F> {
    pub fn new(gkrround: &GKRRound<F>, g: Vec<F>)-> Self {
        Self {
            a_hg: DenseMultilinearExtension::from_evaluations_vec(gkrround.vj.num_vars, vec![F::zero(); 1 << gkrround.vj.num_vars]),
            phase: ProverPhase::Uninitialized,
            f1: gkrround.mult().clone(),
            f2: gkrround.vi.clone(),
            f3: gkrround.vj.clone(),
            g,
        }
    }

    pub fn initialize_phase_one(
        f1: SparseMultilinearExtension<F>,
        f3: DenseMultilinearExtension<F>,
        g: &[F],
    ) -> (DenseMultilinearExtension<F>, SparseMultilinearExtension<F>) {
        let dim = f3.num_vars;
        assert_eq!(f1.num_vars, dim * 3);
        assert_eq!(g.len(), dim);
        let mut a_hg: Vec<F> = (0..(1 << dim)).map(|_| F::zero()).collect();
        let f1_at_g = f1.fix_variables(g);
        for (xy, v) in f1_at_g.evaluations.iter() {
            let x = xy & ((1 << dim) - 1);
            let y = xy >> dim;
            a_hg[x] += *v * f3[y];
        }
        (DenseMultilinearExtension::from_evaluations_vec(dim, a_hg), f1_at_g)
    }
}

pub fn get_first_l_bits_as_usize(bits: usize, l: usize) -> usize {
    // Extract the first l bits (bits 0 to l)
    // For indices structured as (x, y) where each is l bits, this gets x
    get_bits_from_to_as_usize(bits, 0, l)
}

pub fn get_last_l_bits_as_usize(bits: usize, l: usize) -> usize {
    // Extract the last l bits (bits l to 2*l)
    // For indices structured as (x, y) where each is l bits, this gets y
    get_bits_from_to_as_usize(bits, l, 2 * l)
}

pub fn get_middle_l_bits_as_usize(bits: usize, l: usize) -> usize {
    // Extract the middle l bits (bits l to 2*l)
    // For indices structured as (x, b, c) where each is l bits, this gets b
    get_bits_from_to_as_usize(bits, l, 2 * l)
}

pub fn get_bits_from_to_as_usize(bits: usize, from: usize, to: usize) -> usize {
    let mask = (1 << (to - from)) - 1;
    (bits >> from) & mask
}

impl<F: Field> SumCheckProver<F> for Libra<F> {
    fn compute_sum(&mut self) -> F {
        self.handle_phases();
        self.sum_a_hg()
    }

    fn get_verifier_function(&self) -> ark_poly::SparseMultilinearExtension<F> {
        todo!()
    }

    fn fix_variable(&mut self, r: F) {
        let n = self.a_hg.evaluations.len();
        let half = n >> 1;

        let mut new_f2 = Vec::with_capacity(half);
        let mut new_hg = Vec::with_capacity(half);

        for i in 0..half {
            let i0 = i << 1;
            let i1 = i0 | 1;

            let f2_0 = self.f2[i0];
            let f2_1 = self.f2[i1];

            let hg_0 = self.a_hg[i0];
            let hg_1 = self.a_hg[i1];

            new_f2.push(f2_0 + r * (f2_1 - f2_0));
            new_hg.push(hg_0 + r * (hg_1 - hg_0));
        }

        self.f2 = DenseMultilinearExtension::from_evaluations_vec(
            self.f2.num_vars - 1,
            new_f2,
        );

        self.a_hg = DenseMultilinearExtension::from_evaluations_vec(
            self.a_hg.num_vars - 1,
            new_hg,
        );
    }

    fn layer_reduction_message(&self, b_star: &[F], c_star: &[F]) -> crate::gkr::layer::LayerReductionMessage<F> {
        todo!()
    }
}

impl<F: Field> Libra<F> {
    fn handle_phases(&mut self) {
        self.initialise_phase_one_if_uninitialized();
    }

    fn initialise_phase_one_if_uninitialized(&mut self) {
        if self.phase == ProverPhase::Uninitialized {
            let (arr, _) = Libra::initialize_phase_one(self.f1.clone(), self.f3.clone(), &self.g);
            self.a_hg = arr;
        }
    }

    fn sum_a_hg(&mut self) -> F {
        let dim = self.a_hg.num_vars;
        let mut sum = F::zero();
        for i in 0..1 << dim {
            sum += self.f2[i] * self.a_hg[i];
        }
        sum
    }
}

#[cfg(test)]
mod tests {
    use super::Libra;
    use ark_bls12_381::Fr;
    use ark_ff::Field;
    use ark_poly::MultilinearExtension;
    use ark_std::{test_rng, UniformRand};

    use crate::gkr::gkr_round::GKRRound;
    use crate::provers::naive::NaiveProver;
    use crate::structures::data_structures::SumCheckProver;
    use crate::util::random_gate;

    fn naive_phase_one_claim(
        gkr_round: &GKRRound<Fr>,
        g: &[Fr],
    ) -> Fr {
        
        let mut prover = NaiveProver::new(gkr_round.clone(), &g.to_vec());
        prover.compute_sum()
    }

    #[test]
    fn test_libra_initial_claim_matches_naive_reference() {
        let (gkr_round, g) = random_gkr_round_and_gate();
        let mut libra = Libra::new(&gkr_round, g.clone());
        let libra_claim = libra.compute_sum();
        let naive_claim = naive_phase_one_claim(&gkr_round, &g);
        assert_eq!(libra_claim, naive_claim, "Libra sum != naive phase-one reference");
    }

    #[test]
    fn test_libra_first_phase_identical_naive() {
        let (gkr_round, g) = random_gkr_round_and_gate();
        let mut libra = Libra::new(&gkr_round, g.clone());
        let mut naive = NaiveProver::new(gkr_round.clone(), &g.clone());
        for i in 0..libra.a_hg.num_vars {
            let random_field_element = Some(Fr::rand(&mut test_rng()));
            let libra_claim = fix_prover_and_get_new_sum(&mut libra, &random_field_element);
            let naive_claim = fix_prover_and_get_new_sum(&mut naive, &random_field_element);
            assert_eq!(libra_claim, naive_claim, "Libra sum != naive phase-one reference");
        }
    }

    fn random_gkr_round_and_gate() -> (GKRRound<Fr>, Vec<Fr>) {
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand();
        let gate = random_gate::<Fr>(gkr_round.gate_labes());
        let g = gate.clone();
        (gkr_round, g)
    }

    fn fix_prover_and_get_new_sum<F: Field>(prover: &mut impl SumCheckProver<F>, random_field_element: &Option<F>) -> F {
        if let Some(rfe) = random_field_element {
            prover.fix_variable(*rfe);
        }
        prover.compute_sum()
    }
}