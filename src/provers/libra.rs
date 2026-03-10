use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, Polynomial, SparseMultilinearExtension};

use crate::{structures::data_structures::SumCheckProver};
use crate::gkr::gkr_round::GKRRound;
use crate::gkr::layer::LayerReductionMessage;
use crate::provers::prover_phases::ProverPhase;
use crate::util::{interpolate_univariate, restrict_mle_to_line};

pub struct Libra<F: Field> {
    f1: SparseMultilinearExtension<F>,
    f2: DenseMultilinearExtension<F>,
    f3: DenseMultilinearExtension<F>,
    f2_clone: DenseMultilinearExtension<F>,
    a_hg: DenseMultilinearExtension<F>,
    g: Vec<F>,
    fixed_labels: Vec<F>,
    phase: ProverPhase,
    dim: usize,
}

impl<F: Field> Libra<F> {
    /// Initializes the Libra prover for a given GKR round and gate values.
    /// Initializes ´phase one´ by when creating.
    pub fn new(gkrround: &GKRRound<F>, g: Vec<F>)-> Self {
        let mut libra = Self {
            a_hg: DenseMultilinearExtension::from_evaluations_vec(gkrround.vj.num_vars, vec![F::zero(); 1 << gkrround.vj.num_vars]),
            phase: ProverPhase::Uninitialized,
            f1: gkrround.mult().clone(),
            f2: gkrround.vi.clone(),
            f3: gkrround.vj.clone(),
            f2_clone: gkrround.vi.clone(),
            dim: gkrround.vj.num_vars,
            g,
            fixed_labels: Vec::new(),
        };
        // Initialize by default.
        libra.handle_phases();
        libra
    }

    pub fn initialize_phase_one(
        f1: SparseMultilinearExtension<F>,
        f3: DenseMultilinearExtension<F>,
        g: &[F],
    ) -> (DenseMultilinearExtension<F>, SparseMultilinearExtension<F>) {
        let dim = f3.num_vars;
        //assert_eq!(f1.num_vars, dim * 3);
        //assert_eq!(g.len(), dim);
        let mut a_hg: Vec<F> = (0..(1 << dim)).map(|_| F::zero()).collect();
        let f1_at_g = f1.fix_variables(g);
        for (xy, v) in f1_at_g.evaluations.iter() {
            let x = xy & ((1 << dim) - 1);
            let y = xy >> dim;
            a_hg[x] += *v * f3[y];
        }
        (DenseMultilinearExtension::from_evaluations_vec(dim, a_hg), f1_at_g)
    }

    pub fn initialize_phase_two(
        f1: &SparseMultilinearExtension<F>,
        g: &[F],
        u: &[F], // Also the ´x´ variables of f1 that has been fixed in phase one, by the verifier.
    ) -> DenseMultilinearExtension<F> {
        // Fix the gate
        let f1_g = f1.fix_variables(g);
        // First fixed the gate, and then fixed the ´x´ variables from the verifier.
        f1_g.fix_variables(u).to_dense_multilinear_extension()
    }
}

impl<F: Field> SumCheckProver<F> for Libra<F> {
    fn compute_sum(&mut self) -> F {
        self.handle_phases();
        self.compute_inner_product()
    }

    fn get_verifier_function(&mut self) -> SparseMultilinearExtension<F> {
        match self.phase {
            ProverPhase::PhaseOne =>
                self.compute_verifier_sum(self.f2.clone()),
            ProverPhase::PhaseTwo =>
                self.compute_verifier_sum(self.f3.clone()),
            ProverPhase::Uninitialized => panic!("Phase should have been initialized by now"),
        }
    }

    fn fix_variable(&mut self, r: F) {
        assert_ne!(self.phase, ProverPhase::Uninitialized);

        self.fixed_labels.push(r);
        self.handle_phases(); // need to initialize phase two if we just fixed variable s + 1 in x.
        let (new_f2, new_hg) = self.fold_f2_and_ahg(r);
        self.update_correct_mle(new_f2);
        self.a_hg = DenseMultilinearExtension::from_evaluations_vec(self.a_hg.num_vars - 1, new_hg);
    }

    /// Roughly copy-pasted from the fast prover, but it is not tested yet.
    fn layer_reduction_message(&self, b_star: &[F], c_star: &[F]) -> LayerReductionMessage<F> {
        // TODO add a test for it at some point ig.
        let k_ip1 = self.f2_clone.num_vars;
        assert_eq!(b_star.len(), k_ip1);
        assert_eq!(b_star.len(), c_star.len());

        let ts: Vec<F> = (0..=k_ip1).map(|i| F::from(i as u64)).collect();
        let values = restrict_mle_to_line(&self.f2_clone, &b_star, &c_star, &ts);
        let g = interpolate_univariate(&values, &ts);

        LayerReductionMessage::new(g.evaluate(&F::zero()), g.evaluate(&F::one()), g)
    }
}

impl<F: Field> Libra<F> {
    fn compute_verifier_sum(&mut self, f2: DenseMultilinearExtension<F>) -> SparseMultilinearExtension<F> {
        let dim = self.a_hg.num_vars;
        let mut s0 = F::zero();
        let mut s1 = F::zero();
        for i in 0..1 << dim {
            let evaluation = match self.phase {
                ProverPhase::Uninitialized => panic!("Phase should have been initialized by now"),
                ProverPhase::PhaseOne => self.a_hg[i] * f2[i],
                ProverPhase::PhaseTwo => self.compute_phase_two_evaluated_at_i(i),
            };
            if i & 1 == 0 {
                s0 += evaluation;
            } else {
                s1 += evaluation;
            }
        }
        SparseMultilinearExtension::from_evaluations(1, vec![&(0, s0), &(1, s1)])
    }
}

impl<F: Field> Libra<F> {
    fn update_correct_mle(&mut self, new_f2: Vec<F>) {
        match self.phase {
            ProverPhase::PhaseOne =>
                self.f2 = DenseMultilinearExtension::from_evaluations_vec(self.f2.num_vars - 1, new_f2),
            ProverPhase::PhaseTwo =>
                self.f3 = DenseMultilinearExtension::from_evaluations_vec(self.f3.num_vars - 1, new_f2),
            _ => panic!("Phase should have been initialized by now"),
        }
    }
}

impl<F: Field> Libra<F> {
    fn fold_f2_and_ahg(&mut self, r: F) -> (Vec<F>, Vec<F>) {
        let n = self.a_hg.evaluations.len();
        let half = n >> 1;

        let mut new_f2 = Vec::with_capacity(half);
        let mut new_hg = Vec::with_capacity(half);

        for i in 0..half {
            self.fold_pair_at_index(r, &mut new_f2, &mut new_hg, i);
        }
        (new_f2, new_hg)
    }

    fn fold_pair_at_index(&mut self, r: F, new_f2: &mut Vec<F>, new_hg: &mut Vec<F>, index: usize) {
        match self.phase {
            ProverPhase::PhaseOne => {
                self.update_functions_from_folded_variable(r, new_f2, new_hg, index, self.f2.clone());
            }
            ProverPhase::PhaseTwo => {
                self.update_functions_from_folded_variable(r, new_f2, new_hg, index, self.f3.clone());
            }
            _ => panic!("Phase should have been initialized by now"),
        }
    }

    fn update_functions_from_folded_variable(&mut self, r: F, new_f2: &mut Vec<F>, new_hg: &mut Vec<F>, index: usize, phase_mle_func: DenseMultilinearExtension<F>) {
        let i0 = index << 1;
        let i1 = i0 | 1;

        let f2_0 = phase_mle_func[i0];
        let f2_1 = phase_mle_func[i1];

        let hg_0 = self.a_hg[i0];
        let hg_1 = self.a_hg[i1];

        let new_f2_value = f2_0 + r * (f2_1 - f2_0);
        let new_hg_value = hg_0 + r * (hg_1 - hg_0);
        Self::push_folded_values(new_f2, new_hg, new_f2_value, new_hg_value);
    }

    fn push_folded_values(new_f2: &mut Vec<F>, new_hg: &mut Vec<F>, new_f2_value: F, new_hg_value: F) {
        new_f2.push(new_f2_value);
        new_hg.push(new_hg_value);
    }
}

impl<F: Field> Libra<F> {
    fn handle_phases(&mut self) {
        self.initialize_phase_two_if_done_with_phase_one();
        self.initialise_phase_one_if_uninitialized();
    }

    fn initialize_phase_two_if_done_with_phase_one(&mut self) {
        let in_phase_one = self.phase == ProverPhase::PhaseOne;
        let done_with_phase_one = self.f3.num_vars + 1 == self.fixed_labels.len(); // E.g, we have folded all variables of f2. So x is now bounded.
        if in_phase_one && done_with_phase_one {
            self.update_self_phase_two();
        }
    }

    fn update_self_phase_two(&mut self) {
        let x_from_verifier = &self.fixed_labels[0..self.f3.num_vars];
        let arr = Libra::initialize_phase_two(&self.f1, &self.g, x_from_verifier);
        self.a_hg = arr;
        self.phase = ProverPhase::PhaseTwo;
    }

    fn initialise_phase_one_if_uninitialized(&mut self) {
        if self.phase == ProverPhase::Uninitialized {
            let (arr, _) = Libra::initialize_phase_one(self.f1.clone(), self.f3.clone(), &self.g);
            self.a_hg = arr;
            self.phase = ProverPhase::PhaseOne;
        }
    }

    /// Computes the inner product of f2 and a_hg.
    /// Which is the claim we will be working with in the sum-check protocol.
    fn compute_inner_product(&mut self) -> F {
        let dim = self.a_hg.num_vars;
        let mut sum = F::zero();
        for i in 0..1 << dim {
            if self.phase == ProverPhase::PhaseOne {
                sum += self.f2[i] * self.a_hg[i];
            } else {
                sum += self.compute_phase_two_evaluated_at_i(i);
            }
        }
        sum
    }

    fn compute_phase_two_evaluated_at_i(&mut self, i: usize) -> F {
        let u = &self.fixed_labels[0..self.dim];
        let f2_u = self.f2_clone.evaluate(&u.to_vec());
        self.f3[i] * f2_u * self.a_hg[i]
    }
}

#[cfg(test)]
mod tests {
    use super::Libra;
    use ark_bls12_381::Fr;
    use ark_ff::Field;
    use ark_std::{test_rng, UniformRand};

    use crate::gkr::gkr_round::GKRRound;
    use crate::provers::naive::NaiveProver;
    use crate::structures::data_structures::{SumCheckProver, SumCheckVerifier};
    use crate::util::random_gate;
    use crate::verifiers::standard_verifier::StandardVerifier;

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
        let (mut libra, mut naive) = get_libra_and_naive_prover(&gkr_round, &g);
        for _ in 0..libra.a_hg.num_vars {
            assert_libra_claim_identical_to_naive_in_round(&mut libra, &mut naive);
        }
    }

    #[test]
    fn test_libra_identical_to_naive_full() {
        let (gkr_round, g) = random_gkr_round_and_gate();
        let (mut libra, mut naive) = get_libra_and_naive_prover(&gkr_round, &g);
        for _ in 0..libra.f2.num_vars * 2 {
            assert_libra_claim_identical_to_naive_in_round(&mut libra, &mut naive);
        }
    }

    #[test]
    fn test_get_verifier_function_phase_one() {
        let (gkr_round, g) = random_gkr_round_and_gate();
        let mut libra = Libra::new(&gkr_round, g.clone());
        let verifier_func = libra.get_verifier_function();

        let verifier_sum = verifier_func.evaluations.iter().map(|(_, &v)| v).sum();
        let claimed_sum = libra.compute_sum();
        assert_eq!(claimed_sum, verifier_sum);
    }

    #[test]
    fn test_verifier_integration() {
        let (gkr_round, g) = random_gkr_round_and_gate();
        let mut libra = Libra::new(&gkr_round, g.clone());
        let mut verifier = StandardVerifier::new(3, libra.compute_sum(), gkr_round.clone());
        // Simulate the sum-check protocol for phase one
        for _ in 0..libra.a_hg.num_vars * 2 {
            let verifier_func = libra.get_verifier_function();
            let random_field = verifier.handle_round(&verifier_func);
            libra.fix_variable(random_field);
            verifier.set_claim(libra.compute_sum());
        }
    }

    #[test]
    fn test_verifier_integration_phase_one() {
        let (gkr_round, g) = random_gkr_round_and_gate();
        let mut libra = Libra::new(&gkr_round, g.clone());
        let mut verifier = StandardVerifier::new(3, libra.compute_sum(), gkr_round.clone());
        // Simulate the sum-check protocol for phase one
        for _ in 0..libra.a_hg.num_vars {
            let verifier_func = libra.get_verifier_function();
            let random_field = verifier.handle_round(&verifier_func);
            libra.fix_variable(random_field);
            verifier.set_claim(libra.compute_sum());
        }
    }

    fn get_libra_and_naive_prover(gkr_round: &GKRRound<Fr>, g: &Vec<Fr>) -> (Libra<Fr>, NaiveProver<Fr>) {
        let mut libra = Libra::new(&gkr_round, g.clone());
        let mut naive = NaiveProver::new(gkr_round.clone(), &g.clone());
        (libra, naive)
    }

    fn assert_libra_claim_identical_to_naive_in_round(libra: &mut Libra<Fr>, naive: &mut NaiveProver<Fr>) {
        let random_field_element = Some(Fr::rand(&mut test_rng()));
        let libra_claim = fix_prover_and_get_new_sum(libra, &random_field_element);
        let naive_claim = fix_prover_and_get_new_sum(naive, &random_field_element);
        assert_eq!(libra_claim, naive_claim, "Libra sum != naive phase-one reference");
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