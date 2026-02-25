use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, Polynomial, SparseMultilinearExtension};
use crate::data_structures::{GKRRound, Prover};
use crate::util::{index_to_field_element};

pub struct FastProver<F: Field> {
    fixed_mult: SparseMultilinearExtension<F>,
    gate: Vec<F>,
    gkr_round: GKRRound<F>,
    p: Vec<F>,
    q: Vec<F>,
}

impl<F: Field> FastProver<F> {
    pub fn new(
        gkr_round: GKRRound<F>,
        gate: &Vec<F>,
    ) -> Self {
        let should_initialize_phase_one = gkr_round.vi.num_vars > 0;
        let mut temp_res = Self {
            fixed_mult: gkr_round.mult().fix_variables(&*gate), // I don't think i needed to call clone.
            p: Vec::new(),
            q: Vec::new(),
            gate: gate.clone(),
            gkr_round: gkr_round.clone()
        };
        if should_initialize_phase_one {
            temp_res.initialize_phase_one()
        }
        temp_res
    }

    fn create_combined_vec_array(first_arr: &Vec<F>, last_arr: &Vec<F>) -> Vec<F> {
        let mut vec_res = Vec::with_capacity(first_arr.len() + last_arr.len());
        vec_res.extend_from_slice(first_arr);
        vec_res.extend_from_slice(last_arr);
        vec_res
    }

    fn initialize_phase_one(&mut self) {
        let first_half = self.gkr_round.vi().num_vars;

        for i in 0..(1 << first_half) {
            let i_index = index_to_field_element(i, first_half);
            self.q.push(self.gkr_round.vi().evaluate(&i_index));
            // now do the sum for j.
            for j in 0..(1 << self.gkr_round.vj().num_vars) {
                let j_index = index_to_field_element(j, self.gkr_round.vj().num_vars);
                let vj_value = self.gkr_round.vj().evaluate(&j_index);

                let combined_vec = Self::create_combined_vec_array(&i_index, &j_index);

                let mult_val = self.fixed_mult.evaluate(&combined_vec);
                let res = mult_val * vj_value;
                if self.p.len() < (1 << first_half) {
                    self.p.push(res);
                } else {
                    self.p[i] += res;
                }
            }
        }
    }

}

impl<F: Field> Prover<F> for FastProver<F> {
    fn compute_sum(&self) -> F {
        // This currently only works for the first half.
        let mut sum = F::zero();
        for i in 0..self.p.len() {
            sum += self.p[i] * self.q[i];
        }
        sum
    }

    fn get_verifier_function(&self) -> DenseMultilinearExtension<F> {
        let mut s0 = F::zero();
        let mut s1 = F::zero();

        for mask in 0..self.p.len() {

            let value = self.p[mask] * self.q[mask];
            if (mask & 1 == 0) {
                s0 += &value;
            } else {
                s1 += &value;
            }
        }

        DenseMultilinearExtension::from_evaluations_vec(1, vec![s0, s1])
    }

    fn fix_variable(&mut self, r: F) {
        let n = self.p.len();
        assert_eq!(n % 2, 0);
        let half = n >> 1;
        let mut new_p = Vec::with_capacity(half);
        let mut new_q = Vec::with_capacity(half);

        for i in 0..half {
            // LSB bit = 0 hver gang.
            let index_0 = i << 1;
            let index_1 = index_0 | 1;
            let a_0 = self.p[index_0];
            let a_1 = self.p[index_1];
            let b_0 = self.q[index_0];
            let b_1 = self.q[index_1];
            new_p.push(a_0 + r * (a_1 - a_0));
            new_q.push(b_0 + r * (b_1 - b_0));
        }

        self.p = new_p;
        self.q = new_q;
    }
}


mod test {
    use ark_bls12_381::Fr;
    use ark_ff::Zero;
    use ark_std::{test_rng, UniformRand};
    use crate::data_structures::{GKRRound, Prover, Verifier};
    use crate::fast_prover::FastProver;
    use crate::naive_sum_check::NaiveProver;
    use crate::standard_verifier::StandardVerifier;
    use crate::util::{random_gate, random_gkr_round_gates};

    #[test]
    fn first_phase_sum_is_identical_to_naive() {
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand();
        let random_gate = random_gate(gkr_round.gate_labes());
        let fast_prover = FastProver::new(gkr_round.clone(), &random_gate);

        let naive_sum = NaiveProver::ark_compute_sum_naive(&gkr_round.mult(), &gkr_round.vi, &gkr_round.vj, &random_gate);
        let fast_sum = fast_prover.compute_sum();
        assert_eq!(naive_sum, fast_sum);
    }

    #[test]
    fn test_fix_variable_reduces_amount_of_variables() {
        let mut rand = test_rng();
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand();
        let random_gate = random_gate(gkr_round.gate_labes());
        let mut fast_prover = FastProver::new(gkr_round.clone(), &random_gate);

        let r_field = Fr::rand(&mut rand);

        fast_prover.fix_variable(r_field);
        assert_eq!(fast_prover.p.len(), (1 << gkr_round.gate_labes() - 1));
        assert_eq!(fast_prover.q.len(), (1 << gkr_round.gate_labes() - 1));
    }

    #[test]
    fn test_verifier_func_value() {
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand();
        let random_gate = random_gate(gkr_round.gate_labes());
        let mut prover = FastProver::new(gkr_round.clone(), &random_gate);

        let verifier_func = prover.get_verifier_function();
        prover.fix_variable(Fr::zero());
        assert_eq!(prover.compute_sum(), verifier_func[0]);
    }

    #[test]
    fn test_fix_variable_same_as_naive() {
        let mut rand = test_rng();
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand();
        let random_gate = random_gate(gkr_round.gate_labes());
        let mut fast_prover = FastProver::new(gkr_round.clone(), &random_gate);

        let mut naive_prover = NaiveProver::new(gkr_round, &random_gate);

        let r_field = Fr::rand(&mut rand);

        fast_prover.fix_variable(r_field);
        naive_prover.fix_variable(r_field);
        let fast_sum = fast_prover.compute_sum();
        let naive_sum = naive_prover.compute_sum();
        assert_eq!(fast_sum, naive_sum);
    }

    #[test]
    fn test_get_verifier_function() {
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand();
        let random_gate = random_gate(gkr_round.gate_labes());
        let mut fast_prover = FastProver::new(gkr_round.clone(), &random_gate);
        // max degree is 5 for now, i think it should be 2.
        let verifier = StandardVerifier::new(5, fast_prover.compute_sum(), gkr_round);
        let verifier_func = fast_prover.get_verifier_function();
        assert!(verifier.check_claimed_value(&verifier_func));
    }
}