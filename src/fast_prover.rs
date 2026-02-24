use ark_bls12_381::Fr;
use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, Polynomial, SparseMultilinearExtension};
use crate::data_structures::{GKRRound, Prover};
use crate::util::{index_to_field_element, random_gate};

pub struct FastProver<F: Field> {
    fixed_mult: SparseMultilinearExtension<F>,
    gate: Vec<F>,
    gkr_round: GKRRound<F>,
    p: Vec<F>,
    q: Vec<F>,
}

impl<F: Field> FastProver<F> {
    pub fn new(
        mult: &SparseMultilinearExtension<F>,
        vi: &DenseMultilinearExtension<F>,
        vj: &DenseMultilinearExtension<F>,
        gate: &Vec<F>,
    ) -> Self {
        let gkr_round = GKRRound::new(mult, vi, vj);
        let should_initialize_phase_one = vi.num_vars > 0;
        let mut temp_res = Self {
            fixed_mult: mult.fix_variables(&*gate), // I don't think i needed to call clone.
            p: Vec::new(),
            q: Vec::new(),
            gate: gate.clone(),
            gkr_round
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

    fn fix_variable(&mut self, random_field_element: F) {
        // This fixes last bit while naive fixes first bit.
        // TODO rewrite the bit mask and test if it then gives the sum is same as naive after fixing.
        let new_size = self.p.len() >> 1;
        let mut new_p = Vec::with_capacity(new_size);
        let mut new_q = Vec::with_capacity(new_size);

        // the bit being fixed (highest bit)
        let bit = new_size; // == 1 << (num_vars - 1)

        for x_prime in 0..new_size {
            let idx0 = x_prime;        // bit = 0
            let idx1 = x_prime | bit;  // bit = 1

            let p0 = self.p[idx0];
            let p1 = self.p[idx1];

            let q0 = self.q[idx0];
            let q1 = self.q[idx1];

            new_p.push(p0 + random_field_element * (p1 - p0));
            new_q.push(q0 + random_field_element * (q1 - q0));
        }

        self.p = new_p;
        self.q = new_q;
    }
}


mod test {
    use ark_bls12_381::Fr;
    use ark_ff::Zero;
    use ark_std::{test_rng, UniformRand};
    use rand::random;
    use crate::data_structures::{GKRRound, Prover, Verifier};
    use crate::fast_prover::FastProver;
    use crate::naive_sum_check::NaiveProver;
    use crate::standard_verifier::StandardVerifier;
    use crate::util::{random_gate, random_gkr_round_gates};

    #[test]
    fn first_phase_sum_is_identical_to_naive() {
        let mut rng = test_rng();
        let (mult, vi, vj) = random_gkr_round_gates::<Fr, _>(7, &mut rng);
        let random_gate = random_gate(7);
        let mut fast_prover = FastProver::new(&mult, &vi, &vj, &random_gate);

        let naive_sum = NaiveProver::ark_compute_sum_naive(&mult, &vi, &vj, &random_gate);
        let fast_sum = fast_prover.compute_sum();
        assert_eq!(naive_sum, fast_sum);
    }

    #[test]
    fn test_fix_variable_fast_prover() {
        let mut rand = test_rng();
        let gkrr: GKRRound<Fr> = GKRRound::new_rand();
        let random_gate = random_gate(gkrr.gate_labes());
        let mut fast_prover = FastProver::new(gkrr.mult(), gkrr.vi(), gkrr.vj(), &random_gate);

        let r_field = Fr::rand(&mut rand);

        fast_prover.fix_variable(r_field);
        assert_eq!(fast_prover.p.len(), (1 << gkrr.gate_labes() - 1));
        assert_eq!(fast_prover.q.len(), (1 << gkrr.gate_labes() - 1));
    }

    #[test]
    fn test_get_verifier_function() {
        let gkrr: GKRRound<Fr> = GKRRound::new_rand();
        let random_gate = random_gate(gkrr.gate_labes());
        let fast_prover = FastProver::new(gkrr.mult(), gkrr.vi(), gkrr.vj(), &random_gate);
        // max degree is 5 for now, i think it should be 2.
        let verifier = StandardVerifier::new(5, fast_prover.compute_sum());
        let verifier_func = fast_prover.get_verifier_function();
        assert!(verifier.check_claimed_value(&verifier_func));
    }
}
