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
        mult: SparseMultilinearExtension<F>,
        vi: DenseMultilinearExtension<F>,
        vj: DenseMultilinearExtension<F>,
        gate: Vec<F>,
    ) -> Self {
        let gkr_round = GKRRound::new(mult.clone(), vi, vj);
        Self {
            fixed_mult: mult.clone().fix_variables(&*gate),
            p: Vec::new(),
            q: Vec::new(),
            gate,
            gkr_round
        }
    }

    fn create_combined_vec_array(first_arr: &Vec<F>, last_arr: &Vec<F>) -> Vec<F> {
        let mut vec_res = Vec::with_capacity(first_arr.len()+last_arr.len());
        for i in 0..first_arr.len() {
            vec_res.push(first_arr[i]);
        }
        for i in 0..first_arr.len() {
            vec_res.push(last_arr[i]);
        }
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
        for i in 0..(1 << self.gkr_round.vi().num_vars) {
            sum += self.p[i] * self.q[i];
        }
        sum
    }

    fn get_verifier_function(&self) -> DenseMultilinearExtension<F> {
        todo!()
    }

    fn fix_variable(&mut self, random_field_element: F) {
        todo!()
    }
}


mod test {
    use ark_bls12_381::Fr;
    use ark_std::test_rng;
    use rand::random;
    use crate::data_structures::Prover;
    use crate::fast_prover::FastProver;
    use crate::naive_sum_check::NaiveProver;
    use crate::util::{random_gate, random_gkr_round_gates};

    #[test]
    fn first_phase_sum_is_identical_to_naive() {
        let mut rng = test_rng();
        let (mult, vi, vj) = random_gkr_round_gates::<Fr, _>(7, &mut rng);
        let random_gate = random_gate(7);
        let mut fast_prover = FastProver::new(mult.clone(), vi.clone(), vj.clone(), random_gate.clone());
        let naive_prover = NaiveProver::new(mult.clone(), vi.clone(), vj.clone(), random_gate.clone());
        fast_prover.initialize_phase_one();

        let naive_sum = NaiveProver::ark_compute_sum_naive(&mult, &vi, &vj, &*random_gate);
        let fast_sum = fast_prover.compute_sum();
        assert_eq!(naive_sum, fast_sum);
    }
}
