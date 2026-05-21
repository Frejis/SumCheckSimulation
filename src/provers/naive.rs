use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, Polynomial, SparseMultilinearExtension};
use crate::gkr::gkr_round::GKRRound;
use crate::gkr::layer::LayerReductionMessage;
use crate::provers::fast::FastProver;
use crate::structures::data_structures::SumCheckProver;
use crate::util::{index_to_field_element};

pub struct NaiveProver<F: Field> {
    gkr_round: GKRRound<F>,
    // The fixed multiplication predicate for a given gate.
    fixed_mult: SparseMultilinearExtension<F>,
    // The fixed addition gate for a given gate.
    fixed_add: SparseMultilinearExtension<F>,
    fixed_variables: Vec<F>,
    layer_value_mle: DenseMultilinearExtension<F>,
}

impl<F: Field> NaiveProver<F> {
    pub fn new(
        gkr_round: GKRRound<F>,
        gate: &[F],
    ) -> NaiveProver<F> {
        let fixed_mult = gkr_round.mult_predicate().fix_variables(gate);
        let fixed_add = gkr_round.add_predicate().fix_variables(gate);
        NaiveProver {
            layer_value_mle: gkr_round.vi.clone(),
            gkr_round,
            fixed_mult,
            fixed_add,
            fixed_variables: Vec::new(),
        }
    }
    
    fn calculate_sum_naive(
        &self,
    ) -> F {
        let vi_variables = self.gkr_round.vi().num_vars;
        let mut sum_xy = F::zero();
        for x in 0..(1 << vi_variables) {
            let f2_x = self.gkr_round.vi()[x];
            let mult_f1_gx = self.fixed_mult
                .fix_variables(&index_to_field_element(x, vi_variables))
                .to_dense_multilinear_extension();
            let add_f1_gx = self.fixed_add
                .fix_variables(&index_to_field_element(x, vi_variables))
                .to_dense_multilinear_extension();
            for y in 0..(1 << self.gkr_round.vj.num_vars) {
                // Adding ~add_i(g,b,c) ~W_i(b)
                sum_xy += add_f1_gx[y] * f2_x;
                // Adding ~add_i(g,b,c) ~W_i(c)
                sum_xy += add_f1_gx[y] * self.gkr_round.vj[y];
                // Adding the term for mult predicate
                sum_xy += mult_f1_gx[y] * f2_x * self.gkr_round.vj[y];
            }
        }
        sum_xy
    }

    pub fn ark_compute_sum_naive(
        mult_predicate: &SparseMultilinearExtension<F>,
        add_predicate: &SparseMultilinearExtension<F>,
        f2: &DenseMultilinearExtension<F>,
        f3: &DenseMultilinearExtension<F>,
        g: &[F],
    ) -> F {
        let dim = f2.num_vars;
        let mult_predicate_at_gate = mult_predicate.fix_variables(g);
        let add_predicate_at_gate = add_predicate.fix_variables(g);
        let mut sum_xy = F::zero();
        for x in 0..(1 << dim) {
            let f2_x = f2[x];
            let mf1_gx = mult_predicate_at_gate
                .fix_variables(&index_to_field_element(x, dim))
                .to_dense_multilinear_extension();
            let  af1_gx = add_predicate_at_gate
                .fix_variables(&index_to_field_element(x, dim))
                .to_dense_multilinear_extension();
            for y in 0..(1 << dim) {
                let test = af1_gx[y];
                let fst_add_term = test * f2[x];
                let snd_add_term = test * f3[y];
                sum_xy += mf1_gx[y] * f2_x * f3[y] + fst_add_term + snd_add_term;
            }
        }
        sum_xy
    }

    pub fn ark_compute_sum_alr_fixed(
        mult_predicate: &SparseMultilinearExtension<F>,
        add_predicate: &SparseMultilinearExtension<F>,
        f2: &DenseMultilinearExtension<F>,
        f3: &DenseMultilinearExtension<F>,
        g: &[F],
        f: &F,
    ) -> F {
        let mult_predicate_at_gate = mult_predicate.fix_variables(g);
        let add_predicate_at_gate = add_predicate.fix_variables(g);
        let mut sum_xy = F::zero();
        for x in 0..(1 << f2.num_vars - 1) {
            let f2_x = f2.fix_variables(&[*f])[x];
            let mf1_gx = mult_predicate_at_gate
                .fix_variables(&[*f])
                .fix_variables(&index_to_field_element(x, f2.num_vars - 1))
                .to_dense_multilinear_extension();
            let  af1_gx = add_predicate_at_gate
                .fix_variables(&[*f])
                .fix_variables(&index_to_field_element(x, f2.num_vars - 1))
                .to_dense_multilinear_extension();
            for y in 0..(1 << f3.num_vars) {
                let fst_add_term = af1_gx[y] * f2.fix_variables(&[*f])[x];
                let snd_add_term = af1_gx[y] * f3[y];
                sum_xy += mf1_gx[y] * f2_x * f3[y] + fst_add_term + snd_add_term;
            }
        }
        sum_xy
    }
}

impl<F: Field> SumCheckProver<F> for NaiveProver<F> {
    // Needs to be refactored just my last Sum-check which I know works.
    fn compute_sum(&mut self) -> F {
        self.calculate_sum_naive()
    }

    fn  get_verifier_function(&mut self) -> SparseMultilinearExtension<F> {
        let mut s0 = F::zero();
        let mut s1 = F::zero();

        let mult_evaluations =self.fixed_mult.evaluations.iter();
        for (xy, val) in mult_evaluations {
            let dim = self.gkr_round.vi().num_vars;
            let x = xy & ((1 << dim) - 1);
            let y = xy >> dim;
            let value = *val * self.gkr_round.vi()[x] * self.gkr_round.vj()[y];
            if xy & 1 == 0 {
                s0 += value;
            } else {
                s1 += value;
            }
        }
        let add_evaluations =self.fixed_add.evaluations.iter();
        for (xy, val) in add_evaluations {
            let dim = self.gkr_round.vi().num_vars;
            let x = xy & ((1 << dim) - 1);
            let y = xy >> dim;
            let value = *val * self.gkr_round.vi()[x] + *val * self.gkr_round.vj()[y];
            if x & 1 == 0 {
                s0 += value;
            } else {
                s1 += value;
            }
        }
        SparseMultilinearExtension::from_evaluations(1, vec![&(0, s0), &(1, s1)])
    }

    fn fix_variable(&mut self, random_field_element: F) {
        /*
        1. Fix the first set of variables in the predicates for the gate.
        2. Then fix in vi. Once vi has no more variables fix vj.
        */
        let field_packed = &[random_field_element];
        self.fixed_variables.push(random_field_element);
        self.fixed_mult = self.fixed_mult.fix_variables(field_packed);
        self.fixed_add = self.fixed_add.fix_variables(field_packed);

        if self.gkr_round.vi.num_vars > 0 {
            self.gkr_round.vi = self.gkr_round.vi.fix_variables(field_packed);
        } else {
            self.gkr_round.vj = self.gkr_round.vj.fix_variables(field_packed);
        }
    }

    /// TODO refactor this possibly currently just copy+pasted from the fast implementation
    /// I could probably refacter them to be the same just not a priority rn.
    fn layer_reduction_message(&self, s_i_plus_1: usize) -> LayerReductionMessage<F> {
        let mle = &self.layer_value_mle;
        let b_star = self.fixed_variables[0..s_i_plus_1].to_vec();
        let c_star = self.fixed_variables[s_i_plus_1..2*s_i_plus_1].to_vec();
        let poly = FastProver::restrict_poly(&*b_star, &*c_star, mle);
        let z_1 = mle.evaluate(&b_star);
        let z_2 = mle.evaluate(&c_star);
        LayerReductionMessage::new(z_1, z_2, poly)
    }

    fn new(gkr_round: GKRRound<F>, gate: &[F]) -> Self {
        NaiveProver::new(gkr_round, gate)
    }
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Fr;
    use ark_ff::{One, Zero};
    use ark_poly::MultilinearExtension;
    use ark_std::{test_rng, UniformRand};
    use crate::gkr::gkr_round::GKRRound;
    use crate::provers::naive::NaiveProver;
    use crate::structures::data_structures::SumCheckProver;
    use crate::util;
    use crate::util::index_to_field_element;

    #[test]
    fn sanity_check() {
        for mask in 0..30 {
            let mut points = vec![Fr::zero(); 30];
            let field_index: Vec<Fr> = index_to_field_element(mask, 30);
            for j in 0..30 {
                let bit = (mask >> j) & 1 != 0;
                points[j] = if bit { Fr::one() } else { Fr::zero() }
            }
            assert_eq!(field_index, points)
        }
    }
    #[test]
    fn test_get_verifier_function() {
        let gkr_round = GKRRound::new_rand(7);
        let fixed_gate = util::random_gate(gkr_round.gate_labes());

        let mut prover: NaiveProver<Fr> = NaiveProver::new(gkr_round, &fixed_gate);
        // Now we test the g_func gives what we expect
        let verifier_func = prover.get_verifier_function();
        // now we evaluate this function at Fr::zero() and Fr::one() and it has to be equal to the sum it claims.
        // Just as the verifier would do.
        let verifier_sum = verifier_func.evaluations.iter().map(|(_, &v)| v).sum();
        let claimed_sum = prover.compute_sum();
        assert_eq!(claimed_sum, verifier_sum);
    }

    #[test]
    fn test_fixing_a_variable() {
        let mut rng = test_rng();
        let gkr_round = GKRRound::new_rand(7);
        let fixed_gate = util::random_gate(gkr_round.gate_labes());

        let mut prover: NaiveProver<Fr> = NaiveProver::new(gkr_round, &fixed_gate);

        let old_verifier_vars = prover.get_verifier_function().num_vars();

        let r_field = Fr::rand(&mut rng);
        prover.fix_variable(r_field);

        assert_eq!(prover.get_verifier_function().num_vars, old_verifier_vars);
        assert_eq!(old_verifier_vars, 1);
    }

    #[test]
    fn test_naive_is_same_as_arks_initially() {
        let gkr_round = GKRRound::new_rand(7);
        assert_eq!(gkr_round.mult_predicate().num_vars, 21);
        let fixed_gate = util::random_gate(7);
        assert_eq!(fixed_gate.len(), 7);

        let mut prover: NaiveProver<Fr> = NaiveProver::new(gkr_round.clone(), &fixed_gate.clone());
        let ark_sum = NaiveProver::ark_compute_sum_naive(&gkr_round.mult_predicate(), &gkr_round.add_predicate(), &gkr_round.vi, &gkr_round.vj, &*fixed_gate);

        let my_sum = prover.compute_sum();
        assert_eq!(my_sum, ark_sum);
    }

    #[test]
    fn test_fix_variable_value() {
        let gkr_round = GKRRound::new_rand(7);
        let fixed_gate = util::random_gate(gkr_round.gate_labes());

        let mut prover: NaiveProver<Fr> = NaiveProver::new(gkr_round, &fixed_gate);

        let verifier_func = prover.get_verifier_function();
        prover.fix_variable(Fr::zero());
        assert_eq!(prover.compute_sum(), verifier_func[0]);
    }
    #[test]
    fn test_naive_fix_same_as_fix_ark() {
        let gkr_round = GKRRound::new_rand(7);
        let fixed_gate = util::random_gate(gkr_round.gate_labes());

        let mut prover: NaiveProver<Fr> = NaiveProver::new(gkr_round.clone(), &fixed_gate);

        let r_field = Fr::rand(&mut test_rng());
        prover.fix_variable(r_field);

        let naive_sum = prover.compute_sum();
        let ark_sum = NaiveProver::ark_compute_sum_alr_fixed(
                                                            &gkr_round.mult_predicate(),
                                                            gkr_round.add_predicate(),
                                                            gkr_round.vi(),
                                                            gkr_round.vj(),
                                                            &fixed_gate,
                                                            &r_field,
        );
        assert_eq!(ark_sum, naive_sum);
    }
}