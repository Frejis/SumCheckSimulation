use std::arch::x86_64::_mm_aeskeygenassist_si128;
use ark_ff::{Field, SqrtPrecomputation, Zero};
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, Polynomial, SparseMultilinearExtension};
use ark_std::iterable::Iterable;
use crate::data_structures::Prover;
use crate::util::index_to_field_element;

pub struct NaiveProver<F: Field> {
    mult: SparseMultilinearExtension<F>, // mle of mult for round k with (r, i , j)
    fixed_mult: SparseMultilinearExtension<F>, // This has the gate r fixed.
    pub(crate) vi: DenseMultilinearExtension<F>, // mle of v_(k-1)(i)
    pub(crate) vj: DenseMultilinearExtension<F>, // mle of v_(k-1)(i)
    r: Vec<F>, // The gate "r" that is fixed.
}

impl<F: Field> NaiveProver<F> {
    pub fn new(
        mult: SparseMultilinearExtension<F>,
        vi: DenseMultilinearExtension<F>,
        vj: DenseMultilinearExtension<F>,
        r: Vec<F>,
    ) -> NaiveProver<F> {
        assert_eq!(mult.num_vars, vi.num_vars * 3);
        assert_eq!(mult.num_vars, vj.num_vars * 3);
        NaiveProver {
            fixed_mult: mult.clone().fix_variables(&*r),
            mult,
            vi,
            vj,
            r,
        }
    }

    fn get_mult_fixed(&self) -> SparseMultilinearExtension<F> {
        assert!(self.r.len() < self.mult.num_vars);
        self.fixed_mult.clone()
    }

    /// Takes the variables needed to evalute the product of the GKR function.
    ///
    fn eval_g(&self, point: &Vec<F>) -> F {
        assert_eq!(point.len(), self.vi.num_vars + self.vj.num_vars);
        assert_eq!(self.get_mult_fixed().num_vars, self.vi.num_vars + self.vj.num_vars);
        /*
        A bit blank of ideas.
        */
        let u_len = self.vi.num_vars();
        let v_len = self.vj.num_vars();

        assert_eq!(point.len(), u_len + v_len);

        let u = &point[..u_len].to_vec();
        let v = &point[u_len..u_len + v_len].to_vec();

        let vi_val = self.vi.evaluate(u);
        let vj_val = self.vj.evaluate(v);
        let mult_val = self.get_mult_fixed().evaluate(point);

        vi_val * vj_val * mult_val
    }

    fn calculate_sum_naive(
        &self,
    ) -> F {
        let vi_variables = self.vi.num_vars;

        let f1_g = &self.fixed_mult;
        let mut sum_xy = F::zero();
        for x in 0..(1 << vi_variables) {
            let f2_x = self.vi[x];
            let f1_gx = f1_g
                .fix_variables(&index_to_field_element(x, vi_variables))
                .to_dense_multilinear_extension();
            for y in 0..(1 << self.vj.num_vars) {
                sum_xy += f1_gx[y] * f2_x * self.vj[y];
            }
        }
        sum_xy
    }

    pub fn ark_compute_sum_naive(
        f1: &SparseMultilinearExtension<F>,
        f2: &DenseMultilinearExtension<F>,
        f3: &DenseMultilinearExtension<F>,
        g: &[F],
    ) -> F {
        let dim = f2.num_vars;
        let f1_g = f1.fix_variables(g);
        let mut sum_xy = F::zero();
        for x in 0..(1 << dim) {
            let f2_x = f2[x];
            let f1_gx = f1_g
                .fix_variables(&index_to_field_element(x, dim))
                .to_dense_multilinear_extension();
            for y in 0..(1 << dim) {
                sum_xy += f1_gx[y] * f2_x * f3[y];
            }
        }
        sum_xy
    }
}

impl<F: Field> Prover<F> for NaiveProver<F> {
    // Needs to be refactored just my last sumcheck which i know works.
    fn compute_sum(&self) -> F {
        self.calculate_sum_naive()
    }

    fn get_verifier_function(&self) -> DenseMultilinearExtension<F> {
        // clone existing functions.
        /*
        Iterate over all possible assignments of the bits.
        After that take the first variable and set it to 0 and the other one 1.
        */
        // Assume that the gate has been fixed.
        //assert_eq!(self.mult.num_vars, self.vi.num_vars * 3);
        let n = self.vi.num_vars + self.vj.num_vars;

        let total = 1 << n;
        let mut s0 = F::zero();
        let mut s1 = F::zero();

        for mask in 0..total {
            let field_index: Vec<F> = index_to_field_element(mask, n);

            let value = self.eval_g(&field_index);
            if (field_index[0].is_zero()) {
                s0 += &value;
            } else if (field_index[0].is_one()) {
                s1 += &value;
            }
        }

        DenseMultilinearExtension::from_evaluations_vec(1, vec![s0, s1])
    }

    fn fix_variable(&mut self, random_field_element: F) {
        /*
        1. Fix the first variable in mult. Then fix in vi. Once vi has no more variables fix vj.
        */
        let field_packed = &[random_field_element];
        self.fixed_mult = self.fixed_mult.fix_variables(field_packed);

        if (self.vi.num_vars > 0) {
            self.vi = self.vi.fix_variables(field_packed);
        } else {
            self.vj = self.vj.fix_variables(field_packed);
        }
    }
}

mod tests {
    use ark_bls12_381::Fr;
    use ark_ff::{One, Zero};
    use ark_std::{test_rng, UniformRand};
    use crate::data_structures::Prover;
    use crate::naive_sum_check::NaiveProver;
    use crate::util;
    use crate::util::{index_to_field_element, random_gate};

    #[test]
    fn sanity_check() {
        for mask in 0..30 {
            let mut points = vec![ark_bls12_381::Fr::zero(); 30];
            let field_index: Vec<ark_bls12_381::Fr> = index_to_field_element(mask, 30);
            for j in 0..30 {
                let bit = (mask >> j) & 1 != 0;
                points[j] = if bit { ark_bls12_381::Fr::one() } else { ark_bls12_381::Fr::zero() }
            }
            assert_eq!(field_index, points)
        }
    }
    #[test]
    fn test_get_verifier_function() {
        let mut rng = test_rng();

        // The first 5 variables are the gate, hence we also need 5 to fix the gate label.
        let variables = 7;
        let fixed_gate = util::random_gate(7);
        let (mult, vi, vj) = util::random_gkr_round_gates(variables, &mut rng);
        let prover: NaiveProver<Fr> = NaiveProver::new(mult, vi, vj, Vec::from(fixed_gate));
        // Now we test the g_func gives what we expect
        let verifier_func = prover.get_verifier_function();
        println!("{:?}", prover.mult.num_vars);
        println!("{:?}", prover.vi.num_vars);
        println!("{:?}", prover.vj.num_vars);
        // now we evaluate this function at Fr::zero() and Fr::one() and it has to be equal to the sum it claims.
        // Just as the verifier would do.
        let verifier_sum = verifier_func.evaluations[0] + verifier_func.evaluations[1];
        let claimed_sum = prover.compute_sum();
        println!("{:?}", claimed_sum);
        println!("{:?}", verifier_sum);
        assert_eq!(claimed_sum, verifier_sum);
    }

    #[test]
    fn test_fixing_a_variable() {
        let gate_length = 7;
        let fixed_gate: Vec<Fr> = random_gate(gate_length);

        let mut rng = test_rng();
        let (mult, vi, vj) = util::random_gkr_round_gates(gate_length, &mut rng);
        let mut prover = NaiveProver::new(mult, vi, vj, Vec::from(fixed_gate));

        let r_field = Fr::rand(&mut rng);
        prover.fix_variable(r_field);

        assert_eq!(prover.vi.num_vars + prover.vj.num_vars, gate_length * 2 - 1);
    }

    #[test]
    fn test_naive_is_same_as_arks_initially() {
        let gate_length = 7;
        let fixed_gate: Vec<Fr> = random_gate(gate_length);

        let mut rng = test_rng();
        let (mult, vi, vj) = util::random_gkr_round_gates(gate_length, &mut rng);
        let ark_sum = NaiveProver::ark_compute_sum_naive(&mult, &vi, &vj, &*fixed_gate);
        let mut prover = NaiveProver::new(mult, vi, vj, Vec::from(fixed_gate));

        let my_sum = prover.compute_sum();
        println!("My sum = {my_sum}, actual sum = {ark_sum}");
        assert_eq!(my_sum, ark_sum);
    }
}
