use std::arch::x86_64::_mm_aeskeygenassist_si128;
use ark_ff::{Field, Zero, One};
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, Polynomial, SparseMultilinearExtension};
use ark_std::iterable::Iterable;
use crate::circuit_structures::GateType;
use crate::data_structures::{GKRRound, LayerReductionMessage, SumCheckProver};
use crate::gkr_protocol::LayerReductionOracle;
use crate::util::index_to_field_element;

pub struct NaiveProver<F: Field> {
    gkr_round: GKRRound<F>,
    fixed_mult: SparseMultilinearExtension<F>,
    fix_variables: Vec<F>,
    layer_value_mle: DenseMultilinearExtension<F>,
}

impl<F: Field> NaiveProver<F> {
    pub fn new(
        gkr_round: GKRRound<F>,
        gate: &Vec<F>,
    ) -> NaiveProver<F> {
        //assert_eq!(gkr_round.mult().num_vars, gkr_round.vi().num_vars * 3);
        //assert_eq!(gkr_round.mult().num_vars, gkr_round.vj().num_vars * 3);
        //assert_eq!(gate.len(), gkr_round.gate_labes());
        let fixed_mult = gkr_round.mult().fix_variables(gate);
        NaiveProver {
            layer_value_mle: gkr_round.vi.clone(),
            gkr_round,
            fixed_mult,
            fix_variables: Vec::new(),
        }
    }

    fn get_mult_fixed(&self) -> SparseMultilinearExtension<F> {
        self.fixed_mult.clone()
    }

    fn line_point(b_star: &[F], c_star: &[F], t: F) -> Vec<F> {
        b_star
            .iter()
            .zip(c_star.iter())
            .map(|(b, c)| *b + t * (*c - *b))
        .collect()
    }

    fn poly_add_assign(lhs: &mut Vec<F>, rhs: &[F]) {
        if lhs.len() < rhs.len() {
            lhs.resize(rhs.len(), F::zero());
        }
        for (i, coeff) in rhs.iter().enumerate() {
            lhs[i] += *coeff;
        }
    }

    /// Multiply polynomial by (c0 + c1*x), coefficient form.
    fn poly_mul_linear(poly: &[F], c0: F, c1: F) -> Vec<F> {
        let mut out = vec![F::zero(); poly.len() + 1];
        for (i, a) in poly.iter().enumerate() {
            out[i] += *a * c0;
            out[i + 1] += *a * c1;
        }
        out
    }

    /// Interpolate polynomial coefficients from samples (x_i, y_i) using Lagrange basis.
    fn lagrange_interpolate_coeffs(xs: &[F], ys: &[F]) -> Vec<F> {
        assert_eq!(xs.len(), ys.len());
        let n = xs.len();
        let mut result = vec![F::zero(); n];

        for i in 0..n {
            // basis numerator: prod_{j != i} (x - x_j)
            let mut basis = vec![F::one()];
            let mut denom = F::one();

            for j in 0..n {
                if i == j {
                    continue;
                }
                basis = Self::poly_mul_linear(&basis, -xs[j], F::one());
                denom *= xs[i] - xs[j];
            }

            let inv = denom.inverse().expect("distinct interpolation points required");
            let scale = ys[i] * inv;
            for coeff in basis.iter_mut() {
                *coeff *= scale;
            }
            Self::poly_add_assign(&mut result, &basis);
        }

        result
    }

    /// Takes the variables needed to evalute the product of the GKR function.
    ///
    fn eval_g(&self, point: &Vec<F>) -> F {
        assert_eq!(point.len(), self.gkr_round.vi.num_vars + self.gkr_round.vj.num_vars);
        assert_eq!(self.get_mult_fixed().num_vars, self.gkr_round.vi.num_vars + self.gkr_round.vj.num_vars);
        /*
        A bit blank of ideas.
        */
        let u_len = self.gkr_round.vi.num_vars();
        let v_len = self.gkr_round.vj.num_vars();

        assert_eq!(point.len(), u_len + v_len);

        let u = &point[..u_len].to_vec();
        let v = &point[u_len..u_len + v_len].to_vec();

        let vi_val = self.gkr_round.vi().evaluate(u);
        let vj_val = self.gkr_round.vj.evaluate(v);
        let mult_val = self.get_mult_fixed().evaluate(point);

        match self.gkr_round.gate_type {
            GateType::Add => (vi_val + vj_val) * mult_val,
            GateType::Mul => vi_val * vj_val * mult_val,
        }
    }

    fn calculate_sum_naive(
        &self,
    ) -> F {
        let vi_variables = self.gkr_round.vi().num_vars;

        let f1_g = &self.fixed_mult;
        let mut sum_xy = F::zero();
        for x in 0..(1 << vi_variables) {
            let f2_x = self.gkr_round.vi()[x];
            let f1_gx = f1_g
                .fix_variables(&index_to_field_element(x, vi_variables))
                .to_dense_multilinear_extension();
            for y in 0..(1 << self.gkr_round.vj.num_vars) {
                match self.gkr_round.gate_type {
                    GateType::Add => sum_xy += f1_gx[y] * (f2_x + self.gkr_round.vj[y]),
                    GateType::Mul => sum_xy += f1_gx[y] * f2_x * self.gkr_round.vj[y],
                }
            }
        }
        sum_xy
    }

    pub fn ark_compute_sum_naive(
        f1: &SparseMultilinearExtension<F>,
        f2: &DenseMultilinearExtension<F>,
        f3: &DenseMultilinearExtension<F>,
        g: &[F],
        gate_type: &GateType
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
                match gate_type {
                    GateType::Add => sum_xy += f1_gx[y] * (f2_x + f3[y]),
                    GateType::Mul => sum_xy += f1_gx[y] * f2_x * f3[y],
                }
            }
        }
        sum_xy
    }

    pub fn compute_vi_sum(&mut self) -> F {
        //let len = self.fix_variables.len();
        //let fixed_vi_variables = &self.fix_variables[0..len];
        self.gkr_round.vi.evaluations.iter().sum()
    }

    pub fn compute_vj_sum(&mut self) -> F {
        self.gkr_round.vj.evaluations.iter().sum()
    }

    /// This function should be called in iteration i + 1 and will perform reduction to single point
    /// verification.
    pub fn final_get_verifier_func(&mut self) -> DenseMultilinearExtension<F> {
        let z_1 = self.compute_vi_sum();
        let z_2 = self.compute_vj_sum();
        DenseMultilinearExtension::from_evaluations_vec(1, vec![z_1, z_2])
    }

    pub fn eval_univariate(coeffs: &[F], x: F) -> F {
        coeffs.iter().rev().fold(F::zero(), |acc, c| acc * x + c)
    }
}

impl<F: Field> LayerReductionOracle<F> for NaiveProver<F> {
    fn layer_reduction_message(&mut self, b_star: &[F], c_star: &[F]) -> LayerReductionMessage<F> {
        assert_eq!(b_star.len(), c_star.len());

        // q(t) has degree <= k where k = number of variables in W_{i+1}
        let k = self.layer_value_mle.num_vars;
        assert_eq!(k, b_star.len(), "b*/c* dimension must match W_(i+1) arity");

        let mut xs = Vec::with_capacity(k + 1);
        let mut ys = Vec::with_capacity(k + 1);

        for i in 0..=k {
            let t_i = F::from(i as u64);
            let pt = Self::line_point(b_star, c_star, t_i);
            let y_i = self.layer_value_mle.evaluate(&pt);
            xs.push(t_i);
            ys.push(y_i);
        }

        let q_coeffs = Self::lagrange_interpolate_coeffs(&xs, &ys);

        let z1 = ys[0]; // q(0) = W(b*)
        let z2 = ys[1]; // q(1) = W(c*)

        LayerReductionMessage { z1, z2, q_coeffs }
    }
}

impl<F: Field> SumCheckProver<F> for NaiveProver<F> {
    // Needs to be refactored just my last sumcheck which i know works.
    fn compute_sum(&mut self) -> F {
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
        let n = self.gkr_round.vi().num_vars + self.gkr_round.vj.num_vars;

        let total = 1 << n;
        let mut s0 = F::zero();
        let mut s1 = F::zero();

        for mask in 0..total {
            let field_index: Vec<F> = index_to_field_element(mask, n);

            let value = self.eval_g(&field_index);
            if field_index[0].is_zero() {
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
        self.fix_variables.push(random_field_element);
        let field_packed = &[random_field_element];
        self.fixed_mult = self.fixed_mult.fix_variables(field_packed);

        if (self.gkr_round.vi.num_vars > 0) {
            self.gkr_round.vi = self.gkr_round.vi.fix_variables(field_packed);
        } else {
            self.gkr_round.vj = self.gkr_round.vj.fix_variables(field_packed);
        }
    }

    fn compute_z_1(&mut self) -> F {
        self.compute_vj_sum()
    }

    fn compute_z_2(&mut self) -> F {
        self.compute_vi_sum()
    }

}

mod tests {
    use ark_bls12_381::Fr;
    use ark_ff::{Field, Fp256, One, Zero};
    use ark_poly::{MultilinearExtension, Polynomial};
    use ark_std::{test_rng, UniformRand};
    use crate::data_structures::{GKRRound, SumCheckProver};
    use crate::gkr_protocol::LayerReductionOracle;
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
        let gkr_round = GKRRound::new_rand();
        let fixed_gate = util::random_gate(gkr_round.gate_labes());

        let mut prover: NaiveProver<Fr> = NaiveProver::new(gkr_round, &fixed_gate);
        // Now we test the g_func gives what we expect
        let verifier_func = prover.get_verifier_function();
        // now we evaluate this function at Fr::zero() and Fr::one() and it has to be equal to the sum it claims.
        // Just as the verifier would do.
        let verifier_sum = verifier_func.evaluations[0] + verifier_func.evaluations[1];
        let claimed_sum = prover.compute_sum();
        assert_eq!(claimed_sum, verifier_sum);
    }

    #[test]
    fn test_fixing_a_variable() {
        let mut rng = test_rng();
        let gkr_round = GKRRound::new_rand();
        let fixed_gate = util::random_gate(gkr_round.gate_labes());

        let mut prover: NaiveProver<Fr> = NaiveProver::new(gkr_round, &fixed_gate);

        let old_verifier_vars = prover.get_verifier_function().num_vars();

        let r_field = Fr::rand(&mut rng);
        prover.fix_variable(r_field);

        assert_eq!(prover.get_verifier_function().num_vars, old_verifier_vars);
    }

    #[test]
    fn test_naive_is_same_as_arks_initially() {
        let gkr_round = GKRRound::new_rand();
        let fixed_gate = util::random_gate(gkr_round.gate_labes());

        let mut prover: NaiveProver<Fr> = NaiveProver::new(gkr_round.clone(), &fixed_gate);
        let ark_sum = NaiveProver::ark_compute_sum_naive(&gkr_round.mult(), &gkr_round.vi, &gkr_round.vj, &*fixed_gate, &gkr_round.gate_type);

        let my_sum = prover.compute_sum();
        assert_eq!(my_sum, ark_sum);
    }

    #[test]
    fn test_fix_variable_value() {
        let gkr_round = GKRRound::new_rand();
        let fixed_gate = util::random_gate(gkr_round.gate_labes());

        let mut prover: NaiveProver<Fr> = NaiveProver::new(gkr_round, &fixed_gate);

        let verifier_func = prover.get_verifier_function();
        prover.fix_variable(Fr::zero());
        assert_eq!(prover.compute_sum(), verifier_func[0]);
    }

    #[test]
    fn test_layer_reduction_message_endpoints_and_random_point() {
        let mut rng = test_rng();
        let gkr_round = GKRRound::new_rand();
        let k = gkr_round.gate_labes();
        let fixed_gate = random_gate::<Fr>(k);
        let mut prover: NaiveProver<Fr> = NaiveProver::new(gkr_round.clone(), &fixed_gate);

        let b_star = (0..k).map(|_| Fr::rand(&mut rng)).collect::<Vec<_>>();
        let c_star = (0..k).map(|_| Fr::rand(&mut rng)).collect::<Vec<_>>();

        let msg = prover.layer_reduction_message(&b_star, &c_star);

        let q0 = NaiveProver::eval_univariate(&msg.q_coeffs, Fr::zero());
        let q1 = NaiveProver::eval_univariate(&msg.q_coeffs, Fr::one());

        assert_eq!(q0, msg.z1, "q(0) must equal z1");
        assert_eq!(q1, msg.z2, "q(1) must equal z2");

        let r = Fr::rand(&mut rng);
        let line_r = b_star
            .iter()
            .zip(c_star.iter())
            .map(|(b, c)| *b + r * (*c - *b))
            .collect::<Vec<_>>();

        let expected = gkr_round.vi.evaluate(&line_r);
        let got = NaiveProver::eval_univariate(&msg.q_coeffs, r);
        assert_eq!(got, expected, "q(r) must match W(b + r(c-b))");
    }
    #[test]
    fn test_layer_reduction_message_interpolation_points_and_degree_bound() {
        let mut rng = test_rng();
        let gkr_round = GKRRound::new_rand();
        let k = gkr_round.gate_labes();
        let fixed_gate = random_gate::<Fr>(k);
        let mut prover: NaiveProver<Fr> = NaiveProver::new(gkr_round.clone(), &fixed_gate);

        let b_star = (0..k).map(|_| Fr::rand(&mut rng)).collect::<Vec<_>>();
        let c_star = (0..k).map(|_| Fr::rand(&mut rng)).collect::<Vec<_>>();

        let msg = prover.layer_reduction_message(&b_star, &c_star);

        assert!(
            msg.q_coeffs.len() <= k + 1,
            "degree(q) must be <= k, so coeff length <= k+1"
        );

        for i in 0..=k {
            let t = Fr::from(i as u64);
            let line_t = b_star
                .iter()
                .zip(c_star.iter())
                .map(|(b, c)| *b + t * (*c - *b))
                .collect::<Vec<_>>();
            let expected = gkr_round.vi.evaluate(&line_t);
            let got = NaiveProver::eval_univariate(&msg.q_coeffs, t);
            assert_eq!(got, expected, "q(t_i) must match sampled line value");
        }
    }
}