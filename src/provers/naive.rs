use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, Polynomial, SparseMultilinearExtension};
use crate::gkr::gkr_round::GKRRound;
use crate::gkr::layer::LayerReductionMessage;
use crate::structures::circuit_structures::GateType;
use crate::structures::data_structures::SumCheckProver;
use crate::util::{_line_point, index_to_field_element};

pub struct NaiveProver<F: Field> {
    gkr_round: GKRRound<F>,
    fixed_mult: SparseMultilinearExtension<F>,
    // Keep original next-layer values for layer reduction
    original_next_layer: DenseMultilinearExtension<F>,
}

impl<F: Field> NaiveProver<F> {
    pub fn new(
        gkr_round: GKRRound<F>,
        gate: &Vec<F>,
    ) -> NaiveProver<F> {
        let fixed_mult = gkr_round.mult().fix_variables(gate);

        let original_next_layer = gkr_round.vi.clone(); // Store original for layer reduction
        NaiveProver {
            gkr_round,
            fixed_mult,
            original_next_layer,
        }
    }

    fn get_mult_fixed(&self) -> SparseMultilinearExtension<F> {
        self.fixed_mult.clone()
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
        let n = self.gkr_round.vi().num_vars + self.gkr_round.vj.num_vars;
        let mut sum = F::zero();
        for i in 0..(1 << n) {
            let point = index_to_field_element(i, n);
            sum += self.eval_g(&point);
        }
        sum
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
}

impl<F: Field> SumCheckProver<F> for NaiveProver<F> {
    // Needs to be refactored just my last sumcheck which I know works.
    fn compute_sum(&mut self) -> F {
        self.calculate_sum_naive()
    }

    fn get_verifier_function(&self) -> Vec<F> {
        let n = self.gkr_round.vi().num_vars + self.gkr_round.vj.num_vars;
        if n == 0 {
            return vec![F::zero(), F::zero(), F::zero()];
        }

        let total = 1 << (n - 1);
        let mut s0 = F::zero();
        let mut s1 = F::zero();
        let mut s2 = F::zero();

        for mask in 0..total {
            let rest_vars = index_to_field_element(mask, n - 1);
            
            // Evaluate g(0, rest_vars), g(1, rest_vars), g(2, rest_vars)
            let mut p0 = vec![F::zero()];
            p0.extend(&rest_vars);
            s0 += self.eval_g(&p0);
            
            let mut p1 = vec![F::one()];
            p1.extend(&rest_vars);
            s1 += self.eval_g(&p1);
            
            let mut p2 = vec![F::from(2u64)];
            p2.extend(&rest_vars);
            s2 += self.eval_g(&p2);
        }

        vec![s0, s1, s2]
    }

    fn fix_variable(&mut self, random_field_element: F) {
        let field_packed = &[random_field_element];
        // Arkworks' fix_variables(point) fixes variables starting from index 0.
        // If we want it to fix variables in the order we use them in get_verifier_function (which is at the FRONT),
        // then fix_variables([r]) will fix index 0, which is exactly what we want.
        
        self.fixed_mult = self.fixed_mult.fix_variables(field_packed);

        if self.gkr_round.vi.num_vars > 0 {
            self.gkr_round.vi = self.gkr_round.vi.fix_variables(field_packed);
        } else if self.gkr_round.vj.num_vars > 0 {
            self.gkr_round.vj = self.gkr_round.vj.fix_variables(field_packed);
        }
    }

    fn layer_reduction_message(&self, b_star: &[F], c_star: &[F]) -> LayerReductionMessage<F> {
        // Use the original, unmodified next-layer values
        let k_ip1 = self.original_next_layer.num_vars;
        assert_eq!(b_star.len(), k_ip1);
        assert_eq!(b_star.len(), c_star.len());

        let z1 = self.original_next_layer.evaluate(&b_star.to_vec());
        let z2 = self.original_next_layer.evaluate(&c_star.to_vec());

        // Build polynomial q(t) = W(line(t)) where line(t) = (1-t)*b* + t*c*
        // W is multilinear, so W(line(t)) is degree at most k_ip1.
        let mut evaluations = Vec::with_capacity(k_ip1 + 1);
        for i in 0..=k_ip1 {
            let t = F::from(i as u64);
            let point = _line_point(b_star, c_star, t);
            evaluations.push(self.original_next_layer.evaluate(&point));
        }

        LayerReductionMessage::new(z1, z2, evaluations)
    }
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Fr;
    use ark_ff::{One, Zero};
    use ark_poly::Polynomial;
    use ark_std::{test_rng, UniformRand};
    use crate::gkr::gkr_round::GKRRound;
    use crate::structures::data_structures::SumCheckProver;
    use crate::provers::naive::NaiveProver;
    use crate::util;
    use crate::util::index_to_field_element;

    #[test]
    fn test_layer_reduction_message() {
        let gkr_round = GKRRound::new_rand();
        let fixed_gate = util::random_gate(gkr_round.gate_labes());
        let b_star = util::random_gate(gkr_round.vi().num_vars);
        let c_star = util::random_gate(gkr_round.vi().num_vars);

        let prover: NaiveProver<Fr> = NaiveProver::new(gkr_round, &fixed_gate);
        let message = prover.layer_reduction_message(&*b_star, &*c_star);
        
        assert_eq!(message.z1(), prover.original_next_layer.evaluate(&b_star));
        assert_eq!(message.z2(), prover.original_next_layer.evaluate(&c_star));
    }

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
        let gkr_round = GKRRound::new_rand();
        let fixed_gate = util::random_gate(gkr_round.gate_labes());

        let prover: NaiveProver<Fr> = NaiveProver::new(gkr_round, &fixed_gate);
        // Now we test the g_func gives what we expect
        let verifier_func = prover.get_verifier_function();
        // now we evaluate this function at Fr::zero() and Fr::one() and it has to be equal to the sum it claims.
        // Just as the verifier would do.
        let verifier_sum = verifier_func[0] + verifier_func[1];
        let claimed_sum = prover.calculate_sum_naive();
        assert_eq!(claimed_sum, verifier_sum);
    }

    #[test]
    fn test_fixing_a_variable() {
        let mut rng = test_rng();
        let gkr_round = GKRRound::new_rand();
        let fixed_gate = util::random_gate(gkr_round.gate_labes());

        let mut prover: NaiveProver<Fr> = NaiveProver::new(gkr_round, &fixed_gate);

        let old_verifier_vars = prover.gkr_round.vi.num_vars + prover.gkr_round.vj.num_vars;

        let r_field = Fr::rand(&mut rng);
        prover.fix_variable(r_field);

        let new_verifier_vars = prover.gkr_round.vi.num_vars + prover.gkr_round.vj.num_vars;
        assert_eq!(new_verifier_vars, old_verifier_vars - 1);
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
        let q0 = verifier_func[0];
        prover.fix_variable(Fr::zero());
        assert_eq!(prover.calculate_sum_naive(), q0);
    }
}