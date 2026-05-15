use ark_ff::{Field};
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, Polynomial, SparseMultilinearExtension};
use crate::gkr::gkr_round::GKRRound;
use crate::gkr::layer::LayerReductionMessage;
use crate::structures::data_structures::SumCheckProver;
use crate::util::{index_to_field_element};

pub struct FastProver<F: Field> {
    fixed_mult: SparseMultilinearExtension<F>,
    fixed_add: SparseMultilinearExtension<F>,
    gkr_round: GKRRound<F>,
    mult_p: Vec<F>,
    mult_q: Vec<F>,
    add_pred: Vec<F>,
    add_f2: Vec<F>,
    add_pred_f3: Vec<F>,
    fixed_variables: Vec<F>,
    has_phase_two_been_init: bool,
    layer_value_mle: DenseMultilinearExtension<F>,
}

impl<F: Field> FastProver<F> {
    pub(crate) fn layer_reduction_message(&self, s_i_plus_1: usize) -> LayerReductionMessage<F> {
        let mle = &self.layer_value_mle;
        let b_star = self.fixed_variables[0..s_i_plus_1].to_vec();
        let c_star = self.fixed_variables[s_i_plus_1..2*s_i_plus_1].to_vec();
        let poly = FastProver::restrict_poly(&*b_star, &*c_star, mle);
        let z_1 = mle.evaluate(&b_star);
        let z_2 = mle.evaluate(&c_star);
        LayerReductionMessage::new(z_1, z_2, poly)
    }
}

impl<F: Field> FastProver<F> {
    pub fn new(
        gkr_round: GKRRound<F>,
        gate: &[F],
    ) -> Self {
        let should_initialize_phase_one = gkr_round.vi.num_vars > 0;
        let fixed_mult = gkr_round.mult_predicate().fix_variables(gate);
        let fixed_add = gkr_round.add_predicate().fix_variables(gate);
        let mut temp_res = Self {
            fixed_mult,
            fixed_add,
            mult_p: Vec::new(),
            mult_q: Vec::new(),
            add_f2: Vec::new(),
            add_pred: Vec::new(),
            add_pred_f3: Vec::new(),
            gkr_round: gkr_round.clone(),
            fixed_variables: Vec::new(),
            has_phase_two_been_init: false,
            layer_value_mle: gkr_round.vi,
        };
        if should_initialize_phase_one {
            temp_res.initialize_phase_one();
        }
        temp_res
    }

    fn initialize_phase_one(&mut self) {
        self.init_phase_one_mult()
    }

    fn initialize_phase_two(&mut self) {
        self.init_phase_two_mult()
    }

    fn init_phase_two_mult(&mut self) {
        let size = 1 << self.gkr_round.vj.num_vars;
        self.init_p_q_zero(size);
        assert_eq!(self.fixed_variables.len(), self.gkr_round.vj.num_vars);
        let fixed_mult = self.fixed_mult.fix_variables(&self.fixed_variables);
        let fr = self.gkr_round.vi.evaluate(&self.fixed_variables);
        for i in 0..size {
            let field_index: Vec<F> = index_to_field_element(i, self.gkr_round.vj.num_vars);
            self.mult_p[i] = fixed_mult.evaluate(&field_index);
            self.mult_q[i] = fr * self.gkr_round.vj.evaluate(&field_index);
        }
    }

    fn init_p_q_zero(&mut self, size: usize) {
        self.mult_p = vec![F::zero(); size];
        self.mult_q = vec![F::zero(); size];
    }

    fn init_phase_one_mult(&mut self) {
        let dim = self.gkr_round.vi().num_vars;
        let size = 1 << dim;
        self.init_phase_one_mult_arrays(dim, size);
        self.init_phase_one_add_arrays(dim, size);
    }

    fn init_phase_one_mult_arrays(&mut self, dim: usize, size: usize) {
        self.mult_p = vec![F::zero(); size];
        self.mult_q = vec![F::zero(); size];
        let mult_predicate_nonzero = self.fixed_mult.evaluations.iter();
        for (xy, value) in mult_predicate_nonzero {
            let x = xy & ((1 << dim) - 1);
            let y = xy >> dim;
            self.mult_p[x] += *value * self.gkr_round.vj[y];
        }
        for i in 0..size {
            let i_index = index_to_field_element(i, dim);
            let vi_val = self.gkr_round.vi().evaluate(&i_index);
            self.mult_q[i] = vi_val;
        }
    }

    fn init_phase_one_add_arrays(&mut self, dim: usize, size: usize) {
        self.add_pred = vec![F::zero(); size];
        self.add_f2 = vec![F::zero(); size];
        self.add_pred_f3 = vec![F::zero(); size];

        let add_predicate_nonzero = self.fixed_add.evaluations.iter();
        for (xy, value) in add_predicate_nonzero {
            let x = xy & ((1 << dim) - 1);
            let y = xy >> dim;
            self.add_pred_f3[x] += *value * self.gkr_round.vj()[y];
            self.add_pred[x] += *value;
        }

        for i in 0..size {
            let i_index = index_to_field_element(i, dim);
            let vi_val = self.gkr_round.vi().evaluate(&i_index);
            self.add_f2[i] = vi_val;
        }
    }
}

impl<F: Field> SumCheckProver<F> for FastProver<F> {
    fn compute_sum(&mut self) -> F { // This currently only works for the first half.
        if self.fixed_variables.len() == self.gkr_round.vi.num_vars() && !self.has_phase_two_been_init {
            // Now we have to initialize phase two.
            self.has_phase_two_been_init = true;
            self.initialize_phase_two();
        }
        let mut sum = F::zero();
        for i in 0..self.mult_p.len() {
            // Add the multiplication term
            sum += self.mult_p[i] * self.mult_q[i];
        }
        for i in 0..self.add_pred.len() {
            sum += self.add_pred[i] * self.add_f2[i];
            sum += self.add_pred_f3[i]
        }
        sum
    }

    fn get_verifier_function(&mut self) -> Vec<(usize, F)> {
        let mut s0 = F::zero();
        let mut s1 = F::zero();


        for mask in 0..self.mult_p.len() {
            let value = self.mult_p[mask] * self.mult_q[mask];
            if mask & 1 == 0 {
                s0 += value;
            }
            else {
                s1 += value;
            }
        }
        for mask in 0..self.add_pred.len() {
            let mut value = self.add_pred[mask] * self.add_f2[mask];
            value += self.add_pred_f3[mask];
            if mask & 1 == 0 { s0 += value; }
            else { s1 += value; }
        }

        vec![(0, s0), (1, s1)]
    }

    fn fix_variable(&mut self, r: F) {
        self.fixed_variables.push(r);
        self.fix_variable_mult(r);
        self.fix_variable_add(r);
    }

    fn layer_reduction_message(&self, s_i_plus_1: usize) -> LayerReductionMessage<F> {
        FastProver::layer_reduction_message(self, s_i_plus_1)
    }

    fn new(gkr_round: GKRRound<F>, gate: &[F]) -> Self {
        FastProver::new(gkr_round, gate)
    }
}

impl<F: Field> FastProver<F> {

    fn fix_variable_mult(&mut self, r: F) {
        let n = self.mult_p.len();
        //assert_eq!(n % 2, 0);
        let half = n >> 1;
        let mut new_p = Vec::with_capacity(half);
        let mut new_q = Vec::with_capacity(half);

        for i in 0..half {
            // LSB bit = 0 every time.
            let index_0 = i << 1;
            let index_1 = index_0 | 1;
            let a_0 = self.mult_p[index_0];
            let a_1 = self.mult_p[index_1];
            let b_0 = self.mult_q[index_0];
            let b_1 = self.mult_q[index_1];
            new_p.push(a_0 + r * (a_1 - a_0));
            new_q.push(b_0 + r * (b_1 - b_0));
        }

        self.mult_p = new_p;
        self.mult_q = new_q;
    }

    fn fix_variable_add(&mut self, r: F) {
        let n = self.add_pred_f3.len();
        let half = n >> 1;
        let mut new_add_pred = Vec::with_capacity(half);
        let mut new_add_f2 = Vec::with_capacity(half);
        let mut new_add_pred_f3 = Vec::with_capacity(half);

        for i in 0..half {
            // LSB bit = 0 every time.
            let index_0 = i << 1;
            let index_1 = index_0 | 1;
            let pred_0 = self.add_pred[index_0];
            let pred_1 = self.add_pred[index_1];
            let pred_f3_0 = self.add_pred_f3[index_0];
            let pred_f3_1 = self.add_pred_f3[index_1];
            let f2_0 = self.add_f2[index_0];
            let f2_1 = self.add_f2[index_1];
            new_add_pred.push(pred_0 + r * (pred_1 - pred_0));
            new_add_f2.push(f2_0 + r * (f2_1 - f2_0));
            new_add_pred_f3.push(pred_f3_0 + r * (pred_f3_1 - pred_f3_0));
        }

        self.add_pred = new_add_pred;
        self.add_f2 = new_add_f2;
        self.add_pred_f3 = new_add_pred_f3;
    }
}

#[cfg(test)]
mod test {
    use ark_bls12_381::Fr;
    use ark_ff::{One, Zero};
    use ark_poly::univariate::SparsePolynomial;
    use ark_std::{test_rng, UniformRand};
    use crate::structures::data_structures::{SumCheckProver, SumCheckVerifier};
    use crate::gkr::gkr_round::GKRRound;
    use crate::provers::fast::FastProver;
    use crate::provers::naive::NaiveProver;
    use crate::util::random_gate;
    use crate::verifiers::standard_verifier::StandardVerifier;

    #[test]
    fn first_phase_sum_is_identical_to_ark() {
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand(7);
        let random_gate = random_gate(gkr_round.gate_labes());
        let mut fast_prover = FastProver::new(gkr_round.clone(), &random_gate);

        let ark_sum = NaiveProver::ark_compute_sum_naive(&gkr_round.mult_predicate(), &gkr_round.add_predicate(), &gkr_round.vi, &gkr_round.vj, &random_gate);
        let fast_sum = fast_prover.compute_sum();
        assert_eq!(ark_sum, fast_sum);
    }

    #[test]
    fn first_phase_sum_is_identical_to_naive() {
        let (mut fast_prover, mut naive_prover) = create_naive_and_fast_prover();
        assert_eq!(naive_prover.compute_sum(), fast_prover.compute_sum());
    }

    fn create_naive_and_fast_prover() -> (FastProver<Fr>, NaiveProver<Fr>) {
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand(7);
        let random_gate = random_gate(gkr_round.gate_labes());
        let fast_prover = FastProver::new(gkr_round.clone(), &random_gate);
        let naive_prover = NaiveProver::new(gkr_round, &random_gate);
        (fast_prover, naive_prover)
    }

    #[test]
    fn test_fix_variable_reduces_amount_of_variables() {
        let mut rand = test_rng();
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand(7);
        let random_gate = random_gate(gkr_round.gate_labes());
        let mut fast_prover = FastProver::new(gkr_round.clone(), &random_gate);

        let r_field = Fr::rand(&mut rand);

        fast_prover.fix_variable(r_field);
        assert_eq!(fast_prover.mult_p.len(), 1 << gkr_round.gate_labes() - 1);
        assert_eq!(fast_prover.mult_q.len(), 1 << gkr_round.gate_labes() - 1);
    }

    #[test]
    fn test_verifier_func_value() {
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand(7);
        let random_gate = random_gate(gkr_round.gate_labes());
        let mut prover = FastProver::new(gkr_round.clone(), &random_gate);

        let verifier_func = prover.get_verifier_function();
        prover.fix_variable(Fr::zero());
        assert_eq!(prover.compute_sum(), verifier_func[0].1);
    }

    #[test]
    fn test_fix_variable_value_zero() {
        let gkr_round = GKRRound::new_rand(7);
        let fixed_gate = random_gate(gkr_round.gate_labes());

        let mut prover: FastProver<Fr> = FastProver::new(gkr_round, &fixed_gate);

        let verifier_func = prover.get_verifier_function();
        prover.fix_variable(Fr::zero());
        assert_eq!(prover.compute_sum(), verifier_func[0].1);
    }

    #[test]
    fn test_fix_variable_value_one() {
        let gkr_round = GKRRound::new_rand(7);
        let fixed_gate = random_gate(gkr_round.gate_labes());

        let mut prover: FastProver<Fr> = FastProver::new(gkr_round, &fixed_gate);

        let verifier_func = prover.get_verifier_function();
        prover.fix_variable(Fr::one());
        assert_eq!(prover.compute_sum(), verifier_func[1].1);
    }

    #[test]
    fn test_verifier_functions_agree_with_naive() {
        let (mut fast, mut naive) = create_naive_and_fast_prover();
        let fast_verifier_func = fast.get_verifier_function();
        let naive_verifier_func = naive.get_verifier_function();

        assert_eq!(naive_verifier_func[1].1 + naive_verifier_func[0].1,
                   fast_verifier_func[1].1 + fast_verifier_func[0].1
        );
        assert_eq!(naive_verifier_func[1], fast_verifier_func[1]);
        assert_eq!(naive_verifier_func[0], fast_verifier_func[0]);
    }

    #[test]
    fn test_verifier_functions_agree_after_fixing_fixed_value() {
        let (mut fast, mut naive) = create_naive_and_fast_prover();
        let fast_verifier_func = fast.get_verifier_function();
        let naive_verifier_func = naive.get_verifier_function();

        fast.fix_variable(Fr::one());
        naive.fix_variable(Fr::one());

        assert_eq!(fast.compute_sum(), fast_verifier_func[1].1);
        assert_eq!(naive.compute_sum(), naive_verifier_func[1].1);
        assert_eq!(naive_verifier_func[1], fast_verifier_func[1]);
    }

    #[test]
    fn test_verifier_functions_agree_after_fixing_fixed_random() {
        let (mut fast, mut naive) = create_naive_and_fast_prover();

        //let random = Fr::rand(&mut test_rng());
        let random = Fr::one() + Fr::one();

        fast.fix_variable(random);
        naive.fix_variable(random);
        assert_eq!(fast.compute_sum(), naive.compute_sum());
    }

    #[test]
    fn test_fix_variable_same_sum_as_naive() {
        let mut rand = test_rng();
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand(7);
        let r_field = Fr::rand(&mut rand);
        //let r_field = Fr::one();

        let random_gate = random_gate(gkr_round.gate_labes());

        let mut fast_prover = FastProver::new(gkr_round.clone(), &random_gate);
        let mut naive_prover = NaiveProver::new(gkr_round, &random_gate);

        assert_eq!(fast_prover.compute_sum(), naive_prover.compute_sum());

        fast_prover.fix_variable(r_field);
        naive_prover.fix_variable(r_field);

        assert_eq!(fast_prover.compute_sum(), naive_prover.compute_sum());
    }

    #[test]
    fn test_fast_fix_same_as_fix_ark() {
        let gkr_round = GKRRound::new_rand(7);
        let fixed_gate = random_gate(gkr_round.gate_labes());

        let mut prover: FastProver<Fr> = FastProver::new(gkr_round.clone(), &fixed_gate);

        let r_field = Fr::rand(&mut test_rng());
        prover.fix_variable(r_field);

        let fast_sum = prover.compute_sum();
        let ark_sum = NaiveProver::ark_compute_sum_alr_fixed(
            &gkr_round.mult_predicate(),
            gkr_round.add_predicate(),
            gkr_round.vi(),
            gkr_round.vj(),
            &fixed_gate,
            &r_field,
        );
        assert_eq!(ark_sum, fast_sum);
    }

    #[test]
    fn test_get_verifier_function() {
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand(7);
        let random_gate = random_gate(gkr_round.gate_labes());
        let mut fast_prover = FastProver::new(gkr_round.clone(), &random_gate);
        // max degree is 5 for now, I think it should be 2.
        let verifier = StandardVerifier::new(2, fast_prover.compute_sum());
        let points = fast_prover.get_verifier_function();
        let x = verifier.check_claimed_value(&SparsePolynomial::from_coefficients_vec(points));
        assert!(x);
    }

    #[test]
    fn test_naive_agrees_with_fast_forall_variables() {
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand(7);

        let r_gate = random_gate(gkr_round.gate_labes());
        let mut fast_prover = FastProver::new(gkr_round.clone(), &r_gate);
        let mut naive_prover = NaiveProver::new(gkr_round.clone(), &r_gate);
        for i in 0..gkr_round.gate_labes() {
            println!("Round number: {i}");
            let r_variable = Fr::rand(&mut test_rng());
            fast_prover.fix_variable(r_variable);
            naive_prover.fix_variable(r_variable);
            assert_eq!(fast_prover.compute_sum(), naive_prover.compute_sum());
        }
        assert_eq!(fast_prover.fixed_variables.len(), gkr_round.gate_labes());
    }
}