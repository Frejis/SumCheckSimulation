use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, Polynomial, SparseMultilinearExtension};
use crate::gkr::gkr_round::GKRRound;
use crate::gkr::layer::LayerReductionMessage;
use crate::structures::circuit_structures::GateType;
use crate::structures::data_structures::SumCheckProver;
use crate::util::{index_to_field_element, interpolate_univariate, restrict_mle_to_line};

pub struct FastProver<F: Field> {
    fixed_mult: SparseMultilinearExtension<F>,
    fixed_add: SparseMultilinearExtension<F>,
    gkr_round: GKRRound<F>,
    mult_p: Vec<F>,
    mult_q: Vec<F>,
    add: Vec<F>,
    fixed_variables: Vec<F>,
    has_phase_two_been_init: bool,
    layer_value_mle: DenseMultilinearExtension<F>,
}

impl<F: Field> FastProver<F> {
    pub fn new(
        gkr_round: GKRRound<F>,
        gate: &[F],
    ) -> Self {
        if gkr_round.gate_type == GateType::Add {
            panic!();
        }
        let should_initialize_phase_one = gkr_round.vi.num_vars > 0;
        let fixed_mult = gkr_round.mult_predicate().fix_variables(gate);
        let fixed_add = gkr_round.add_predicate().fix_variables(gate);
        let mut temp_res = Self {
            fixed_mult,
            fixed_add,
            mult_p: Vec::new(),
            mult_q: Vec::new(),
            add: Vec::new(),
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
        match self.gkr_round.gate_type {
            GateType::Add => todo!(),
            GateType::Mul => self.init_phase_one_mult()
        }
    }

    fn initialize_phase_two(&mut self) {
        match self.gkr_round.gate_type {
            GateType::Add => { panic!() // I don't even want to bother at this point
                }
            GateType::Mul => self.init_phase_two_mult(),
        }
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
        self.add = vec![F::zero(); size];
        let add_predicate_nonzero = self.fixed_add.evaluations.iter();
        for (xy, value) in add_predicate_nonzero {
            let x = xy & ((1 << dim) - 1);
            let y = xy >> dim;
            let left_value = self.gkr_round.vi[x];
            let right_value = self.gkr_round.vj[y];
            self.add[x] += *value * left_value + *value * right_value;
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
        for i in 0..self.add.len() {
            sum += self.add[i]
        }
        sum
    }

    fn get_verifier_function(&mut self) -> SparseMultilinearExtension<F> {
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
        for mask in 0..self.add.len() {
            let mut value = self.add[mask];
            if mask & 1 == 0 { s0 += value; }
            else { s1 += value; }
        }

        SparseMultilinearExtension::from_evaluations(1, vec![&(0, s0), &(1, s1)])
    }

    fn fix_variable(&mut self, r: F) {
        self.fixed_variables.push(r);
        self.fix_variable_mult(r);
        self.fix_variable_add(r);
    }
    fn layer_reduction_message(&self, b_star: &[F], c_star: &[F]) -> LayerReductionMessage<F> {
        let k_ip1 = self.layer_value_mle.num_vars;
        assert_eq!(b_star.len(), k_ip1);
        assert_eq!(b_star.len(), c_star.len());

        let ts: Vec<F> = (0..=k_ip1).map(|i| F::from(i as u64)).collect();
        let values = restrict_mle_to_line(&self.layer_value_mle, b_star, c_star, &ts);
        let g = interpolate_univariate(&values, &ts);

        LayerReductionMessage::new(g.evaluate(&F::zero()), g.evaluate(&F::one()), g)
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
            // LSB bit = 0 hver gang.
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
        let n = self.add.len();
        let half = n >> 1;
        let mut new_add = Vec::with_capacity(half);

        for i in 0..half {
            // LSB bit = 0 hver gang.
            let index_0 = i << 1;
            let index_1 = index_0 | 1;
            let a_0 = self.add[index_0];
            let a_1 = self.add[index_1];
            new_add.push(a_0 + r * (a_1 - a_0));
        }

        self.add = new_add;
    }
}

#[cfg(test)]
mod test {
    use ark_bls12_381::Fr;
    use ark_ff::{One, Zero};
    use ark_poly::Polynomial;
    use ark_std::{test_rng, UniformRand};
    use crate::structures::circuit_structures::GateType;
    use crate::structures::data_structures::{SumCheckProver, SumCheckVerifier};
    use crate::gkr::gkr_round::GKRRound;
    use crate::provers::fast::FastProver;
    use crate::provers::naive::NaiveProver;
    use crate::util;
    use crate::util::random_gate;
    use crate::verifiers::standard_verifier::StandardVerifier;

    #[test]
    fn first_phase_sum_is_identical_to_ark() {
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand();
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
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand();
        let random_gate = random_gate(gkr_round.gate_labes());
        let mut fast_prover = FastProver::new(gkr_round.clone(), &random_gate);
        let mut naive_prover = NaiveProver::new(gkr_round, &random_gate);
        (fast_prover, naive_prover)
    }

    #[test]
    fn test_fix_variable_reduces_amount_of_variables() {
        let mut rand = test_rng();
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand();
        let random_gate = random_gate(gkr_round.gate_labes());
        let mut fast_prover = FastProver::new(gkr_round.clone(), &random_gate);

        let r_field = Fr::rand(&mut rand);

        fast_prover.fix_variable(r_field);
        assert_eq!(fast_prover.mult_p.len(), 1 << gkr_round.gate_labes() - 1);
        assert_eq!(fast_prover.mult_q.len(), 1 << gkr_round.gate_labes() - 1);
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
    fn test_fix_variable_value_zero() {
        let gkr_round = GKRRound::new_rand();
        let fixed_gate = util::random_gate(gkr_round.gate_labes());

        let mut prover: FastProver<Fr> = FastProver::new(gkr_round, &fixed_gate);

        let verifier_func = prover.get_verifier_function();
        prover.fix_variable(Fr::zero());
        assert_eq!(prover.compute_sum(), verifier_func[0]);
    }

    #[test]
    fn test_fix_variable_value_one() {
        let gkr_round = GKRRound::new_rand();
        let fixed_gate = util::random_gate(gkr_round.gate_labes());

        let mut prover: FastProver<Fr> = FastProver::new(gkr_round, &fixed_gate);

        let verifier_func = prover.get_verifier_function();
        prover.fix_variable(Fr::one());
        assert_eq!(prover.compute_sum(), verifier_func[1]);
    }

    #[test]
    fn test_verifier_functions_agree_with_naive() {
        let (mut fast, mut naive) = create_naive_and_fast_prover();
        let fast_verifier_func = fast.get_verifier_function();
        let naive_verifier_func = naive.get_verifier_function();

        assert_eq!(naive_verifier_func[1] + naive_verifier_func[0],
                   fast_verifier_func[1] + fast_verifier_func[0]
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

        assert_eq!(fast.compute_sum(), fast_verifier_func[1]);
        assert_eq!(naive.compute_sum(), naive_verifier_func[1]);
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
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand();
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
        let gkr_round = GKRRound::new_rand();
        let fixed_gate = util::random_gate(gkr_round.gate_labes());

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
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand();
        let random_gate = random_gate(gkr_round.gate_labes());
        let mut fast_prover = FastProver::new(gkr_round.clone(), &random_gate);
        // max degree is 5 for now, i think it should be 2.
        let verifier = StandardVerifier::new(5, fast_prover.compute_sum(), gkr_round);
        let verifier_func = fast_prover.get_verifier_function();
        assert!(verifier.check_claimed_value(&verifier_func));
    }

    #[test]
    fn test_naive_agrees_with_fast_forall_variables() {
        let mut gkr_round: GKRRound<Fr> = GKRRound::new_rand();
        gkr_round.gate_type = GateType::Mul;

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