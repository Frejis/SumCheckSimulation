use std::cmp::PartialEq;
use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, Polynomial, SparseMultilinearExtension};
use crate::circuit_structures::GateType;
use crate::data_structures::{GKRRound, Prover};
use crate::util::{index_to_field_element};

pub struct FastProver<F: Field> {
    fixed_mult: SparseMultilinearExtension<F>,
    gate: Vec<F>,
    gkr_round: GKRRound<F>,
    p: Vec<F>,
    q: Vec<F>,
    fixed_variables: Vec<F>,
    has_phase_two_been_init: bool,
}

impl<F: Field> FastProver<F> {
    pub fn new(
        gkr_round: GKRRound<F>,
        gate: &Vec<F>,
    ) -> Self {
        if gkr_round.gate_type == GateType::Add {
            panic!();
        }
        let should_initialize_phase_one = gkr_round.vi.num_vars > 0;
        let mut temp_res = Self {
            fixed_mult: gkr_round.mult().fix_variables(&*gate), // I don't think i needed to call clone.
            p: Vec::new(),
            q: Vec::new(),
            gate: gate.clone(),
            gkr_round: gkr_round.clone(),
            fixed_variables: Vec::new(),
            has_phase_two_been_init: false,
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
        match self.gkr_round.gate_type {
            GateType::Add => self.init_phase_one_add(),
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
        self.p = vec![F::zero(); size];
        self.q = vec![F::zero(); size];

        assert_eq!(self.fixed_variables.len(), self.gkr_round.vj.num_vars);
        let fixed_mult = self.gkr_round.mult().fix_variables(&self.fixed_variables);

        let fr = self.gkr_round.vi.evaluate(&self.fixed_variables);

        for i in 0..size {
            let field_index: Vec<F> = index_to_field_element(i, self.gkr_round.vj.num_vars);
            let combined_vec = Self::create_combined_vec_array(&self.fixed_variables, &field_index);
            self.p[i] = fixed_mult.evaluate(&combined_vec);
            self.q[i] = fr * self.gkr_round.vj.evaluate(&field_index);
        }
    }

    fn init_phase_one_add(&mut self) {
        let first_half = self.gkr_round.vi().num_vars;
        let second_half = self.gkr_round.vj().num_vars;
        let size = 1 << first_half;

        self.p = vec![F::zero(); size];
        self.q = vec![F::zero(); size];

        for i in 0..size {
            let i_index = index_to_field_element(i, first_half);
            let vi_val = self.gkr_round.vi().evaluate(&i_index);

            for j in 0..(1 << second_half) {
                let j_index = index_to_field_element(j, second_half);
                let vj_value = self.gkr_round.vj().evaluate(&j_index);

                let combined_vec = Self::create_combined_vec_array(&i_index, &j_index);
                let mult_val = self.fixed_mult.evaluate(&combined_vec);

                self.p[i] += mult_val * (vj_value + vi_val);
            }
        }
    }

    fn init_phase_one_mult(&mut self) {
        let first_half = self.gkr_round.vi().num_vars;
        let second_half = self.gkr_round.vj().num_vars;
        let size = 1 << first_half;

        self.p = vec![F::zero(); size];
        self.q = vec![F::zero(); size];

        for i in 0..size {
            let i_index = index_to_field_element(i, first_half);
            let vi_val = self.gkr_round.vi().evaluate(&i_index);
            self.q[i] = vi_val;

            for j in 0..(1 << second_half) {
                let j_index = index_to_field_element(j, second_half);
                let vj_value = self.gkr_round.vj().evaluate(&j_index);

                let combined_vec = Self::create_combined_vec_array(&i_index, &j_index);
                let mult_val = self.fixed_mult.evaluate(&combined_vec);

                self.p[i] += mult_val * vj_value;
            }
        }
    }
}

impl<F: Field> Prover<F> for FastProver<F> {
    fn compute_sum(&mut self) -> F { // This currently only works for the first half.
        if self.fixed_variables.len() - 1 == self.gkr_round.vi.num_vars() && !self.has_phase_two_been_init {
            // Now we have to initialize phase two.
            self.has_phase_two_been_init = true;
            self.initialize_phase_two();
        }
        let mut sum = F::zero();
        for i in 0..self.p.len() {
            match self.gkr_round.gate_type {
                GateType::Add => sum += self.p[i],
                GateType::Mul => sum += self.p[i] * self.q[i],
            }
        }
        sum
    }

    fn get_verifier_function(&self) -> DenseMultilinearExtension<F> {
        let mut s0 = F::zero();
        let mut s1 = F::zero();
        for mask in 0..self.p.len() {
            let value = match self.gkr_round.gate_type {
                GateType::Mul => self.p[mask] * self.q[mask],
                GateType::Add => self.p[mask],
            };
            if mask & 1 == 0 { s0 += value; }
            else { s1 += value; } }
        DenseMultilinearExtension::from_evaluations_vec(1, vec![s0, s1])
    }

    fn fix_variable(&mut self, r: F) {
        self.fixed_variables.push(r);
        match self.gkr_round.gate_type {
            GateType::Add => self.fix_variable_add(r), // This was more like a sanity check.
            GateType::Mul => self.fix_variable_mult(r),
        }
    }
}


impl<F: Field> FastProver<F> {
    fn fix_variable_add(&mut self, r: F) {
        let n = self.p.len();
        assert_eq!(n % 2, 0);
        let half = n >> 1;
        let mut new_p = Vec::with_capacity(half);
        let mut new_q = Vec::with_capacity(half);

        for i in 0..half {

            let index_0 = i << 1;
            let index_1 = index_0 | 1;
            let a_0 = self.p[index_0];
            let a_1 = self.p[index_1];

            new_p.push(a_0 + r * (a_1 - a_0));
            new_q.push(F::zero());
        }

        self.p = new_p;
        self.q = new_q;
    }
}

impl<F: Field> FastProver<F> {
    fn fix_variable_mult(&mut self, r: F) {
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
    use crate::circuit_structures::GateType;
    use crate::data_structures::{GKRRound, Prover, Verifier};
    use crate::fast_prover::FastProver;
    use crate::naive_sum_check::NaiveProver;
    use crate::standard_verifier::StandardVerifier;
    use crate::util::{random_gate, random_gkr_round_gates};

    #[test]
    fn first_phase_sum_is_identical_to_naive() {
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand();
        let random_gate = random_gate(gkr_round.gate_labes());
        let mut fast_prover = FastProver::new(gkr_round.clone(), &random_gate);

        let naive_sum = NaiveProver::ark_compute_sum_naive(&gkr_round.mult(), &gkr_round.vi, &gkr_round.vj, &random_gate, &gkr_round.gate_type);
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
    fn test_fix_variable_same_as_naive_mult_gate() {
        let mut rand = test_rng();
        let mut gkr_round: GKRRound<Fr> = GKRRound::new_rand();
        gkr_round.gate_type = GateType::Mul;
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
    fn test_fix_variable_same_as_naive_add_gate() {
        let mut rand = test_rng();
        let mut gkr_round: GKRRound<Fr> = GKRRound::new_rand();
        gkr_round.gate_type = GateType::Add;


        //let random_gate = random_gate(gkr_round.gate_labes());
        let random_gate = vec![Fr::zero(); gkr_round.gate_labes()];

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