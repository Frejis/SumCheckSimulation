//! The linear-time sum-check prover.

use ark_ff::Field;
use ark_poly::{
    DenseMultilinearExtension, MultilinearExtension, Polynomial, SparseMultilinearExtension,
};

use crate::sumcheck::{GKRRound, SumCheckProver};
use crate::util::index_to_field_element;

/// Linear-time sum-check prover following the prefix/suffix technique
/// (Thaler, ch. 4). Instead of re-evaluating the predicates every round it
/// maintains bookkeeping tables over the remaining hypercube:
///
/// - `mult_p[x] = Σ_y mult(g,x,y)·vj(y)` and `mult_q[x] = vi(x)`
/// - `add_pred[x] = Σ_y add(g,x,y)`, `add_f2[x] = vi(x)` and
///   `add_pred_f3[x] = Σ_y add(g,x,y)·vj(y)`
///
/// so the current sum is `Σ_x mult_p·mult_q + add_pred·add_f2 + add_pred_f3`.
/// Fixing a variable folds each table in half. Once all `x` variables are
/// fixed, phase two re-initializes the tables over the `y` variables.
pub struct FastProver<F: Field> {
    fixed_mult: SparseMultilinearExtension<F>,
    fixed_add: SparseMultilinearExtension<F>,
    round: GKRRound<F>,
    mult_p: Vec<F>,
    mult_q: Vec<F>,
    add_pred: Vec<F>,
    add_f2: Vec<F>,
    add_pred_f3: Vec<F>,
    fixed_variables: Vec<F>,
    has_phase_two_been_init: bool,
    /// Untouched copy of `vi`, kept for the layer-reduction message.
    layer_value_mle: DenseMultilinearExtension<F>,
}

impl<F: Field> FastProver<F> {
    pub fn new(round: GKRRound<F>, gate: &[F]) -> Self {
        let should_initialize_phase_one = round.vi.num_vars > 0;
        let fixed_mult = round.mult_predicate().fix_variables(gate);
        let fixed_add = round.add_predicate().fix_variables(gate);
        let mut prover = Self {
            fixed_mult,
            fixed_add,
            mult_p: Vec::new(),
            mult_q: Vec::new(),
            add_f2: Vec::new(),
            add_pred: Vec::new(),
            add_pred_f3: Vec::new(),
            fixed_variables: Vec::new(),
            has_phase_two_been_init: false,
            layer_value_mle: round.vi.clone(),
            round,
        };
        if should_initialize_phase_one {
            prover.initialize_phase_one();
        }
        prover
    }

    fn initialize_phase_one(&mut self) {
        let dim = self.round.vi().num_vars;
        let size = 1 << dim;

        self.mult_p = vec![F::zero(); size];
        self.mult_q = vec![F::zero(); size];
        self.add_f2 = vec![F::zero(); size];
        self.add_pred = vec![F::zero(); size];
        self.add_pred_f3 = vec![F::zero(); size];

        for (xy, value) in self.fixed_mult.evaluations.iter() {
            let x = xy & ((1 << dim) - 1);
            let y = xy >> dim;
            self.mult_p[x] += *value * self.round.vj[y];
        }
        for (xy, value) in self.fixed_add.evaluations.iter() {
            let x = xy & ((1 << dim) - 1);
            let y = xy >> dim;
            self.add_pred_f3[x] += *value * self.round.vj()[y];
            self.add_pred[x] += *value;
        }
        for i in 0..size {
            let i_index = index_to_field_element(i, dim);
            let vi_val = self.round.vi().evaluate(&i_index);
            self.mult_q[i] = vi_val;
            self.add_f2[i] = vi_val;
        }
    }

    /// Re-initializes the tables over the `y` variables once every `x`
    /// variable has been fixed at `b*`. The remaining sum is
    /// `Σ_y mult(b*,y)·vi(b*)·vj(y) + add(b*,y)·(vi(b*) + vj(y))`,
    /// where each factor is multilinear in `y`, so it maps onto the same
    /// two-table fold machinery as phase one (with `add_pred_f3` unused).
    fn initialize_phase_two(&mut self) {
        let size = 1 << self.round.vj.num_vars;
        self.mult_p = vec![F::zero(); size];
        self.mult_q = vec![F::zero(); size];
        self.add_pred = vec![F::zero(); size];
        self.add_f2 = vec![F::zero(); size];
        self.add_pred_f3 = vec![F::zero(); size];

        assert_eq!(self.fixed_variables.len(), self.round.vj.num_vars);
        let fixed_mult = self.fixed_mult.fix_variables(&self.fixed_variables);
        let fixed_add = self.fixed_add.fix_variables(&self.fixed_variables);
        let vi_at_fixed = self.round.vi.evaluate(&self.fixed_variables);
        for i in 0..size {
            let field_index: Vec<F> = index_to_field_element(i, self.round.vj.num_vars);
            let vj_i = self.round.vj.evaluate(&field_index);
            self.mult_q[i] = vi_at_fixed * vj_i;
            self.add_f2[i] = vi_at_fixed + vj_i;
        }
        for (j, value) in fixed_mult.evaluations {
            self.mult_p[j] += value;
        }
        for (j, value) in fixed_add.evaluations {
            self.add_pred[j] += value;
        }
    }

    /// Folds the multiplication tables in half by fixing their lowest
    /// variable to `r`.
    fn fix_variable_mult(&mut self, r: F) {
        let half = self.mult_p.len() >> 1;
        let mut new_p = Vec::with_capacity(half);
        let mut new_q = Vec::with_capacity(half);

        for i in 0..half {
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

    /// Folds the addition tables in half by fixing their lowest variable to `r`.
    fn fix_variable_add(&mut self, r: F) {
        let half = self.add_pred_f3.len() >> 1;
        let mut new_add_pred = Vec::with_capacity(half);
        let mut new_add_f2 = Vec::with_capacity(half);
        let mut new_add_pred_f3 = Vec::with_capacity(half);

        for i in 0..half {
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

impl<F: Field> SumCheckProver<F> for FastProver<F> {
    fn new(round: GKRRound<F>, gate: &[F]) -> Self {
        FastProver::new(round, gate)
    }

    fn compute_sum(&mut self) -> F {
        if self.fixed_variables.len() == self.round.vi.num_vars && !self.has_phase_two_been_init {
            self.has_phase_two_been_init = true;
            self.initialize_phase_two();
        }
        let mut sum = F::zero();
        for i in 0..self.mult_p.len() {
            sum += self.mult_p[i] * self.mult_q[i];
        }
        for i in 0..self.add_pred.len() {
            sum += self.add_pred[i] * self.add_f2[i];
            sum += self.add_pred_f3[i];
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
            } else {
                s1 += value;
            }
        }
        for mask in 0..self.add_pred.len() {
            let value = self.add_pred[mask] * self.add_f2[mask] + self.add_pred_f3[mask];
            if mask & 1 == 0 {
                s0 += value;
            } else {
                s1 += value;
            }
        }

        vec![(0, s0), (1, s1)]
    }

    fn fix_variable(&mut self, r: F) {
        self.fixed_variables.push(r);
        self.fix_variable_mult(r);
        self.fix_variable_add(r);
    }

    fn fixed_variables(&self) -> &[F] {
        &self.fixed_variables
    }

    fn layer_value_mle(&self) -> &DenseMultilinearExtension<F> {
        &self.layer_value_mle
    }
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Fr;
    use ark_std::{test_rng, UniformRand};

    use crate::sumcheck::{FastProver, GKRRound, SumCheckProver};
    use crate::util::random_gate;

    /// Needs access to the private bookkeeping tables, so it lives here
    /// rather than in the shared sum-check test module.
    #[test]
    fn fixing_a_variable_halves_the_tables() {
        let round: GKRRound<Fr> = GKRRound::new_rand(7);
        let gate = random_gate(round.gate_labels());
        let mut prover = FastProver::new(round.clone(), &gate);

        prover.fix_variable(Fr::rand(&mut test_rng()));

        assert_eq!(prover.mult_p.len(), 1 << (round.gate_labels() - 1));
        assert_eq!(prover.mult_q.len(), 1 << (round.gate_labels() - 1));
    }
}
