//! A sum-check prover that recomputes the full sum from scratch every round.

use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, SparseMultilinearExtension};

use crate::sumcheck::{GKRRound, SumCheckProver};
use crate::util::index_to_field_element;

/// Naive sum-check prover: fixes variables directly in the multilinear
/// extensions and recomputes the whole sum whenever it is queried, giving the
/// baseline the linear-time [`FastProver`](crate::sumcheck::FastProver) is
/// benchmarked against.
pub struct NaiveProver<F: Field> {
    round: GKRRound<F>,
    /// The multiplication predicate with the gate label already fixed.
    fixed_mult: SparseMultilinearExtension<F>,
    /// The addition predicate with the gate label already fixed.
    fixed_add: SparseMultilinearExtension<F>,
    fixed_variables: Vec<F>,
    /// Untouched copy of `vi`, kept for the layer-reduction message.
    layer_value_mle: DenseMultilinearExtension<F>,
}

impl<F: Field> NaiveProver<F> {
    pub fn new(round: GKRRound<F>, gate: &[F]) -> Self {
        let fixed_mult = round.mult_predicate().fix_variables(gate);
        let fixed_add = round.add_predicate().fix_variables(gate);
        Self {
            layer_value_mle: round.vi.clone(),
            round,
            fixed_mult,
            fixed_add,
            fixed_variables: Vec::new(),
        }
    }

    /// Evaluates the GKR layer sum
    /// `Σ_{x,y} add(g,x,y)·(vi(x) + vj(y)) + mult(g,x,y)·vi(x)·vj(y)`
    /// over the remaining hypercube by brute force.
    fn calculate_sum_naive(&self) -> F {
        let vi_variables = self.round.vi().num_vars;
        let mut sum_xy = F::zero();
        for x in 0..(1 << vi_variables) {
            let vi_x = self.round.vi()[x];
            let x_index = index_to_field_element(x, vi_variables);
            let mult_at_x = self
                .fixed_mult
                .fix_variables(&x_index)
                .to_dense_multilinear_extension();
            let add_at_x = self
                .fixed_add
                .fix_variables(&x_index)
                .to_dense_multilinear_extension();
            for y in 0..(1 << self.round.vj.num_vars) {
                // ~add_i(g,b,c) · ~W_i(b)
                sum_xy += add_at_x[y] * vi_x;
                // ~add_i(g,b,c) · ~W_i(c)
                sum_xy += add_at_x[y] * self.round.vj[y];
                // ~mult_i(g,b,c) · ~W_i(b) · ~W_i(c)
                sum_xy += mult_at_x[y] * vi_x * self.round.vj[y];
            }
        }
        sum_xy
    }
}

impl<F: Field> SumCheckProver<F> for NaiveProver<F> {
    fn new(round: GKRRound<F>, gate: &[F]) -> Self {
        NaiveProver::new(round, gate)
    }

    fn compute_sum(&mut self) -> F {
        self.calculate_sum_naive()
    }

    fn get_verifier_function(&mut self) -> Vec<(usize, F)> {
        let mut s0 = F::zero();
        let mut s1 = F::zero();
        let dim = self.round.vi().num_vars;

        for (xy, val) in self.fixed_mult.evaluations.iter() {
            let x = xy & ((1 << dim) - 1);
            let y = xy >> dim;
            let value = *val * self.round.vi()[x] * self.round.vj()[y];
            if xy & 1 == 0 {
                s0 += value;
            } else {
                s1 += value;
            }
        }
        for (xy, val) in self.fixed_add.evaluations.iter() {
            let x = xy & ((1 << dim) - 1);
            let y = xy >> dim;
            let value = *val * self.round.vi()[x] + *val * self.round.vj()[y];
            if xy & 1 == 0 {
                s0 += value;
            } else {
                s1 += value;
            }
        }
        vec![(0, s0), (1, s1)]
    }

    fn fix_variable(&mut self, r: F) {
        // 1. Fix the next variable in the predicates.
        // 2. Then fix it in vi; once vi has no more variables, fix vj.
        let packed = &[r];
        self.fixed_variables.push(r);
        self.fixed_mult = self.fixed_mult.fix_variables(packed);
        self.fixed_add = self.fixed_add.fix_variables(packed);

        if self.round.vi.num_vars > 0 {
            self.round.vi = self.round.vi.fix_variables(packed);
        } else {
            self.round.vj = self.round.vj.fix_variables(packed);
        }
    }

    fn fixed_variables(&self) -> &[F] {
        &self.fixed_variables
    }

    fn layer_value_mle(&self) -> &DenseMultilinearExtension<F> {
        &self.layer_value_mle
    }
}
