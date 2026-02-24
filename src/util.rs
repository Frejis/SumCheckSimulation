use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, SparseMultilinearExtension};
use ark_std::rand::RngCore;
use ark_std::test_rng;
use crate::data_structures::{GKRRound, Prover, Verifier};
use crate::fast_prover::FastProver;
use crate::naive_sum_check::NaiveProver;
use crate::standard_verifier::StandardVerifier;

/// Taken from arkworks sumcheck protocol.
/// Can be seen [here](https://github.com/arkworks-rs/sumcheck/blob/master/src/gkr_round_sumcheck/test.rs)
/// This should only really be used for testing.
pub fn random_gkr_round_gates<F: Field, R: RngCore>(
    dim: usize,
    rng: &mut R,
) -> (
    SparseMultilinearExtension<F>,
    DenseMultilinearExtension<F>,
    DenseMultilinearExtension<F>,
) {
    (
        SparseMultilinearExtension::rand_with_config(dim * 3, 1 << dim, rng),
        DenseMultilinearExtension::rand(dim, rng),
        DenseMultilinearExtension::rand(dim, rng),
    )
}
/// Also taken from [arkworks](https://github.com/arkworks-rs/sumcheck/blob/master/src/gkr_round_sumcheck/test.rs)
/// Takes as input ${0,1}^n$ and returns $\mathbb{F}^n$
pub fn index_to_field_element<F: Field>(mut index: usize, mut nv: usize) -> Vec<F> {
    let mut ans = Vec::with_capacity(nv);
    while nv != 0 {
        ans.push(((index & 1) as u64).into());
        index >>= 1;
        nv -= 1;
    }
    ans
}

/// Again this is a testing function since it relies on ark_std::test_rng()
pub fn random_gate<F: Field>(label_length: usize) -> Vec<F> {
    let mut rng = test_rng();
    /*
    let mut res = Vec::with_capacity(label_length);
    for i in 0..label_length {
        res[i] = F::rand(&mut rng);
    }
    res
    */
    let mut res = Vec::with_capacity(label_length);
    for _ in 0..label_length {
        res.push(F::rand(&mut rng));
    }
    res
}


struct GKRRoundSumCheckSimulator<F: Field, P: Prover<F>, V: Verifier<F>> {
    gkr_round: GKRRound<F>,
    prover: P,
    verifier: V,
}

impl<F: Field, P: Prover<F>, V: Verifier<F>> GKRRoundSumCheckSimulator<F, P, V> {

    pub fn new(prover: P, verifier: V) -> Self {
        Self {
            gkr_round: GKRRound::new_rand(),
            prover,
            verifier,
        }
    }

    pub fn new_fast_prover_std_verifier() -> Self {
        let gkr_round = GKRRound::new_rand();
        let random_gate = random_gate(gkr_round.gate_labes());

        let prover = FastProver::new(
            gkr_round.mult(),
            gkr_round.vi(),
            gkr_round.vj(),
            &random_gate,
        );

        let verifier = StandardVerifier::new(3, F::zero());

        Self {
            gkr_round,
            prover,
            verifier,
        }
    }

    pub fn new_naive_prover_std_verifier() -> Self {
        let gkr_round = GKRRound::new_rand();
        let random_gate = random_gate(gkr_round.gate_labes());
        Self {
            prover: &NaiveProver::new(*gkr_round.mult(), *gkr_round.vi(), *gkr_round.vj(), random_gate),
            verifier: &StandardVerifier::new(3, F::zero()),
            gkr_round,
        }
    }
}