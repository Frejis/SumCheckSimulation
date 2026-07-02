//! The sum-check protocol.
//!
//! Contains the prover/verifier traits ([`SumCheckProver`], [`SumCheckVerifier`]),
//! the per-layer protocol instance ([`GKRRound`]), a [`NaiveProver`] that
//! recomputes the sum from scratch every round, a linear-time [`FastProver`]
//! and the [`StandardVerifier`].

mod fast;
mod instance;
mod naive;
mod protocol;
mod verifier;

#[cfg(test)]
mod tests;

pub use fast::FastProver;
pub use instance::GKRRound;
pub use naive::NaiveProver;
pub use protocol::{restrict_poly, LayerReductionMessage, SumCheckProver, SumCheckVerifier};
pub use verifier::StandardVerifier;
