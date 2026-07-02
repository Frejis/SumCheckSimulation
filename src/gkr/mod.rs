//! The GKR protocol over layered arithmetic circuits.
//!
//! [`GKRCircuit`] describes the wiring, [`compute_predicates`] turns it into
//! per-layer add/mult predicates, and [`GKRDriver`] runs the full protocol —
//! one sum-check per layer — between a [`GKRProver`] and a [`GKRVerifier`].

mod circuit;
mod driver;
mod predicates;
mod prover;
mod verifier;

pub use circuit::{
    EvaluatedGKRCircuit, EvaluatedLayer, GKRCircuit, Gate, GateType, InputLayer, Layer,
};
pub use driver::{GKRDriver, LayerConnection};
pub use predicates::{compute_predicates, LayerPredicates};
pub use prover::GKRProver;
pub use verifier::GKRVerifier;
