# SumCheckSimulation
This is a simulation of the sum-check protocol. It will naturally follow my own implementation of the naive implementation compared with the prefix suffix sum-check protocol.

Initially this was meant as  a pure comparison between my own implementation to see how the difference would be between 
a fast Prover using the ideas uncovered in Libra versus the "Naive" sum-check prover.

To realize (realise?) this i implemented a GKR protocol. I had originally thought to compare my speed to the 
implementation done in Arkworks (the framework i am using).
However upon further investigation it seems their implementation only supports multiplication predicates.
As such when generating circuits that have consists of both addition and multiplication gates the Arkworks impl falls
short. To combat this however I modified their naive computation of the sum to handle addition as well.

It has to be said that their simulation is difficult to read so I may not have a full grasp on what they are doing so
the following will be my impression of it in an informal debate sense.

I am unsure whether they implemented Shamir? Something to make the protocol non-interactive. Also they directly seem
to have implemented the Libra functions which I only realized had pseudo algorithms after implementing and setting 
up my own strucutre, yet the pseudo algorithms do not describe fully what is expected to be implemented either imo.

Anyways, the support for addition gate in the fast Sum-check prover stems from the description done in Libra.
I tried doing a smarter way of it (only using one vector and other combinations) described further in the report as a
remark (if i remember to put it in :pray:).

Currently the plan is what I would like to add. The idea is to make it easier to illustrate what the prover is actually
doing. Nothing beats a figure with good old descriptions, however atm i am unsure how to generate such a figure.

# Plan
1. Automatic figures generated that illustrates the benchmark so it is easy to compare.
2. Create a figure of the circuits that the benchmarks are running on (generate a figure from a gkr circuit struct).
