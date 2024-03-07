We estimate the survivor average causal effect (SACE) in cluster-randomized trials using a generalized linear mixed effects model approach.
The method extends the approach in the article "A mixed model approach to estimate the survivor average causal effect in 
cluster-randomized trials" which studies a continuous outcome. The current implemention focuses on a binary outcome variable. 
We model the outcome measure and the principal strata membership both with random intercepts to account for intracluster correlations 
due to cluster-level randomization. We test performance of the developed method through a simulation study. 

R functions

bme: Implementation of the EM algorithm to estimate model parameters.

mefit: SACE estimation for a binary outcome with random intercepts in both the membership models and the outcome models.

simfun: Functions used in the simulations.

simdata: Data simulation.

simfit: Evaluation of the performance in terms of the bias, the MSE, and the proportion of coverage.
