# Examples of the article

These six scripts are the ones presented in *TensorMixedStates: A Julia library for
simulating pure and mixed quantum states using matrix product states*,
[SciPost Phys. Codebases **72** (2026)](https://doi.org/10.21468/SciPostPhysCodeb.72).

They run at the sizes the article published — twenty to fifty sites, bond dimensions up to
six hundred, up to four hundred time steps — so they are long, and the continuous
integration does not run them. To try one out quickly, cut its sizes down: `MAXDIM` sits
near the top of each script, and the number of sites, the duration and the time step near
the end.

The shorter examples of [`../high_level`](../high_level) are the ones to read first to see
what the code looks like.
