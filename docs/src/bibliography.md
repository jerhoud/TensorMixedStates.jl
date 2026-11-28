# Bibliography

The articles the methods of TMS come from, each with what it gives TMS. To cite TMS itself, see
[References](index.md#References).

## Matrix product states and operators

- M. Fishman, S. R. White and E. M. Stoudenmire, *The ITensor Software Library for Tensor
  Network Calculations*, SciPost Phys. Codebases 4 (2022),
  [doi:10.21468/SciPostPhysCodeb.4](https://doi.org/10.21468/SciPostPhysCodeb.4),
  [arXiv:2007.14822](https://arxiv.org/abs/2007.14822): the tensors, the MPS and MPO and the
  solvers TMS is built on.
- S. R. White, *Density matrix formulation for quantum renormalization groups*, Phys. Rev.
  Lett. 69, 2863 (1992),
  [doi:10.1103/PhysRevLett.69.2863](https://doi.org/10.1103/PhysRevLett.69.2863): DMRG, which
  `Dmrg` runs through ITensorMPS.
- C. Hubig, I. P. McCulloch and U. Schollwöck, *Generic construction of efficient matrix
  product operators*, Phys. Rev. B 95, 035129 (2017),
  [doi:10.1103/PhysRevB.95.035129](https://doi.org/10.1103/PhysRevB.95.035129),
  [arXiv:1611.02498](https://arxiv.org/abs/1611.02498): the terms of an operator compacted into
  as few channels as their rank, the sums and products of `PreMPO`, the factors rescaled by
  powers of two.
- D. E. Parker, X. Cao and M. P. Zaletel, *Local matrix product operators: Canonical form,
  compression, and control theory*, Phys. Rev. B 102, 035147 (2020),
  [doi:10.1103/PhysRevB.102.035147](https://doi.org/10.1103/PhysRevB.102.035147),
  [arXiv:1909.06341](https://arxiv.org/abs/1909.06341): the truncation of an MPO by
  `mpo_limits`, on an orthonormal basis of each site, the identity apart.

## Mixed states

- F. Verstraete, J. J. García-Ripoll and J. I. Cirac, *Matrix Product Density Operators:
  Simulation of Finite-Temperature and Dissipative Systems*, Phys. Rev. Lett. 93, 207204
  (2004), [doi:10.1103/PhysRevLett.93.207204](https://doi.org/10.1103/PhysRevLett.93.207204),
  [arXiv:cond-mat/0406426](https://arxiv.org/abs/cond-mat/0406426), and M. Zwolak and G. Vidal,
  *Mixed-State Dynamics in One-Dimensional Quantum Lattice Systems: A Time-Dependent
  Superoperator Renormalization Algorithm*, Phys. Rev. Lett. 93, 207205 (2004),
  [doi:10.1103/PhysRevLett.93.207205](https://doi.org/10.1103/PhysRevLett.93.207205),
  [arXiv:cond-mat/0406440](https://arxiv.org/abs/cond-mat/0406440): the mixed representation,
  an MPS of the vectorized density matrix acted on by superoperators, and `thermalize`, its
  evolution in imaginary time.
- A. H. Werner, D. Jaschke, P. Silvi, M. Kliesch, T. Calarco, J. Eisert and S. Montangero,
  *Positive Tensor Network Approach for Simulating Open Quantum Many-Body Systems*, Phys. Rev.
  Lett. 116, 237201 (2016),
  [doi:10.1103/PhysRevLett.116.237201](https://doi.org/10.1103/PhysRevLett.116.237201),
  [arXiv:1412.5746](https://arxiv.org/abs/1412.5746), and C. D. White, M. Zaletel,
  R. S. K. Mong and G. Refael, *Quantum dynamics of thermalizing systems*, Phys. Rev. B 97,
  035127 (2018), [doi:10.1103/PhysRevB.97.035127](https://doi.org/10.1103/PhysRevB.97.035127),
  [arXiv:1707.01506](https://arxiv.org/abs/1707.01506): what the `cutoff` of a density matrix
  bounds, see [Truncation](@ref).

## Time evolution

- J. Haegeman, J. I. Cirac, T. J. Osborne, I. Pižorn, H. Verschelde and F. Verstraete,
  *Time-Dependent Variational Principle for Quantum Lattices*, Phys. Rev. Lett. 107, 070601
  (2011), [doi:10.1103/PhysRevLett.107.070601](https://doi.org/10.1103/PhysRevLett.107.070601),
  [arXiv:1103.0936](https://arxiv.org/abs/1103.0936), and J. Haegeman, C. Lubich, I. Oseledets,
  B. Vandereycken and F. Verstraete, *Unifying time evolution and optimization with matrix
  product states*, Phys. Rev. B 94, 165116 (2016),
  [doi:10.1103/PhysRevB.94.165116](https://doi.org/10.1103/PhysRevB.94.165116),
  [arXiv:1408.5056](https://arxiv.org/abs/1408.5056): the time dependent variational principle,
  which `Tdvp` runs through ITensorMPS.
- M. Yang and S. R. White, *Time-dependent variational principle with ancillary Krylov
  subspace*, Phys. Rev. B 102, 094315 (2020),
  [doi:10.1103/PhysRevB.102.094315](https://doi.org/10.1103/PhysRevB.102.094315),
  [arXiv:2005.06104](https://arxiv.org/abs/2005.06104): the global Krylov expansion `Tdvp`
  makes before a step from a small bond dimension.
- M. P. Zaletel, R. S. K. Mong, C. Karrasch, J. E. Moore and F. Pollmann, *Time-evolving a
  matrix product state with long-ranged interactions*, Phys. Rev. B 91, 165112 (2015),
  [doi:10.1103/PhysRevB.91.165112](https://doi.org/10.1103/PhysRevB.91.165112),
  [arXiv:1407.1832](https://arxiv.org/abs/1407.1832): the approximations WI and WII of
  `ApproxW`, and its order 2.
- K. Bidzhiev and G. Misguich, *Out-of-equilibrium dynamics in a quantum impurity model:
  Numerics for particle transport and entanglement entropy*, Phys. Rev. B 96, 195117 (2017),
  [doi:10.1103/PhysRevB.96.195117](https://doi.org/10.1103/PhysRevB.96.195117),
  [arXiv:1707.06111](https://arxiv.org/abs/1707.06111), appendix A: the complex steps of
  `ApproxW` of orders 3 and 4.
- C. Van Loan, *Computing integrals involving the matrix exponential*, IEEE Trans. Autom.
  Control 23, 395 (1978), [doi:10.1109/TAC.1978.1101743](https://doi.org/10.1109/TAC.1978.1101743):
  the integrals of exponentials that the steps of `ApproxW` and those of a mixed `GaussianState`
  need, read in the exponential of a block triangular matrix.
- S. Blanes and P. C. Moan, *Fourth- and sixth-order commutator-free Magnus integrators for
  linear and non-linear dynamical systems*, Appl. Numer. Math. 56, 1519 (2006),
  [doi:10.1016/j.apnum.2005.11.004](https://doi.org/10.1016/j.apnum.2005.11.004): the
  commutator-free Magnus integrator of order 4 of the evolvers depending on time (eq. (43)).
- E. M. Stoudenmire and S. R. White, *Minimally entangled typical thermal state algorithms*,
  New J. Phys. 12, 055026 (2010),
  [doi:10.1088/1367-2630/12/5/055026](https://doi.org/10.1088/1367-2630/12/5/055026),
  [arXiv:1002.1305](https://arxiv.org/abs/1002.1305), and C. Camaño, E. N. Epperly and
  J. A. Tropp, *Successive randomized compression: A randomized algorithm for the compressed
  MPO-MPS product*, Quantum 10, 2022 (2026),
  [doi:10.22331/q-2026-03-10-2022](https://doi.org/10.22331/q-2026-03-10-2022),
  [arXiv:2504.06475](https://arxiv.org/abs/2504.06475): the zip-up product of
  `apply_algo = "zipup"`, and how the products of an MPO by a state compare, `"densitymatrix"`
  included.
- S. Paeckel, T. Köhler, A. Swoboda, S. R. Manmana, U. Schollwöck and C. Hubig, *Time-evolution
  methods for matrix-product states*, Ann. Phys. 411, 167998 (2019),
  [doi:10.1016/j.aop.2019.167998](https://doi.org/10.1016/j.aop.2019.167998),
  [arXiv:1901.05824](https://arxiv.org/abs/1901.05824): the truncations of the zip-up product, a
  tenth of the cutoff and twice the bond dimension, which ITensorMPS recommends after them.

## Searches

- J. Cui, J. I. Cirac and M. C. Bañuls, *Variational Matrix Product Operators for the Steady
  State of Dissipative Quantum Systems*, Phys. Rev. Lett. 114, 220601 (2015),
  [doi:10.1103/PhysRevLett.114.220601](https://doi.org/10.1103/PhysRevLett.114.220601),
  [arXiv:1501.06786](https://arxiv.org/abs/1501.06786): `steady_state`, by dmrg on
  ``L^\dagger L``.

## Gaussian states

- S. Bravyi and D. Gosset, *Complexity of quantum impurity problems*, Commun. Math. Phys. 356,
  451 (2017), [doi:10.1007/s00220-017-2976-9](https://doi.org/10.1007/s00220-017-2976-9),
  [arXiv:1609.00735](https://arxiv.org/abs/1609.00735): the conventions of the Majoranas of
  `GaussianState`, the overlap of two Gaussian states with its phase (eq. (27)) and their
  transition matrix (eq. (28)).
- S. Bravyi and R. König, *Classical simulation of dissipative fermionic linear optics*, Quantum
  Inf. Comput. 12, 925 (2012), [arXiv:1112.2184](https://arxiv.org/abs/1112.2184): the
  probability of an occupation in a mixed Gaussian state (eq. (20)), the covariance after
  measuring it (eq. (22)), its evolution under a Lindbladian of quadratic H and linear jumps,
  an equation of Lyapunov (eq. (32)), and its steady state (eq. (33)), used for the mixed states
  of `GaussianState`.
- T. Prosen, *Third quantization: a general method to solve master equations for quadratic open
  Fermi systems*, New J. Phys. 10, 043026 (2008),
  [doi:10.1088/1367-2630/10/4/043026](https://doi.org/10.1088/1367-2630/10/4/043026),
  [arXiv:0801.1257](https://arxiv.org/abs/0801.1257): the jumps of the last site of a spin chain,
  a linear operator behind the parity of all modes, which has the dissipator of the linear
  operator (below eq. (55)), accepted on the mixed states of `GaussianState`.
- D. Karevski and T. Platini, *Quantum Nonequilibrium Steady States Induced by Repeated
  Interactions*, Phys. Rev. Lett. 102, 207207 (2009),
  [doi:10.1103/PhysRevLett.102.207207](https://doi.org/10.1103/PhysRevLett.102.207207),
  [arXiv:0904.3527](https://arxiv.org/abs/0904.3527): the exact steady state of the XX chain driven
  at its ends, a flat profile inside, which the example of `GaussianState` in the manual shows and
  its tests check.
- J. Surace and L. Tagliacozzo, *Fermionic Gaussian states: an introduction to numerical
  approaches*, SciPost Phys. Lect. Notes 54 (2022),
  [doi:10.21468/SciPostPhysLectNotes.54](https://doi.org/10.21468/SciPostPhysLectNotes.54),
  [arXiv:2111.08343](https://arxiv.org/abs/2111.08343): the imaginary-time evolution of a mixed
  Gaussian state by the product rule of covariances (eqs. (86), (153)), which `thermalize` of a mixed
  `GaussianState` replaces by a purification, the product rule needing `I - M` invertible.
- B. Dias and R. König, *Classical simulation of non-Gaussian fermionic circuits*, Quantum 8,
  1350 (2024), [doi:10.22331/q-2024-05-21-1350](https://doi.org/10.22331/q-2024-05-21-1350),
  [arXiv:2307.12912](https://arxiv.org/abs/2307.12912): a Gaussian state held with its phase by
  its covariance and its amplitude on a Fock state, and the measurement of an occupation that
  keeps that phase.
- T. Guaita, L. Hackl and T. Quella, *Representation theory of Gaussian unitary
  transformations for bosonic and fermionic systems*, Ann. Phys. 495, 170757 (2026),
  [doi:10.1016/j.aop.2026.170757](https://doi.org/10.1016/j.aop.2026.170757),
  [arXiv:2409.11628](https://arxiv.org/abs/2409.11628): the random reference state that carries
  the relative phases of a superposition of Gaussian states through Gaussian unitaries.
- G. P. Chen and G. E. Scuseria, *Robust formulation of Wick's theorem for computing matrix
  elements between Hartree–Fock–Bogoliubov wavefunctions*, J. Chem. Phys. 158, 231102 (2023),
  [doi:10.1063/5.0156124](https://doi.org/10.1063/5.0156124),
  [arXiv:2304.13780](https://arxiv.org/abs/2304.13780): the cross terms of components of a
  superposition orthogonal by structure, their zero singular values kept apart without dividing
  by them.
- X. Cao, A. Tilloy and A. De Luca, *Entanglement in a fermion chain under continuous monitoring*,
  SciPost Phys. 7, 024 (2019),
  [doi:10.21468/SciPostPhys.7.2.024](https://doi.org/10.21468/SciPostPhys.7.2.024),
  [arXiv:1804.04638](https://arxiv.org/abs/1804.04638), and Y. Le Gal, X. Turkeshi and M. Schirò,
  *Entanglement dynamics in monitored systems and the role of quantum jumps*, PRX Quantum 5, 030329
  (2024), [doi:10.1103/PRXQuantum.5.030329](https://doi.org/10.1103/PRXQuantum.5.030329),
  [arXiv:2312.13419](https://arxiv.org/abs/2312.13419), appendix F: the trajectories of
  `GaussianState` under jumps that are not unitary, the state carried by its annihilators and
  orthonormalised by a QR at each step, the factor on its norm read from that QR (their eq. (F12),
  here extended to states that do not conserve the number of fermions).
- E. Knill, *Fermionic linear optics and matchgates*,
  [arXiv:quant-ph/0108033](https://arxiv.org/abs/quant-ph/0108033), section 6, and X. Lyu and
  K. Bu, *Displaced fermionic Gaussian states and their classical simulation*,
  [arXiv:2411.18517](https://arxiv.org/abs/2411.18517): the extra mode that makes the odd
  operators of `GaussianState` even, its Majorana multiplying them.
- M. Wimmer, *Algorithm 923: Efficient numerical computation of the Pfaffian for dense and
  banded skew-symmetric matrices*, ACM Trans. Math. Softw. 38, 30 (2012),
  [doi:10.1145/2331130.2331138](https://doi.org/10.1145/2331130.2331138),
  [arXiv:1102.3440](https://arxiv.org/abs/1102.3440): the Parlett-Reid elimination of the
  pfaffians of `GaussianState`.

## States and measurements

- M. T. Fishman and S. R. White, *Compression of correlation matrices and an efficient method
  for forming matrix product states of fermionic Gaussian states*, Phys. Rev. B 92, 075132
  (2015), [doi:10.1103/PhysRevB.92.075132](https://doi.org/10.1103/PhysRevB.92.075132),
  [arXiv:1504.07701](https://arxiv.org/abs/1504.07701): the circuit of rotations that builds
  `slater_state` and `fermi_sea`.
- A. J. Ferris and G. Vidal, *Perfect sampling with unitary tensor networks*, Phys. Rev. B 85,
  165146 (2012), [doi:10.1103/PhysRevB.85.165146](https://doi.org/10.1103/PhysRevB.85.165146),
  [arXiv:1201.3974](https://arxiv.org/abs/1201.3974): the perfect sampling of `sample` on a pure
  state, through ITensorMPS.
- H. Shapourian, K. Shiozaki and S. Ryu, *Partial time-reversal transformation and
  entanglement negativity in fermionic systems*, Phys. Rev. B 95, 165101 (2017),
  [doi:10.1103/PhysRevB.95.165101](https://doi.org/10.1103/PhysRevB.95.165101),
  [arXiv:1611.07536](https://arxiv.org/abs/1611.07536): the fermionic partial transpose, which
  `log_negativity` does not compute, refusing two parts holding fermions.

## Samples and statistics

- G. L. Steele, D. Lea and C. H. Flood, *Fast splittable pseudorandom number generators*, in
  Proceedings of the 2014 ACM International Conference on Object Oriented Programming Systems
  Languages & Applications (OOPSLA), 453 (2014),
  [doi:10.1145/2660193.2660195](https://doi.org/10.1145/2660193.2660195): the finalizer of
  SplitMix64, which mixes the seed of a run with the number of a sample into the seed of that
  sample.
- D. H. D. West, *Updating mean and variance estimates: an improved method*, Commun. ACM 22, 532
  (1979), [doi:10.1145/359146.359153](https://doi.org/10.1145/359146.359153), and T. F. Chan,
  G. H. Golub and R. J. LeVeque, *Algorithms for computing the sample variance: analysis and
  recommendations*, Am. Stat. 37, 242 (1983),
  [doi:10.1080/00031305.1983.10483115](https://doi.org/10.1080/00031305.1983.10483115): the
  weighted means and variances over the samples, accumulated one sample at a time, each value
  shifted by that of the first sample.
- W. G. Cochran, *Sampling Techniques*, 3rd ed., Wiley (1977), chapter 6: the error on a weighted
  mean of the samples, that of a ratio of means.
