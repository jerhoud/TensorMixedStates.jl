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
- S. Blanes and P. C. Moan, *Fourth- and sixth-order commutator-free Magnus integrators for
  linear and non-linear dynamical systems*, Appl. Numer. Math. 56, 1519 (2006),
  [doi:10.1016/j.apnum.2005.11.004](https://doi.org/10.1016/j.apnum.2005.11.004): the
  commutator-free Magnus integrator of order 4 of the evolvers depending on time.
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

## Searches

- J. Cui, J. I. Cirac and M. C. Bañuls, *Variational Matrix Product Operators for the Steady
  State of Dissipative Quantum Systems*, Phys. Rev. Lett. 114, 220601 (2015),
  [doi:10.1103/PhysRevLett.114.220601](https://doi.org/10.1103/PhysRevLett.114.220601),
  [arXiv:1501.06786](https://arxiv.org/abs/1501.06786): `steady_state`, by dmrg on
  ``L^\dagger L``.

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
