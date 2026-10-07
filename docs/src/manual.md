# Manual

```@contents
Pages = ["manual.md"]
Depth = 3
```

## How to read this manual

This manual presents the objects of TMS one after the other, sites, states, operators, then
the algorithms and the measurements, with small examples run directly from Julia, and ends
with how to check the accuracy of the results. The pages
that follow build on it: [High Level Interface](@ref) shows how to write a whole simulation as
a list of phases, which writes its results to files and can be stopped and resumed, the way
most programs are written; [Conserving a quantity](@ref) tells how the conservation of a
quantity is declared, what it forbids and what it saves; [Threads and performance](@ref) is for
when a computation gets large; and [Extending TMS](@ref) is for sites, phases, algorithms and
representations of one's own. The reference pages give every function in detail.

The examples assume TMS has been imported:

```@example manual
using TensorMixedStates
```

## Sites and Systems

The first step in using TMS is the definition of your quantum system. In TMS, a system is composed of a finite
number of sites numbered from 1 (1, 2, ..., N). These sites may be all identical or not.

There are eight different predefined types of site: `Qubit`, `Qudit`, `Fermion`, `Boson`, `Spin`, `Electron`, `Tj` and `Qboson`.

To use each of these sites and the corresponding predefined operators you need first to import the corresponding module.
For example to use qubits, you need to write

```@example manual
using .Qubits
```

Note the "." before the name and the "s" at the end.

To define a site, call the corresponding constructor, for example

```@example manual
s = Qubit()
```

Four site constructors need an argument: `Qudit(dim)` and `Boson(dim)` for the dimension of
the local Hilbert space, `Spin(s)` for the spin, and `Qboson(q, dim)` for the deformation
parameter and the dimension. For example

```@example manual
using .Bosons, .Spins

s = Boson(4)
```

```@example manual
s = Spin(3/2)
```

You can now define a quantum system by declaring the sites it contains:

```@example manual
system1 = System(10, Qubit())
nothing # hide
```

gives you a system with 10 qubits. Systems may have sites of different types, in which case `System` is given an array of sites

```@example manual
using .Fermions

system2 = System([Qubit(), Boson(4), Fermion()])
nothing # hide
```

gives you a three site system.

## States

A state is either a wave function, its pure representation, written `Pure`, or a density
matrix, its mixed representation, written `Mixed`, which open systems need. Both are stored as
a matrix product state (MPS): a chain of tensors, one per site, joined by links whose size,
the bond dimension, measures how much entanglement, or for a density matrix how many
correlations, the state can hold. A density matrix is stored as the vector of its entries, so
each site has dimension ``d^2`` instead of ``d``, its norm is the Hilbert-Schmidt norm
``\sqrt{\mathrm{tr}\,\rho^\dagger\rho}`` rather than its trace, and it costs more than a wave
function of the same system.

To create a state, we call the `State` constructor

```@example manual
state1 = State{Pure}(system1, "Up")
nothing # hide
```

returns a pure up state in a 10 qubit system. Predefined local states are designated by their name. Here `"Up"` is a predefined state of site `Qubit`.

Sites need not be all in the same local state, in which case we give State an array of local states

```@example manual
state2 = State{Mixed}(system2, ["+", "2", "Occ"])
nothing # hide
```

Here we choose a mixed representation.

States may be added or multiplied by a number (they need to be based on the same system). For example

```@example manual
ghz = (State{Pure}(system1, "Up") + State{Pure}(system1, "Dn")) / sqrt(2)
nothing # hide
```

We can transform a pure representation into a mixed representation by

```@example manual
mixedstate = mix(state1)
nothing # hide
```

For mixed states there is a local mixed state `"FullyMixed"` which corresponds to a density matrix proportional to the identity matrix (that is the infinite temperature state).

If you need a local state which is not predefined, it is possible to pass its vector (or matrix for mixed states) directly. For example, we could also define `state1` by

```@example manual
state1 = State{Pure}(system1, [1., 0.])
nothing # hide
```

### [Prepared states](@id manual-prepared-states)

States that are not products of local states, but have an exact MPS of small bond dimension,
are built by functions of their own, which take the local states they are made of:

```@example manual
myneel = superposition(System(4, Qubit()), [1 => ["Up", "Dn", "Up", "Dn"], 1 => ["Dn", "Up", "Dn", "Up"]])
myghz = ghz_state(System(6, Qubit()), "Up", "Dn")          # (|↑↑…↑⟩ + |↓↓…↓⟩)/√2
mydicke = dicke_state(System(6, Qubit()), 2, "Up", "Dn")   # two sites down, in every way
mydimers = dimer_state(System(6, Spin(1/2)), [(1, 2), (3, 4), (5, 6)], "1/2", "-1/2")
myaklt = aklt_state(System(6, Spin(1)))
mysector = fully_mixed(System(6, Fermion()), N => 3)       # 3 fermions at infinite temperature
measure(mydicke, Z)
```

`superposition` superposes product states with the coefficients given, here the two Néel
states, as `mixture` mixes them with weights in a mixed representation, `w_state` is the Dicke
state of one site, `dimer_state` puts a singlet on each pair of sites
given, `aklt_state` is the ground state of the AKLT chain of spins one, and `fully_mixed` the
fully mixed state of a sector, here of 3 fermions. The Slater determinants and the Fermi seas
of free fermions are built by `slater_state` and `fermi_sea`, see [Fermions](@ref), and the
thermal states by `thermal_state`, see [Thermal states](@ref manual-thermal-states). Each is a `State`, which
`CreateState` takes as its `state`, see [Prepared states](@ref) for the details.

## Conserved quantities

A site can be told that a quantity is conserved, for instance the number of fermions:

```@example manual
mysystem = System(10, Fermion(conserve = N))
mystate = State{Pure}(mysystem, [isodd(i) ? "Occ" : "Emp" for i in 1:10])
nothing # hide
```

The tensors are then block sparse, which makes a large computation smaller and faster, and in
exchange the state stays in the sector it was built in: a local state such as `[1, 1] / √2`,
which superposes two numbers of particles, is refused, and so is an operator of no definite charge.
For a mixed state, a quantity may also be conserved strongly, `Fermion(conserve = strong(N))`.
All of this is described in [Conserving a quantity](@ref).

## Limits

TMS uses Matrix Product States to represent quantum states internally. It is important to control the parameters of this approximation, in particular the maximum bond dimension and the cutoff, the total weight of the singular values a
truncation of the tensors may drop. To achieve this, many functions accept a `Limits` object as keyword argument containing those parameters. It is built thus

```@example manual
lim = Limits(cutoff = 1e-10, maxdim = 50)
```

Any of the arguments may be omitted: `maxdim` and `mindim` then impose nothing, and the
default `cutoff`, `eps()`, only discards what rounding leaves, the singular values below about
`1.5e-8` of the norm, while `cutoff = 0` discards nothing. The truncation is a choice to be
checked, see [Checking the accuracy](@ref).

A third parameter, `mindim`, sets the bond dimension the truncation is not allowed to go below, as in `Limits(cutoff = 1e-10, maxdim = 50, mindim = 4)`. `maxdim` takes precedence when the two conflict.

To apply the constraints on a state, one uses

```@example manual
newstate = truncate(ghz; limits = lim)
nothing # hide
```

Many functions accept such an argument. For example, when adding states instead of

```@example manual
state3 = (state1 + ghz) / 2
nothing # hide
```

one can write

```@example manual
state3 = +(state1, ghz; limits = lim) / 2
nothing # hide
```

## Operators

In TMS, there are two kinds of operators: generic operators and indexed operators. For example,

```@example manual
X
```

represents the ``\sigma_x`` Pauli operator for qubits. This is a *generic operator*, it is not applied to a specific site.

```@example manual
X(3)
```

represents the ``\sigma_x`` Pauli operator applied to the system site number 3. This is an *indexed operator*.
 
Note that all predefined operator names start with a capital letter, so it is better to keep your own identifiers lowercase to prevent name collisions with them.

Lowercase is not a safe harbour by itself, though. `using TensorMixedStates` also brings in
seventy lowercase names, among them `state`, `output`, `measure`, `trace`, `dim`,
`matrix`, `tensor`, `apply`, `norm` and `sample`. Assigning to one of them at the top level
of your program shadows the function for the rest of the file, and if you happen to have
used it before assigning to it, Julia 1.10 and 1.11 refuse the assignment outright with
`cannot assign a value to imported variable`. There is no such issue inside a function,
where `state = ...` is an ordinary local variable. This manual prefixes its own variables with
`my`, as in `mystate` and `myop`, which is one way of staying clear.

The operator system is very rich and flexible. For example, if you want to use this Hamiltonian

```math
H = \sum_{i=1}^{n-1} \sigma_x(i) \sigma_x(i+1)
```

you will simply write

```@example manual
n = 10
h = sum(X(i)X(i+1) for i in 1:n-1)
```

Many operations are defined on generic operators:

- addition, multiplication and power by a number
- tensor product: `X⊗X` is a two site operator (`⊗` is usually obtained by typing \otimes in your editor, just in case, one can also write `tensor(X, X)`)
- `dag` represents the adjoint operator, for example `C` is the `c` operator for fermions and `dag(C)` is ``c^\dagger``.
- `Dissipator` represents a Lindblad dissipator, for example `Dissipator(Sp)` is the dissipator whose jump operator, `Sp`, the ``S^+`` operator, flips a qubit toward up
- `Gate` represents an operator to be applied as a gate on a mixed state. It is useful to define noisy gate operators, for example `0.9Gate(Id) + 0.1Gate(X)` is a noisy gate operator that will apply a ``\sigma_x`` gate 10 percent of the time.
- `Proj` represents an operator that projects on the given state, for example `Proj("Up")`
  projects qubits on the up state, or on an eigenspace, `Proj(Sz => 0)` projecting a spin 1 on
  the state of zero ``S_z``.
- `Basis` numbers the basis states of any site, ``\mathrm{diag}(0, 1, \dots)``: its eigenbasis
  is the basis of the site.
- the functions `exp` and `sqrt`: for example `sqrt(Swap)`
- `controlled` for qubits makes controlled gates: `CX = controlled(X)`

For example one can define the Rxy 2-site operator by

```@example manual
Rxy(t) = exp(-im * t * (X⊗X + Y⊗Y) / 4)
```

If this is not enough to define your favorite operator you can create new ones with `named`,
from an expression, a matrix or a function of the sites, see
[Operators and measurements of one's own](@ref):

```@example manual
myop = named([1 1 ; 1 -1] / √2, "MyOp")
```

Finally from generic operators, we define indexed operators by simply applying them to the corresponding sites

```@example manual
Rxy(0.2)(2, 5)
```

```@example manual
myop(3)
```

An operator of several sites that is a function of an operator, as the exponential in `Rxy`,
can only be applied as a gate, unless it is created with the sites it acts on, see
[Operators and measurements of one's own](@ref).

### Fermions

A fermionic operator is written as it is: `C(i)` destroys a fermion on site `i` and
`dag(C)(j)` creates one on site `j`, the Jordan-Wigner strings being inserted for you wherever
the operator is used, in expectation values, in MPOs and in gates. A hopping term is thus
simply

```@example manual
hopping = sum(dag(C)(i)C(i+1) + dag(C)(i+1)C(i) for i in 1:n-1)
nothing # hide
```

The strings follow the order of the sites: `C(j)` is the matrix of `C` on site `j` with `F`
on every site before it. A product of two fermionic operators placed in increasing order of
the sites, `A(i) * B(j)` with `i < j`, is thus `A * F` on site `i`, `F` on every site in
between and `B` on site `j`, and in the other order it takes the sign of the swap. A matrix
given for an operator of several sites, see [Operators and measurements of one's own](@ref),
is read in the basis of the sites as it is, with no string.

A function of an even fermionic operator of several sites, as the exponential of a hopping
term, `exp(-0.1im * (dag(C) ⊗ C + dag(dag(C) ⊗ C)))(2, 5)`, is applied as a gate with its
strings, on sites in any order and apart. A function of an odd operator mixes the two
parities and is refused.

The ground state of free fermions, the Fermi sea of a quadratic hamiltonian, is built at once
by `fermi_sea`, here with three fermions, and any Slater determinant by `slater_state`, see
[Prepared states](@ref):

```@example manual
sea = fermi_sea(System(n, Fermion()), -hopping, 3)
measure(sea, N)
```

## Algorithms

We can now work with states and operators.

### Gates

Gates are applied with `apply`. In a product of gates the rightmost acts first: here the
Hadamard gate on the first qubit, then the controlled not, which turns two up qubits into a
Bell pair

```@example manual
mybell = apply(controlled(X)(1, 2) * H(1), State{Pure}(System(2, Qubit()), "Up"))
measure(mybell, [Z(1)Z(2), X(1)X(2)])
```

Applying all the gates in a single call is much more efficient than one by one, and the
keyword argument `limits` constrains the truncations made along the way.

### Ground states

The ground state of a hamiltonian is computed by `dmrg`, from a starting state, here a random
state of bond dimension 8. It returns the energy and the ground state, which can then be
measured. For the Ising chain in a transverse field

```@example manual
mysystem = System(10, Qubit())
hamiltonian = -sum(Z(i)Z(i + 1) for i in 1:9) - sum(X(i) for i in 1:10)
energy, ground = dmrg(hamiltonian, RandomState{Pure}(mysystem, 8); nsweeps = 10,
                      limits = Limits(maxdim = 20))
energy
```

```@example manual
measure(ground, X)
```

`nsweeps` is the number of sweeps, passes of dmrg along the chain and back, and `limits` may
give one value per sweep, as in `Limits(maxdim = [10, 20, 50])`.

### Time evolution

Time evolution is done by `tdvp` or `approx_W`, which take an evolver, the time to evolve for
and the state. The evolver follows the convention

```julia
evolver = -im * hamiltonian + dissipators
```

dissipators being accepted on a mixed state only. On a pure state, it integrates
``\frac{d}{dt}|\psi\rangle = -iH|\psi\rangle``, and on a mixed state the Lindblad equation

```math
\frac{d\rho}{dt} = -i[H, \rho]
    + \sum_k \left(L_k \rho L_k^\dagger - \tfrac{1}{2}\{L_k^\dagger L_k, \rho\}\right)
```

each term `Dissipator(L)` giving one jump operator ``L_k``, with ``\hbar = 1``. A rate
``\gamma`` is written `γ * Dissipator(L)`, or equally `Dissipator(sqrt(γ) * L)`, since
`Dissipator(c * L)` is `abs2(c) * Dissipator(L)`. Under ``H = \sum_i \sigma_z^i``, qubits
starting in `"+"` precess, with ``\langle \sigma_x \rangle = \cos 2t``, which is
``\cos 1 \approx 0.5403`` at ``t = 0.5``:

```@example manual
myevolved = tdvp(-im * sum(Z(i) for i in 1:4), 0.5, State{Pure}(System(4, Qubit()), "+");
                 nsteps = 10)
measure(myevolved, X)
```

`nsteps` is the number of steps, each of length `t / nsteps`, and `limits` constrains the
state as before. With no hamiltonian and `Dissipator(Sm)` on each site, the qubits of a mixed
state decay from up to down, with ``\langle \sigma_z \rangle = 2e^{-t} - 1``, which is
``-0.2642`` at ``t = 1``:

```@example manual
myrho = State{Mixed}(System(4, Qubit()), "Up")
mydecayed = tdvp(sum(Dissipator(Sm)(i) for i in 1:4), 1.0, myrho; nsteps = 10)
measure(mydecayed, Z)
```

`approx_W` approximates the exponential of the evolver over a step by MPOs, built from the
approximations WI or WII of Zaletel et al. It is called the same way, with in addition the
order of its approximation, from 1 to 4, and `w`, 1 or 2 for WI or WII, for which
`order = 4, w = 2` is usually a good choice:

```@example manual
mydecayed = approx_W(sum(Dissipator(Sm)(i) for i in 1:4), 1.0, myrho; order = 4, w = 2,
                     nsteps = 10)
measure(mydecayed, Z)
```

An evolver may also depend on time, see [Time dependent evolvers](@ref).

### [Noise](@id manual-noise)

Noise acts on a mixed state, either as a channel applied between gates, or as a term of an
evolver. Three kinds are ready for sites of any type, each as a channel `…_gate(p, …)`, taking
place with probability `p`, and as a generator `…_dissipator(γ, …)`, of rate `γ`:

- [`relaxing_gate`](@ref)`(p, state)` resets a site to `state`, see `SetState`; with `n`,
  `relaxing_gate(p, n, state)` resets `n` sites together, and `relaxing_gate(p, [s1, s2])`
  each site to its state;
- [`depolarizing_gate`](@ref)`(p, n)` is the relaxation towards `"FullyMixed"`;
- [`dephasing_gate`](@ref)`(p, A)` erases the coherences between the eigenspaces of `A`, or
  in the basis of the site without `A`, see `Dephase`.

Evolving under a generator for a time `t` is the channel of probability ``1 - e^{-\gamma t}``.
Dephasing qubits in `"+"` at rate 1 gives ``\langle \sigma_x \rangle = e^{-t}``, which is
``0.6065`` at ``t = 0.5``, both ways:

```@example manual
myplus = State{Mixed}(System(2, Qubit()), "+")
mydephased = tdvp(sum(dephasing_dissipator(1.0)(i) for i in 1:2), 0.5, myplus; nsteps = 10)
measure(mydephased, X)
```

```@example manual
measure(apply(prod(dephasing_gate(1 - exp(-0.5))(i) for i in 1:2), myplus), X)
```

Sites relaxed together are not sites relaxed each on its own: `depolarizing_gate(p, 2)`
depolarizes both sites with probability `p`, where
`depolarizing_gate(p) ⊗ depolarizing_gate(p)` depolarizes each with probability `p`. And
`dephasing_dissipator(γ, A)` damps every coherence between two eigenspaces of `A` at the rate
`γ`, where `Dissipator(A)` damps it at ``(\lambda_a - \lambda_b)^2/2``, faster for distant
eigenvalues; they agree on a qubit, `Dissipator(Z)` being `dephasing_dissipator(2, Z)`.

The relaxation, and so the depolarization, changes the charges and is refused on sites
conserving something strongly, where a dephasing in the eigenbasis of an operator commuting
with the charge, as `N`, is accepted. Any other channel is a sum of gates, `Gate(K)` being
``\rho \mapsto K \rho K^\dagger``.

On a qubit, `"Up"` is ``|0\rangle`` and `"Dn"` ``|1\rangle``, the excitation `N` counts: the
decay of an excitation, from `"Dn"` to `"Up"`, is `Dissipator(Sp)`, and `Dissipator(Sm)`, which
takes `"Up"` to `"Dn"` in the evolution above, excites. The module `Qubits` gives the decay as
[`amplitude_damping_gate`](@ref) and [`amplitude_damping_dissipator`](@ref), and the relaxation
of times ``T_1`` and ``T_2`` of a qubit, its decay and its dephasing, as
[`thermal_relaxation_gate`](@ref) and [`thermal_relaxation_dissipator`](@ref).

### Steady states

`steady_state` looks for the state ``\rho`` with ``L\rho = 0``, ``L`` being the evolver,
starting from a mixed state, by running dmrg on ``L^\dagger L``. It returns two things: the
residual ``\|L\rho\|^2``, for ``\rho`` of unit Hilbert-Schmidt norm, which is zero at a true
steady state and tells how well the search converged, and the steady state itself, normalized
to trace one. A Lindbladian may have several steady states, see
[Checking the accuracy](@ref). Qubits decaying toward down at rate 1 and pumped toward up at rate
0.5 settle at ``\langle \sigma_z \rangle = -1/3``:

```@example manual
value, mysteady = steady_state(sum(Dissipator(Sm)(i) + 0.5Dissipator(Sp)(i) for i in 1:4),
                               myrho; nsweeps = 10)
measure(mysteady, Z)
```

### [Thermal states](@id manual-thermal-states)

`thermal_state` takes a mixed state ``\rho`` to ``e^{-\beta H/2} \rho \, e^{-\beta H/2}``,
normalized, by tdvp in imaginary time. From `"FullyMixed"`, the state at infinite temperature,
it gives the thermal state ``e^{-\beta H}/Z``. It returns the logarithm of the trace the state
would have without being normalized, here ``\log Z - 4 \log 2``, and the state. Qubits in the
field ``H = -\sum_i \sigma_z^i`` have ``\langle \sigma_z \rangle = \tanh \beta``, which is
``0.4621`` at ``\beta = 0.5``, and ``\log Z - 4 \log 2 = 4 \log \cosh \beta \approx 0.4805``:

```@example manual
mylog, mythermal = thermal_state(-sum(Z(i) for i in 1:4), 0.5,
                                 State{Mixed}(System(4, Qubit()), "FullyMixed"); nsteps = 10)
mylog
```

```@example manual
measure(mythermal, Z)
```

For more details, see the [Algorithms](algorithms.md) page of the reference or the inline
help.

## [Measurements](@id manual-measurements)

Once we have created a state, we may want to measure it. Take for example a three qubit state

```@example manual
mystate = State{Pure}(System(3, Qubit()), ["+", "Up", "+"])
nothing # hide
```

```@example manual
result = measure(mystate, X(1)X(3))
```

will give ``\langle \psi | \sigma_x^1 \sigma_x^3 | \psi \rangle``

```@example manual
result = measure(mystate, X)
```

will give the array of the ``\langle \psi | \sigma_x^i | \psi \rangle``

```@example manual
result = measure(mystate, (X, Y))
```

will give the matrix of the ``\langle \psi | \sigma_x^i \sigma_y^j | \psi \rangle``, complex
since its diagonal is ``\langle \sigma_x \sigma_y \rangle = i \langle \sigma_z \rangle``

Fermionic operators follow two rules. A single one given for every site, as `C`, is refused by
`measure` and `expect1`, since its value vanishes on every state of definite parity; placed on
a site, as `C(3)`, it is accepted, for a state that superposes parities. In a correlation, as
`(dag(C), C)`, the two operators must be both fermionic or both not.

We can also measure properties of the state as a whole, with state functions such as
`Trace`, `Purity`, `EntanglementEntropy(cut)` or `Fidelity(ref)`: the
[Measurements](measurements.md) page has the table of them all.

A few sites can be looked at together: `ReducedDensityMatrix(positions)` measures their
density matrix, `VonNeumannEntropy(positions)` their entropy and `LogNegativity(a, b)` the
entanglement between two parts of them, which on a mixed state the entropy no longer tells.
Both are ``\log 2`` for one qubit of a Bell pair:

```@example manual
mybell = apply(controlled(X)(1, 2) * H(1), State{Pure}(System(3, Qubit()), "Up"))
measure(mybell, [VonNeumannEntropy([1]), LogNegativity([1], [2])])
```

We can also ask for several measurements at the same time

```@example manual
results = measure(mystate, [X, X(2)Z(3), (X, Y), Trace, MemoryUsage])
```

A measurement can also be drawn, as in an experiment: `probabilities(state, op)` gives the
results of measuring an operator placed on sites with their probabilities, `sample(state, op)`
draws one, and `collapse(state, op)` draws one and gives the state projected onto it. Any
Hermitian operator of a few sites can be measured, and on any number of sites a string of
Pauli operators or the number of particles of a region. The parity of two qubits in `"+"` is
±1 with probability 1/2, and either result leaves an entangled pair:

```@example manual
mytwo = State{Pure}(System(2, Qubit()), "+")
probabilities(mytwo, Z(1) * Z(2))
```

```@example manual
myparity, mypair = collapse(mytwo, Z(1) * Z(2))
myparity, measure(mypair, X(1)X(2))
```

For more details see the reference or the inline help.

## Checking the accuracy

Three main approximations set the error of a result, and the user chooses all three: the
truncation of the states, given by `limits`, the time step of an evolution or of
`thermal_state`, and the convergence of a search, the sweeps of `dmrg` and `steady_state`. A
result can be trusted once it no longer moves when each of them is improved: larger limits, a
smaller time step, more sweeps. The measurements taken along a run tell which of them limits
the accuracy. Two more have defaults: the tolerance of the Krylov method of each local step,
see `Krylov`, and the truncation `mpo_limits` of the ``L^\dagger L`` of `steady_state`, none by
default; truncated, it makes the value `steady_state` returns the residual of the truncated
operator.

### Truncation

`MaxLinkdim` tells when the bond dimension reaches `maxdim`, from which point the results
depend on it, and `EntanglementEntropy(cut)`, the operator space entanglement entropy on a
mixed state, shows the entanglement that makes the bond dimension grow. A truncation by
`cutoff` shows in neither, and only running again tells it. `Limits()` discards the singular
values below about 1.5e-8 of the norm, which keeps rounding from making the bond dimension
grow and may leave an error of that order: the truncation that makes a computation affordable
is a choice, checked by running again with a larger `maxdim` and a smaller `cutoff`.

### Time step

An evolution is checked by halving its time step, and by comparing `Tdvp` with `ApproxW`,
whose errors have different origins. `ApproxW` of order `k` makes an error of order
``\tau^k`` in the time step ``\tau`` when the evolver does not depend on time. When it does,
`tdvp` and `ApproxW` of order 1 or 2 take its functions at the middle of each step, which caps
them at order 2, while `ApproxW` of order 3 or 4 keeps the order of its approximation, see
[Time dependent evolvers](@ref). `tdvp` projects the evolution on the states of the bond
dimension the state has, its steps on two sites letting that dimension grow: from a state of
small bond dimension, as a product state, or for `thermal_state`, compare with
`expand_period = 1`, which enlarges the bond dimension by a global Krylov expansion before
each step. A thermal state has checks of its own, see [Thermal states](@ref).

### Indicators of a density matrix

The exact density matrix is hermitian under a Lindbladian, under gates and noisy gates, which
act on both of its sides, and for a thermal or a steady state, and its trace stays one when the
evolution preserves the trace, as a Lindbladian does. What deviates from that is error, which
three state functions measure along a run:

- `TraceError`, ``1 - \mathrm{tr}\,\rho``, the drift of the trace. Expectation values are
  divided by the trace, which hides the drift but not the error it comes from.
- `HermiticityError`, the squared norm of the anti-hermitian part of ``\rho``, relative to
  that of ``\rho``. That part being error, the square root of `HermiticityError` is a lower
  bound of the error relative to the norm of ``\rho``, in the Hilbert-Schmidt norm: a value of
  `1e-6` proves an error of `1e-3` at least. It bounds the error from below only, an error
  without an anti-hermitian part going unseen. `hermitianize_period` removes the
  anti-hermitian part, and with it what `HermiticityError` shows: measure it before deciding
  to hermitianize.
- `Purity`, ``\mathrm{tr}(\rho^\dagger \rho) / (\mathrm{tr}\,\rho)^2``, lies between ``1/D``
  and 1 for a density matrix of dimension ``D``. Above 1, the state is no longer a density
  matrix: it has lost its positivity, which nothing in a matrix product state guarantees,
  unless `HermiticityError` shows that it has lost its hermiticity.

### Conserved quantities

A quantity declared conserved on the sites, see [Conserving a quantity](@ref), keeps the
errors out of the sectors the exact state never reaches: the tensors have no blocks there, so
that neither truncation nor rounding can put anything in them. A pure state keeps its charge
exactly, a mixed state under a weak symmetry keeps no coherence between two sectors, and
under a strong one stays in its sector. A quantity the exact evolution conserves without its
being declared, as the energy under a Hamiltonian, which no site can declare, is still a
check: its drift along a run is error, though its constancy proves nothing.

### Searches

`dmrg` runs all its sweeps and, in a `GroundState` or a `SteadyState` phase, stops earlier
once the energy changes by less than `tol` between two sweeps, which a search stuck in a
metastable state does as well. `Variance(hamiltonian)`, zero exactly for an eigenstate, checks
the state found, see `variance`, but does not tell the ground state from another eigenstate,
which searches from other states, or in other sectors, do. Extrapolating the energy to
zero variance as `maxdim` grows estimates the exact energy, and the error of the last one.

The value `steady_state` returns, the residual ``\|L\rho\|^2`` for ``\rho`` of unit
Hilbert-Schmidt norm, is zero for a steady state. The state it finds is hermitian when the
steady state is unique and the search has converged, and it warns when its `HermiticityError`
exceeds `1e-6`: a state that is not hermitian comes from a search that has not converged, or
from a Lindbladian with several steady states, of which dmrg returns any combination: one in
each sector of a quantity that its Hamiltonian and its jump operators all commute with, for
instance. The one a system reaches depends on the state it starts from: it is the limit of
its evolution in time, when that converges, which `tdvp` gives. When the quantity is one a
site can conserve, declaring it `strong` keeps `steady_state` in the sector of the starting
state, see [Conserving a quantity](@ref). A hermitian state does not prove the steady state
unique.
