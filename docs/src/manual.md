# Manual

## Import

To use TMS, you must first import it with

    using TensorMixedStates

## Sites and Systems

The first step in using TMS is the definition of your quantum system. In TMS, a system is composed of a finite
number of sites numbered from 1 (1, 2, ..., N). These sites may be all identical or not.

There are seven different predefined types of site: `Qubit`, `Fermion`, `Boson`, `Spin`, `Electron`, `Tj` and `Qboson`.

To use each of these sites and the corresponding predefined operators you need first to import the corresponding module.
For example to use qubits, you need to write

    using .Qubits

Note the "." before the name and the "s" at the end.

To define a site just call the corresponding creator for example

    s = Qubit()

Some site creators need arguments: `Boson` (for the maximum occupancy) and `Spin`, for example

    s = Boson(4)
    s = Spin(3/2)

You can now define a quantum system by declaring the sites it contains:

    system1 = System(10, Qubit())

gives you a system with 10 qubits. Systems may have different types of site, in this case you must feed `System` with an array of sites

    system2 = System([Qubit(), Boson(4), Fermion()])

gives you a three site system.

## States

States may be in pure or mixed representation, these two possibilities are represented in TMS, by `Pure()` or `Mixed()` (note the parenthesis).

To create a state, we call the State creator

    state1 = State(Pure(), system1, "Up")

returns a pure up state in a 10 qubit system. Predefined local states are designated by there name. Here `"Up"` is a predefined state of site `Qubit`.

All sites need not be all in the same local states, in which case we give State an array of local states

    state2 = State(Mixed(), system2, ["+", "2", "Occ"])

Here we chose a mixed representation.

States may be added or multiplied by a number (they need to be based on the same system). For example

    ghz = (State(Pure(), system1, "Up") + State(Pure(), system1, "Dn")) / 2

We can transform a pure representation into a mixed representation by

    mixedstate = mix(purestate)

For mixed states there is a local mixed state `"FullyMixed"` which correspond to a density matrix proportional to the identity matrix.

If you need a local state which is not predefined, it is possible to pass its vector (or matrix for mixed states) directly, For example, we could also define `state1` by

    state1 = State(Pure(), system1, [1., 0.])

## Limits

TMS uses Matrix Product State to internally represent quantum states. It is important to control the parameters of this approximation, in particular the maximum bond dimension and the cutoff on singular values. To achieve this, many functions accept a `Limits` object as keyword argument containing those parameters. It is build thus

    lim = Limits(cutoff = 1e-10, maxdim = 50)

each (or both) of the arguments may be omitted in which case it corresponds to an absence of constraint for this parameter. In particular, `Limits()` represents no constraint.

To apply the constraints on a state, one uses

    newstate = truncate(oldstate; limits = lim)

Many functions accept such an argument. For example, when adding states instead of

    state = (state1 + state2) / 2

One can write

    state = +(state1, state2; limits = lim) / 2

## Operators

In TMS, there are two kinds of operators: generic operators and indexed operators. For example,

    X

represents the ``\sigma_x`` Pauli operator for qubits. This is a *generic operator*, it is not applied to a specific site.

    X(3)

represents the ``\sigma_x`` Pauli operator applied to the system site number 3. This is an *indexed operator*.
 
Note that all predefined operator names start with a capital letter, so it is better to keep your own identifiers lowercase to prevent name collisions.

The operator system is very rich and flexible. For example, if you want to use this Hamiltonian

```math
H = \sum_{i=1}^{n-1} \sigma_x(i) \sigma_x(i+1)
```

you will simply write

    h = sum(X(i)X(i+1) for i in 1:n-1)

Many operations are defined on generic operators:

- addition, multiplication and power by a number
- tensor product: `X⊗X` is a two site operator (one can also write `tensor(X, X)`)
- `dag` represents the adjoint operator, for example `C` is the `c` operator for fermions and `dag(C)` is ``c^\dagger``.
- `Dissipator` represents a Lindblad dissipator, for example `Dissipator(Sp)` is the jump operator that may flip a qubit toward up (`Sp` is the ``S^+`` operator)
- `Gate` represents an operator to be applied as a gate on a mixed state. It is useful to define noisy gate operators, for example `0.9 Gate(Id) + 0.1 Gate(X)` is a noisy gate operator that will apply an ``\sigma_x`` gate 10 percent of the time.
- `Proj` represents an operator that projects on the given state, for example `Proj("Up")` projects qubits on the up state.
- the functions `exp` and `sqrt`: for example `sqrt(Swap)`
- `controlled` for qubits makes controlled gates: `CX = controlled(X)`

For example one can define the Rxy 2-site operator by

    Rxy(t) = exp(-im * t * (X⊗X + Y⊗Y) / 4)

If this is not enough to define your favorite operator you can create new ones by specifying their matrix

    myop = Operator("MyOp", [1 1 ; 1 -1] / √2)

For multiple site or mixed operators you must specify the type

    Swap = Operator{Pure, 2}("Swap", [1 0 0 0 ; 0 0 1 0 ; 0 1 0 0 ; 0 0 0 1])

Finally from generic operators, we define indexed operators by simply applying them to the corresponding sites

    Rxy(0.2)(2, 5)
    myop(3)
    Swap(4, 7)

In the case of Hamiltonian or Lindbladian evolution the Hamiltonian part is to be multiplied by -im:

    evolver = -im * hamiltonian + dissipators

## Algorithms

We can now work with states and operators.

We can apply gates with `apply`

    newstate = apply(gates, oldstate; limits)

the `gates` argument is an indexed operator representing the gates to apply

the keyword argument limits fixes the constraints to apply

We can compute ground states with `dmrg`

    groundstate, energy = dmrg(hamiltonian, startstate; options...)

the options are `limits` to set constraints and `nsweeps` to fix the number of sweeps among others.

We can do time evolution with `tdvp` and `approx_W`

    newstate = tdvp(evolver, time, oldstate; options...)
    newstate = approx_W(evolver, time, oldstate; options...)

the options are `limits` for the constraints, `nsweeps` for the number of step to do and for `approx_W`, `order` and `w` for the parameters of the algorithm (`order = 4, w = 2` are usually good)

For more details, see the reference or the inline help.

## Measurements

Once we have created a state, we may want to measure it.

    result = measure(state, X(1)X(3))

will give ``\langle \psi | \sigma_x^1 \sigma_x^3 | \psi \rangle``

    result = measure(state, X)

will give the array of the ``\langle \psi | \sigma_x^i | \psi \rangle``

    result = measure(state, (X, Y))

will give the matrix of the ``\langle \psi | \sigma_x^i \sigma_y^j | \psi \rangle``

We can also measure other properties like `Trace`, `TraceError`, `Trace2`, `Purity`, `Linkdim` and `MemoryUsage`.

We can also ask for several measurements at the same time

    results = measure(state, [X, X(2)Z(3), (X, Y), Trace, MemoryUsage])

For more details see the reference or the inline help.

## High Level Interface

Most simulations follow the same pattern: start from some simple state, make some evolution and make measurements during or after the evolution and save the results to file. For these simple cases, TMS proposes a simpler interface.

A simple simulation follows a single state through a certain number of phases which act in a simple way on the state and make measurements during and/or after the evolution and save the results to file.

The following phases are available:

- `CreateState` : create a simple state
- `DMRG` : compute the ground state (requires a pure state)
- `ToMixed` : go from pure representation to mixed representation
- `Evolve` : do Hamiltonian or Lindbladian evolution
- `Gates` : apply some gates
- `Partial_trace` : trace the system over some sites (requires a mixed state)

with these phases we define a `SimData` object that describes the simulation and finally, we call

    runTMS(simdata)

which executes the simulation.

As an example, here is the complete code for such a simple simulation:

```
using TensorMixedStates, .Fermions

hamiltonian(n) = -sum(dag(C)(i)C(i+1)+dag(C)(i+1)C(i) for i in 1:n-1)
dissipators(n, gamma) = sum(Dissipator(sqrt(4gamma) * N)(i) for i in 1:n)

sim_data(n, gamma, step) = SimData(
    name = "Fermion tight-binding chain with dephasing noise",
    phases = [
        CreateState(
            type = Mixed(),
            sytem = System(n, Fermion()),
            state = [ iseven(i) ? "Occ" : "Emp" for i in 1:n ]),
        Evolve(
            duration = 4,
            time_step = step,
            algo = Tdvp(),
            evolver = -im*hamiltonian(n) + dissipators(n, gamma),
            limits = Limits(cutoff = 1e-30, maxdim = 100),
            measures = [
                "density.dat" => N,
                "OSEE.dat" => EE(div(n, 2))
            ]
        )
    ]
)

runTMS(sim_data(40, 1., 0.05))
```

For more information, see the reference or inline help for each phase, SimData and runTMS.

