# The gates of quantum circuits, gathered in the module Circuits: those of the qubits and the
# qudits, the rotations, of half the angle in the exponent as in OpenQASM and Qiskit, and the
# controlled gates.

export Circuits

"""
    H

the Hadamard gate of a qubit, ``(X + Z)/\\sqrt{2}``, and on a qudit the Fourier operator,
``|j\\rangle \\mapsto \\sum_k \\omega^{jk} |k\\rangle / \\sqrt{d}``, ``\\omega = e^{2i\\pi/d}``
"""
H

"""
    S

the phase gate of a qubit, ``\\mathrm{diag}(1, i)``, and on a qudit the phase operator, of order
``d`` for odd ``d`` and ``2d`` for even ``d``, see `Qudit`
"""
S

"""
    T

the gate ``\\mathrm{diag}(1, e^{i\\pi/4})`` of a qubit, the square root of `S`
"""
T

"""
    Rx(θ)

the rotation by the angle `θ` about the axis `x`, ``e^{-i\\theta S_x}``: on a qubit
``e^{-i\\theta X/2}``, as in OpenQASM, Qiskit or Nielsen and Chuang, and on a spin the rotation
of the spin. `Ry` and `Rz` turn about `y` and `z`.

# Examples

    Rx(π / 2)(1)
    controlled(Ry(0.3))(1, 2)
"""
Rx(θ) = Operator{1}("Rx($θ)", exp(-im * θ * Sx))

"""
    Ry(θ)

the rotation by the angle `θ` about the axis `y`, ``e^{-i\\theta S_y}``, see `Rx`
"""
Ry(θ) = Operator{1}("Ry($θ)", exp(-im * θ * Sy))

"""
    Rz(θ)

the rotation by the angle `θ` about the axis `z`, ``e^{-i\\theta S_z}``, see `Rx`. On a qubit it
is `Phase(θ)` times the global phase ``e^{-i\\theta/2}``, which a controlled gate makes relative.
"""
Rz(θ) = Operator{1}("Rz($θ)", exp(-im * θ * Sz))

"""
    Rxx(θ)

the Ising coupling of two sites by the angle `θ`, ``e^{-2i\\theta S_x \\otimes S_x}``: on two
qubits ``e^{-i\\theta X \\otimes X/2}``, as in Qiskit or PennyLane. `Ryy` and `Rzz` couple `Sy`
and `Sz`.

# Examples

    Rzz(π / 2)(1, 2)        # the square root of Z ⊗ Z, up to a global phase
"""
Rxx(θ) = Operator{2}("Rxx($θ)", exp(-2im * θ * (Sx ⊗ Sx)))

"""
    Ryy(θ)

the Ising coupling of two sites by the angle `θ`, ``e^{-2i\\theta S_y \\otimes S_y}``, see `Rxx`
"""
Ryy(θ) = Operator{2}("Ryy($θ)", exp(-2im * θ * (Sy ⊗ Sy)))

"""
    Rzz(θ)

the Ising coupling of two sites by the angle `θ`, ``e^{-2i\\theta S_z \\otimes S_z}``, see `Rxx`
"""
Rzz(θ) = Operator{2}("Rzz($θ)", exp(-2im * θ * (Sz ⊗ Sz)))

"""
    CX

the controlled not, `controlled(X)`, also named `CNOT`: on two qubits, `X` on the second when
the first is in `"Dn"`, `"1"`.

# Examples

    apply(CX(1, 2) * H(1), State{Pure}(System(2, Qubit()), "Up"))      # a Bell pair
"""
const CX = controlled(X)

"""
    CNOT

the controlled not, the same gate as `CX`
"""
const CNOT = CX

"""
    CY

the controlled `Y`, `controlled(Y)`
"""
const CY = controlled(Y)

"""
    CZ

the controlled `Z`, `controlled(Z)`
"""
const CZ = controlled(Z)

"""
    CCX

the Toffoli gate, `controlled(CX)`, also named `Toffoli`: on three qubits, `X` on the third when
the first two are in `"Dn"`.
"""
const CCX = controlled(CX)

"""
    Toffoli

the Toffoli gate, the same gate as `CCX`
"""
const Toffoli = CCX

"""
    using .Circuits

module of the gates of quantum circuits: `H`, `S` and `T`, which `H` and `S` extend to qudits,
`Phase`, the gate `Sum` of two qudits, the rotations `Rx`, `Ry`, `Rz`, `Rxx`, `Ryy` and `Rzz`,
and the controlled gates `CX` or `CNOT`, `CY`, `CZ`, and `CCX` or `Toffoli`.
"""
module Circuits
import ..H, ..S, ..T, ..Phase, ..Sum, ..Rx, ..Ry, ..Rz, ..Rxx, ..Ryy, ..Rzz, ..CX, ..CNOT, ..CY,
       ..CZ, ..CCX, ..Toffoli
export H, S, T, Phase, Sum, Rx, Ry, Rz, Rxx, Ryy, Rzz, CX, CNOT, CY, CZ, CCX, Toffoli
end
