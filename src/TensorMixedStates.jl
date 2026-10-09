# The module TensorMixedStates: what it imports from its dependencies, and its source files in
# the order they are read. Version 2 is written from its lowest layer up, the code of 1.x
# staying on the branch dev.

"""
    TensorMixedStates

a library for simulating closed and open quantum systems with matrix product states. A state
is either pure, its wave function as an MPS, or mixed, its density matrix as an MPS, on which
Lindbladian evolution and noisy gates act.
"""
module TensorMixedStates

using ITensors, ITensorMPS, LinearAlgebra, Random

end
