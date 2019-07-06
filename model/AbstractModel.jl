export AbstractModel, IsingModel
export heisenberg_ij, hamiltonian, heisenberg_term, ground_state, energy, energy_exact, get_bonds, energy, heisenberg_2d, nspins

abstract type AbstractModel{D} end

nspins(model::AbstractModel) = prod(size(model)) #get spins count

"""
    energy_exact(tc, model::AbstractModel) -> Float64

Exact ground state energy.
"""
function energy_exact(tc, model::AbstractModel)
    nqubit = nspins(tc)
    expect(hamiltonian(model), state_exact(tc)) |> real #get expectation value of hamiltonian
end

"""
    ground_state(model::AbstractModel) -> DefaultRegister

Get the exact ground state of a model.
"""
function ground_state(model::AbstractModel)
    # get the ground state
    hami = hamiltonian(model)
    E, v = eigsolve(mat(hami), 1, :SR)
    E[1], ArrayReg(v[1])
end

"""
    get_bonds(model::AbstractHeisenberg) -> Vector

Get the weighted bonds of a Heisenberg model in the form Tuple(i,j,w_ij).
"""
function get_bonds end

include("IsingModel.jl")
