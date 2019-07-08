export AbstractModel, IsingModel
export hamiltonian, ground_state, energy_current, get_Hsingle, get_Hcoupling, nspins

abstract type AbstractModel{D} end

nspins(model::AbstractModel) = prod(size(model)) #get spins count

"""
    energy_exact(tc, model::AbstractModel) -> Float64

Exact ground state energy.
"""
function energy_current(reg, hami)
    expect(hami, reg) |> real #get expectation value of hamiltonian
end

"""
    ground_state(model::AbstractModel) -> DefaultRegister

Get the EXACT ground state of a model.
"""
function ground_state(model::AbstractModel)
    # get the ground state
    hami = hamiltonian(model)
    E, v = eigsolve(mat(hami), 1, :SR)
    E[1], ArrayReg(v[1])
end

"""
    get_Hsingle(model::Ising) -> Vector

Get the weighted spin terms of a Ising model in the form Tuple(i, h_i).
"""
function get_Hsingle end

"""
    get_Hcoupling(model::Ising) -> Vector

Get the weighted spin-spin coupling terms of a Ising model in the form Tuple(i, j, J_ij).
"""
function get_Hcoupling end

include("IsingModel.jl")
