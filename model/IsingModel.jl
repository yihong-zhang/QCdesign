struct IsingModel{D} <: AbstractModel{D}
    size::NTuple{D, Int}
    periodic::Bool
    J_coupling::Float64
    h::Float64
    IsingModel(size::Int...; periodic::Bool) = new{length(size)}(size, periodic)
end

Base.size(model::IsingModel) = model.size

single_spin(nqubit::Int, i::Int) = put(nqubit, i => Z)
coupling_spins(nqubit::Int, i::Int, j::Int = i + 1) = put(nqubit, i => Z) * put(nqubit, j => Z)

function hamiltonian(model::IsingModel)
    Sum([x[3] * coupling_spins(nspins(model), x[1], x[2]) for x in get_bonds(model)]...)
end

function get_bonds(model::IsingModel{2})
    m, n = model.size
    cis = LinearIndices(model.size)
    bonds = Tuple{Int, Int, Float64}[]
    for i = 1:m, j = 1:n
        (i != m || model.periodic) && push!(bonds, (cis[i,j], cis[i%m + 1, j], 1.0))
        (j != n || model.periodic) && push!(bonds, (cis[i,j], cis[i, j%n + 1], 1.0))
    end
    bonds
end

function get_bonds(model::IsingModel{1})
    nqubit, = model.size
    [(i, i%nqubit + 1, 1.0) for i in 1:(model.periodic ? nqubit : nqubit-1)]
end

"""
    energy(config, model::AbstractHeisenberg; nbatch) -> Float64

Ground state energy by sampling Quantum circuit.
The hamiltonian is limited to Heisenberg and J1J2 Type.
"""
function energy(qpeps::QPEPSMachine, model::AbstractHeisenberg)
    local eng = 0.0
    for basis in [X, Y, Z]
        mres = gensample(qpeps, basis)
        for (i,j,w) in get_bonds(model)
            eng += w*(1-2*mean(mres[i] .⊻ mres[j]))
        end
    end
    eng/=4
end
