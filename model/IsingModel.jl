struct Ising{D} <: AbstractModel{D}
    size::NTuple{D, Int}
    periodic::Bool
    h::Float64
    J_coupling::Float64
    Ising(size::Int...; periodic::Bool, h::Real, J_coupling::Real) = new{length(size)}(size, periodic, Float64(h), Float64(J_coupling))
end

Base.size(model::Ising) = model.size

single_spin(nqubit::Int, i::Int) = put(nqubit, i => Z)
coupling_spins(nqubit::Int, i::Int, j::Int = i + 1) = put(nqubit, i => Z) * put(nqubit, j => Z)

function hamiltonian(model::Ising)
    Hsingle = Sum([x[2] * single_spin(nspins(model), x[1]) for x in get_Hsingle(model)]...)
    Hcoupling = Sum([x[3] * coupling_spins(nspins(model), x[1], x[2]) for x in get_Hcoupling(model)]...)
    Sum(Hsingle..., Hcoupling...)
end

function get_Hsingle(model::Ising{2})
    m, n = model.size
    cis = LinearIndices(model.size)
    H_single = Tuple{Int, Float64}[]
    for i = 1:m, j = 1:n
        push!(H_single, (cis[i, j], model.h))
    end
    H_single
end

get_Hsingle(model::Ising{1}) = [(i, model.h) for i in 1:model.size[1]]

function get_Hcoupling(model::Ising{2})
    m, n = model.size
    cis = LinearIndices(model.size)
    H_coupling = Tuple{Int, Int, Float64}[]
    for i = 1:m, j = 1:n
        (i != m || model.periodic) && push!(H_coupling, (cis[i, j], cis[i%m + 1, j], model.J_coupling))
        (j != n || model.periodic) && push!(H_coupling, (cis[i, j], cis[i, j%n + 1], model.J_coupling))
    end
    H_coupling
end

function get_Hcoupling(model::Ising{1})
    nqubit, = model.size
    [(i, i%nqubit + 1, model.J_coupling) for i in 1:(model.periodic ? nqubit : nqubit-1)]
end
