struct Heisenberg{D} <: AbstractModel{D}
    size::NTuple{D, Int}
    periodic::Bool
    J_coupling::Float64
    Heisenberg(size::Int...; periodic::Bool, J_coupling::Real) = new{length(size)}(size, periodic, Float64(J_coupling))
end

Base.size(model::Heisenberg) = model.size
Heisenberg_ij(nqubit::Int, i::Int, j::Int = i+1) = put(nqubit, i=>X) * put(nqubit, j=>X) + put(nqubit, i=>Y) * put(nqubit, j=>Y) + put(nqubit, i=>Z) * put(nqubit, j=>Z)
#const Heisenberg_term = repeat(2, X, 1:2) + repeat(2, Y, 1:2) + repeat(2, Z, 1:2)



function get_HeisenbergCoupling(model::Heisenberg{2})
    m, n = model.size
    cis = LinearIndices(model.size)
    HeisenbergCoupling = Tuple{Int, Int, Float64}[]
    for i=1:m, j=1:n
        (i!=m || model.periodic) && push!(HeisenbergCoupling, (cis[i,j], cis[i%m + 1, j], model.J_coupling))
        (j!=n || model.periodic) && push!(HeisenbergCoupling, (cis[i,j], cis[i, j%n + 1], model.J_coupling))
    end
    HeisenbergCoupling
end

function get_HeisenbergCoupling(model::Heisenberg{1})
    nqubit, = model.size
    [(i, i%nqubit + 1, model.J_coupling) for i in 1:(model.periodic ? nqubit : nqubit-1)]
end

hamiltonian(model::Heisenberg) = Sum([x[3] * Heisenberg_ij(nspins(model), x[1], x[2]) for x in get_HeisenbergCoupling(model)]...)
