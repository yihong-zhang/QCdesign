export FullConnectStruct, TwoLevelStruct, NearestStruct
export CZ, paraEntangler, get_struct

abstract type AbsCircuitStructure end
struct FullConnectStruct <: AbsCircuitStructure
    nqubit::Int
    nrepeat::Int
end

struct TwoLevelStruct <: AbsCircuitStructure
    nqubit::Int
    nrepeat::Int
end

struct NearestStruct <: AbsCircuitStructure
    nqubit::Int
    nrepeat::Int
end

@const_gate CZ = mat(control(2, 1, 2=>Z))

rots(n::Int, loc::Int) = chain(n, put(n, loc => rot(G, 0.0)) for G in [X, Y,Z])

function paraEntangler(mode, nqubit::Int, loc1::Int, loc2::Int)
    put(nqubit, (loc1, loc2) => rot(mode, 0.0))
end

"""get full connect para-circuit structure"""
function get_struct(c::AbsCircuitStructure, mode)
    circ = chain(c.nqubit)
    add!(block) = push!(circ, block)
    if typeof(c) == FullConnectStruct
        for i = 1:c.nrepeat
            for j = 1:(c.nqubit - 1)
                for m = j+1 : c.nqubit
                    rots(c.nqubit, j) |> add!
                    rots(c.nqubit, m) |> add!
                    paraEntangler(mode, c.nqubit, j, m) |> add!
                end
            end
        end
    elseif typeof(c) == TwoLevelStruct
        for j = 1:c.nrepeat
#            chain(c.nqubit, put(c.nqubit, loc => chain(rot(Z, 0.0), rot(X, 0.0), rot(Z, 0.0))) for loc = 1:c.nqubit) |> add!
            for i = 1:c.nqubit
                if i%2 == 1
                    chain(c.nqubit, chain(rots(c.nqubit, l1), rots(c.nqubit, l1+1), paraEntangler(mode, c.nqubit, l1, l1 + 1)) for l1 = 1 : 2 : c.nqubit - 2 * (c.nqubit%2)) |> add!
                else
                    chain(c.nqubit, chain(rots(c.nqubit, l2), rots(c.nqubit, l2+1), paraEntangler(mode, c.nqubit, l2, l2 + 1)) for l2 = 2: 2: c.nqubit - 2 * ((c.nqubit+1)%2)) |> add!
                end
            end
        end
    elseif typeof(c) == NearestStruct
        for i = 1:c.nrepeat
            chain(c.nqubit, rots(c.nqubit, loc) for loc = 1:c.nqubit) |> add!
            for j = 1:(c.nqubit-1)
                chain(c.nqubit, paraEntangler(mode, c.nqubit, j, j + 1)) |> add!
            end
        end
    else
        println("Datatype wrong!")
    end
    return circ
end
