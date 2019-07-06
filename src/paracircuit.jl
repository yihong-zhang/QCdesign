export FullConnectStruct, TwoLevelStruct
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

@const_gate CZ = mat(control(2, 1, 2=>Z))

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
                chain(c.nqubit, paraEntangler(mode, c.nqubit, j, m) for m = j+1 : c.nqubit) |> add!
            end
        end
    elseif typeof(c) == TwoLevelStruct
        for i = 1:c.nqubit
            if i%2 == 1
                chain(c.nqubit, paraEntangler(mode, c.nqubit, l1, l1 + 1) for l1 = 1 : 2 : c.nqubit - 2 * (c.nqubit%2)) |> add!
            else
                chain(c.nqubit, paraEntangler(mode, c.nqubit, l2, l2 + 1) for l2 = 2: 2: c.nqubit - 2 * ((c.nqubit+1)%2)) |> add!
            end
        end
    else
        println("Datatype wrong!")
    end
    return circ
end
