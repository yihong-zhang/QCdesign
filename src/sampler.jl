export CircuitRunTime, basis_rotor

struct CircuitRunTime
    circuit
    rotblocks
    mbasis
    mblocks
end


basis_rotor(::ZGate) = I2Gate()
basis_rotor(::XGate) = Ry(-0.5π)
basis_rotor(::YGate) = Rx(0.5π)

basis_rotor(basis::PauliGate, nqubit, locs) = repeat(nqubit, basis_rotor(basis), locs)
