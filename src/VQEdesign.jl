using Yao, Yao.ConstGate, BitBasis
using Statistics, Plots
using KrylovKit: eigsolve


include("paracircuit.jl")
include("../model/AbstractModel.jl")
include("train.jl")

ising1d = Ising(6; periodic = false, h = 2, J_coupling = 0.8)
repeat = 1
ising1dfull = FullConnectStruct(nspins(ising1d), repeat)
ising1dfullcircuit = get_struct(ising1dfull, CNOT)
ising1dmachine = TrainMachine(hamiltonian(ising1d), ising1dfullcircuit)

machine1d, his1d = train(ising1dmachine, ising1d; maxiter=1000, optimizer = Optimise.ADAM(0.1))
Plots.plot(his1d)
paras = parameters(machine1d.circuit)
Plots.plot(paras)


ising2d = Ising(3, 3; periodic = false, h = 2, J_coupling = 0.8)
ising2dfull = FullConnectStruct(nspins(ising2d), repeat)
ising2dfullcircuit = get_struct(ising2dfull, CNOT)
ising2dmachine = TrainMachine(hamiltonian(ising2d), ising2dfullcircuit)

machine2d, his2d = train(ising2dmachine, ising2d; maxiter=1000, optimizer = Optimise.ADAM(0.1))
Plots.plot(his2d)
paras = parameters(machine2d.circuit)
Plots.plot(paras)



reg0 = reduce(⊗, rand_state(1) for i = 1:6)
para0 = copy(paras)
