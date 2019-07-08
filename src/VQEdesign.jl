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

ising1dnear = NearestStruct(nspins(ising1d), repeat)
ising1dnearcircuit = get_struct(ising1dnear, CNOT)
ising1dmachine = TrainMachine(hamiltonian(ising1d), ising1dnearcircuit)


reg0, machine1d, his1d = train(ising1dmachine, ising1d; maxiter=500, optimizer = Optimise.ADAM(0.1))
Plots.plot(his1d)
#paras = parameters(machine1d.circuit)
#Plots.plot(paras)

Eg, vg = ground_state(ising1d)
reg = copy(reg0) |> machine1d.circuit
fidelity(reg, vg)




ising2d = Ising(3, 3; periodic = false, h = 2, J_coupling = 0.8)
repeat = 1
ising2dstruct = FullConnectStruct(nspins(ising2d), repeat)
ising2dcircuit = get_struct(ising2dstruct, CNOT)
ising2dmachine = TrainMachine(hamiltonian(ising2d), ising2dcircuit)
reg0, machine2d, his2d = train(ising2dmachine, ising2d; maxiter=200, optimizer = Optimise.ADAM(0.1))
Plots.plot(his2d)
Eg, vg = ground_state(ising2d)
reg = copy(reg0) |> machine2d.circuit
fidelity(reg, vg)



ising2dfull = FullConnectStruct(nspins(ising2d), repeat)
ising2dfullcircuit = get_struct(ising2dfull, CNOT)
ising2dmachine = TrainMachine(hamiltonian(ising2d), ising2dfullcircuit)

machine2d, his2d = train(ising2dmachine, ising2d; maxiter=1000, optimizer = Optimise.ADAM(0.1))
Plots.plot(his2d)
paras = parameters(machine2d.circuit)
Plots.plot(paras)


"""
Heisenberg model

"""

heisenberg2d = Heisenberg(3, 3; periodic = false, J_coupling = 0.8)
repeat = 1
heisenberg2dstruct = TwoLevelStruct(nspins(heisenberg2d), repeat)
heisenberg2dcircuit = get_struct(heisenberg2dstruct, CNOT)
heisenberg2dmachine = TrainMachine(hamiltonian(heisenberg2d), heisenberg2dcircuit)
reg0, machine2d, his2d = train(heisenberg2dmachine, ising2d; maxiter=5000, optimizer = Optimise.ADAM(0.1))
Plots.plot(his2d)
Eg, vg = ground_state(heisenberg2d)
fidelity(copy(reg0) |> machine2d.circuit, vg)
