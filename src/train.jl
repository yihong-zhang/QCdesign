export TrainMachine, get_gradients, train
using Flux: Optimise

struct TrainMachine
    hami #hamiltonian of the model
    circuit
end

function get_gradients(machine::TrainMachine, reg0::AbstractRegister)
    machine_rotor = collect_blocks(RotationGate, machine.circuit)
    return map(machine_rotor) do r
        r.theta += π/2
        reg = copy(reg0) |> machine.circuit
        E₊ = energy_current(reg, machine.hami)
        r.theta -= π
        reg = copy(reg0) |> machine.circuit
        E₋ = energy_current(reg, machine.hami)
        r.theta += π/2
        0.5(E₊ - E₋)
    end
end

function train(machine::TrainMachine, model::AbstractModel; maxiter = 200, optimizer = Optimise.ADAM(0.1))
    @assert nqubits(machine.hami) == nspins(model)
    reg0 = reduce(⊗, rand_state(1) for i = 1:nspins(model))
    dispatch!(machine.circuit, :random)
    println("E0 = $(energy_current(reg0, machine.hami))")
    flush(stdout)

    history = Float64[]
    params = parameters(machine.circuit)
    println("Number of parameters is $(length(params))")
    for i in 1:maxiter
        grad = get_gradients(machine, reg0)
        Optimise.update!(optimizer, params, grad)
        dispatch!(machine.circuit, params)
        reg = copy(reg0)
        push!(history, energy_current(reg |> machine.circuit, machine.hami))
        println("Iter $i, E = $(history[end])")
        flush(stdout)
    end
    reg0, machine, history
end
