include("stoch_model.jl")

# for dose-response
n = 10^3 + 1 # odd number for integer medians

# save V_released_wt and V_released_dip at t=24h
output_func(sol,i) = (sol[end][[19,21]], false)

# init conditions
Vs = 1:10
DIPs = 0:10

probabilities = []
V_means = []
V_medians = []
DIP_means = []
DIP_medians = []

using Statistics

for (V0, DIP0) in Base.product(Vs, DIPs)
    u0_ = Int64.(copy(u0))
    u0_[1] = V0
    u0_[5] = DIP0

    jump_prob_ = remake(jump_prob, u0=u0_)
    ens_prob_ = EnsembleProblem(jump_prob_,
                                output_func=output_func)

    time_elapsed = @elapsed sim_ = solve(ens_prob_, FunctionMap(),
                                         EnsembleThreads(),
                                         saveat=24.0,
                                         trajectories=n)

    simV = [u[1] for u in sim_.u]
    simDIP = [u[2] for u in sim_.u]

    p_prod = 1.0 - sum(simV .== 0)/n
    V_mean = mean(simV)
    V_median = median(simV)
    DIP_mean = mean(simDIP)
    DIP_median = median(simDIP)

    push!(probabilities, p_prod)
    push!(V_means, V_mean)
    push!(V_medians, V_median)
    push!(DIP_means, DIP_mean)
    push!(DIP_medians, DIP_median)

    println("V0 = $V0, DIP0 = $DIP0 (elapsed $(time_elapsed)s):")
    println("p = $p_prod, V_mean = $V_mean, V_median = $V_median, DIP_mean = $DIP_mean, DIP_median = $DIP_median")
end

Vv = [V0 for (V0, DIP0) in Base.product(Vs, DIPs)]
Dd = [DIP0 for (V0, DIP0) in Base.product(Vs, DIPs)]
CSV.write("result.csv", DataFrame(V0=Vv[:], DIP0=Dd[:], p=probabilities,
                                  V_mean=V_means, V_median=V_medians,
                                  DIP_mean=DIP_means, DIP_median=DIP_medians))
