include("stoch_model.jl")

# for histograms
n = 10^4 + 1

output_func(sol,i) = (sol[end][[11,13,19,21]], false)
ens_prob = EnsembleProblem(jump_prob, output_func=output_func)

@time sim = solve(ens_prob, FunctionMap(),
                  EnsembleThreads(),
                  saveat=24.0,
                  trajectories=n)

CSV.write("hist_stat.csv", DataFrame(gRNA_wt = [u[1] for u in sim.u],
                                     gRNA_dip =[u[2] for u in sim.u],
                                     V_released_wt = [u[3] for u in sim.u],
                                     V_released_dip =[u[4] for u in sim.u]))
