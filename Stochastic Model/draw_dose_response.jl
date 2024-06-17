using Plots
using CSV, DataFrames
using LaTeXStrings
using Measures

data = CSV.read("result.csv", DataFrame)

plot(legend=:bottomright)
for DIP0 in 0:10
    df = filter(x -> x.DIP0 == DIP0, data)
    plot!(df.V0, df.p,
          label=L"V_{free}^{dip}(0) = " * "$DIP0",
          marker=:circle,
          markersize=1.5)
    xaxis!(xticks=df.V0)
    xlabel!(L"[V_{free}^{wt}](0)")
    ylabel!(L"p([V_{released}^{wt}(24)] > 0)")
end
title!("Probability of productive infection")
savefig("result_p.pdf")

# p1 = plot(legend=nothing)
# for DIP0 in 0:10
#     df = filter(x -> x.DIP0 == DIP0, data)
#     plot!(df.V0, df.V_mean,
#           label=L"V_{free}^{dip}(0) = " * "$DIP0",
#           marker=:circle,
#           markersize=2)
#     xaxis!(xticks=df.V0)
#     xlabel!(L"[V_{free}^{wt}](0)")
#     ylabel!(L"$[V_{released}^{wt}](24),$ mean value")
# end
# #title!("Mean number of released at 24 h.p.i. WT virions")
# #savefig("result_means.pdf")
# 
# p2 = plot(legend=:topleft)
# for DIP0 in 0:10
#     df = filter(x -> x.DIP0 == DIP0, data)
#     plot!(df.V0, df.V_median,
#           label=L"V_{free}^{dip}(0) = " * "$DIP0",
#           marker=:circle,
#           markersize=2)
#     xaxis!(xticks=df.V0)
#     xlabel!(L"[V_{free}^{wt}](0)")
#     ylabel!(L"$[V_{released}^{wt}](24),$ median value")
# end
# #title!("Median number of released at 24 h.p.i. WT virions")
# #savefig("result_medians.pdf")
# 
# p = plot(p1, p2,
#          plot_title = "Mean and median numbers of released at 24 h.p.i. WT virions",
#          plot_titlefontsize=10,
#          margins=2mm)
# plot!(size=(212*3.3,120*3.3))
# savefig("result_means_medians_wt.pdf")
# 
# p1 = plot(legend=nothing)
# for DIP0 in 0:10
#     df = filter(x -> x.DIP0 == DIP0, data)
#     plot!(df.V0, df.DIP_mean,
#           label=L"V_{free}^{dip}(0) = " * "$DIP0",
#           marker=:circle,
#           markersize=2)
#     xaxis!(xticks=df.V0)
#     xlabel!(L"[V_{free}^{wt}](0)")
#     ylabel!(L"$[V_{released}^{dip}](24),$ mean value")
# end
# #title!("Mean number of released at 24 h.p.i. WT virions")
# #savefig("result_means.pdf")
# 
# p2 = plot(legend=:topleft)
# for DIP0 in 0:10
#     df = filter(x -> x.DIP0 == DIP0, data)
#     plot!(df.V0, df.DIP_median,
#           label=L"V_{free}^{dip}(0) = " * "$DIP0",
#           marker=:circle,
#           markersize=2)
#     xaxis!(xticks=df.V0)
#     xlabel!(L"[V_{free}^{wt}](0)")
#     ylabel!(L"$[V_{released}^{dip}](24),$ median value")
# end
# #title!("Median number of released at 24 h.p.i. WT virions")
# #savefig("result_medians.pdf")
# 
# p = plot(p1, p2,
#          plot_title = "Mean and median numbers of released at 24 h.p.i. DIPs",
#          plot_titlefontsize=10,
#          margins=2mm)
# plot!(size=(212*3.3,120*3.3))
# savefig("result_means_medians_dip.pdf")
# 
# p1 = plot(legend=nothing)
# for V0 in 1:10
#     df = filter(x -> x.V0 == V0, data)
#     plot!(df.DIP0, df.DIP_mean,
#           label=L"V_{free}^{wt}(0) = " * "$V0",
#           marker=:circle,
#           markersize=2)
#     xaxis!(xticks=df.DIP0)
#     xlabel!(L"[V_{free}^{dip}](0)")
#     ylabel!(L"$[V_{released}^{dip}](24),$ mean value")
# end
# 
# p2 = plot(legend=:topleft)
# for V0 in 1:10
#     df = filter(x -> x.V0 == V0, data)
#     plot!(df.DIP0, df.DIP_median,
#           label=L"V_{free}^{wt}(0) = " * "$V0",
#           marker=:circle,
#           markersize=2)
#     xaxis!(xticks=df.DIP0)
#     xlabel!(L"[V_{free}^{dip}](0)")
#     ylabel!(L"$[V_{released}^{dip}](24),$ median value")
# end
# 
# p = plot(p1, p2,
#          plot_title = "Mean and median numbers of released at 24 h.p.i. DIPs",
#          plot_titlefontsize=10,
#          margins=2mm)
# plot!(size=(212*3.3,120*3.3))
# savefig("result_means_medians_dip_v2.pdf")

# heatmaps

clims = extrema([data.V_mean; data.DIP_mean])
p1 = heatmap(0:10, 1:10, reshape(data.V_mean, 10, 11),
             clims=clims, colorbar=false, c=:viridis)
title!("Total WT released (mean value at t=24 h)", titlefontsize=8)
xlabel!(L"[V_{free}^{dip}](0)", xticks=0:10)
ylabel!(L"[V_{free}^{wt}](0)", yticks=1:10)
p2 = heatmap(0:10, 1:10, reshape(data.DIP_mean, 10, 11),
             clims=clims, colorbar=false, c=:viridis)
title!("Total DIP released (mean value at t=24 h)",
       titlefontsize=8)
xlabel!(L"[V_{free}^{dip}](0)", xticks=0:10)
ylabel!(L"[V_{free}^{wt}](0)", yticks=1:10)
pc = scatter([0,0], [0,1], zcolor=[0,3], clims=clims,
             xlims=(1,1.1), label="",
             c=:viridis, framestyle=:none)
l = @layout [grid(1, 2) a{0.035w}]
p = plot(p1, p2, pc, layout=l)
savefig("result_heatmap_means.pdf")

clims = extrema([data.V_median; data.DIP_median])
p1 = heatmap(0:10, 1:10, reshape(data.V_median, 10, 11),
             clims=clims, colorbar=false, c=:viridis)
title!("Total WT released (median value at t=24 h)", titlefontsize=8)
xlabel!(L"[V_{free}^{dip}](0)", xticks=0:10)
ylabel!(L"[V_{free}^{wt}](0)", yticks=1:10)
p2 = heatmap(0:10, 1:10, reshape(data.DIP_median, 10, 11),
             clims=clims, colorbar=false, c=:viridis)
title!("Total DIP released (median value at t=24 h)",
       titlefontsize=8)
xlabel!(L"[V_{free}^{dip}](0)", xticks=0:10)
ylabel!(L"[V_{free}^{wt}](0)", yticks=1:10)
pc = scatter([0,0], [0,1], zcolor=[0,3], clims=clims,
             xlims=(1,1.1), label="",
             c=:viridis, framestyle=:none)
l = @layout [grid(1, 2) a{0.035w}]
p = plot(p1, p2, pc, layout=l)
savefig("result_heatmap_medians.pdf")

# deterministic heatmaps
data = CSV.read("result_det.csv", DataFrame)
Vs = 1:10
DIPs = 0:10

clims = extrema([data.V_wt; data.V_dip])
p1 = heatmap(DIPs, Vs, reshape(data.V_wt, length(Vs), length(DIPs)),
             clims=clims, colorbar=false, c=:viridis)
title!("Total WT released (determ. sol. at t=24 h)", titlefontsize=8)
xlabel!(L"[V_{free}^{dip}](0)", xticks=DIPs) #0:2:20)
ylabel!(L"[V_{free}^{wt}](0)", yticks=Vs)
p2 = heatmap(DIPs, Vs, reshape(data.V_dip, length(Vs), length(DIPs)),
             clims=clims, colorbar=false, c=:viridis)
title!("Total DIP released (determ. sol. at t=24 h)",
       titlefontsize=8)
xlabel!(L"[V_{free}^{dip}](0)", xticks=DIPs) #0:2:20)
ylabel!(L"[V_{free}^{wt}](0)", yticks=Vs)
pc = scatter([0,0], [0,1], zcolor=[0,3], clims=clims,
             xlims=(1,1.1), label="",
             c=:viridis, framestyle=:none)
l = @layout [grid(1, 2) a{0.035w}]
p = plot(p1, p2, pc, layout=l)
savefig("result_heatmap_det.pdf")
