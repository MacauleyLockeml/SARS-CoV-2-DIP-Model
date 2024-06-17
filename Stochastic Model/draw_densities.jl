using Plots
using StatsPlots
using CSV, DataFrames
using LaTeXStrings
using Measures

data = CSV.read("hist_stat.csv", DataFrame)

p1 = histogram(data.gRNA_wt, label=L"[gRNA^{wt}]", alpha=0.5, normalize=:pdf)
histogram!(data.gRNA_dip, label=L"[gRNA^{dip}]", alpha=0.5, normalize=:pdf)
xaxis!(xlims=(0,4e4), xformatter=x -> string(round(Int, x/1e4)))
annotate!([(4e4 * 0.9, -6e-5, Plots.text("x10⁴", 8, :black, :center))])
ylabel!(L"p")
xlabel!("molecules")

p2 = histogram(data.V_released_wt, label=L"[V_{released}^{wt}]", alpha=0.5, normalize=:pdf)
histogram!(data.V_released_dip, label=L"[V_{released}^{dip}]", alpha=0.5, normalize=:pdf)
xaxis!(xlims=(0,4e2))
ylabel!(L"p")
xlabel!("virions")

p = plot(p1, p2,
        margins=3mm,
        labelfontsize=8,
        plot_title="Probability densities of viral gRNAs and released virions at 24 h.p.i.",
        plot_titlefontsize=10)

savefig("hist.pdf")


p1 = density(data.gRNA_wt, label=L"[gRNA^{wt}]")
density!(data.gRNA_dip, label=L"[gRNA^{dip}]")
xaxis!(xlims=(0,4e4), xformatter=x -> string(round(Int, x/1e4)))
annotate!([(4e4 * 0.9, -6e-5, Plots.text("x10⁴", 8, :black, :center))])
ylabel!(L"p")
xlabel!("molecules")

p2 = density(data.V_released_wt, label=L"[V_{released}^{wt}]")
density!(data.V_released_dip, label=L"[V_{released}^{dip}]")
xaxis!(xlims=(0,4e2))
xlabel!("virions")

#xaxis!(xscale=:log10, xlims=(1,1e5))

p = plot(p1, p2,
        margins=3mm,
        labelfontsize=8,
        plot_title="Probability densities of viral gRNAs and released virions at 24 h.p.i.",
        plot_titlefontsize=10)

savefig("densities.pdf")
