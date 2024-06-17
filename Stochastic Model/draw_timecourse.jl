using PyPlot
PyPlot.ioff()
using LaTeXStrings
using CSV, DataFrames
using Measures

data = CSV.read("kinetics_means.csv", DataFrame)

# draw numerical solution
sL = ["V_{free}^{wt}", "V_{bound}^{wt}", "V_{endosome}^{wt}", "gRNA_{(+)}^{wt}",
      "V_{free}^{dip}", "V_{bound}^{dip}",
      "V_{endosome}^{dip}", "gRNA_{(+)}^{dip}",
      "NSP", "gRNA_{(-)}^{wt}", "gRNA^{wt}", "gRNA_{(-)}^{dip}", "gRNA^{dip}",
      "N", "SP", "N-gRNA^{wt}", "N-gRNA^{dip}",
      "V_{assembled}^{wt}", "V_{released}^{wt}",
      "V_{assembled}^{dip}", "V_{released}^{dip}"]

L = [latexstring("[{\\rm $s}]") for s in sL]
L[16] = L"$[{\rm N}$-${\rm gRNA}^{\rm wt}]$"
L[17] = L"$[{\rm N}$-${\rm gRNA}^{\rm dip}]$"


function draw_sol(data, filename)
    i2mm = 0.039370077777777776
    function mydraw(disp_vars, tspan; loc="upper left", ythresh=1.0, ythscale=0.4,
                                      lstyle="-", colors=nothing)
        if isnothing(colors)
              pl = ax.plot(data.t, Array(data[:,disp_vars]), linestyle=lstyle, linewidth=2.0)
        else
            for (var,col) in zip(disp_vars,colors)
                  pl = ax.plot(data.t, Array(data[:,var]), linestyle=lstyle, color=col, linewidth=2.0)
            end
        end
        ax.set_xlim(tspan)
        ax.set_xlabel("t, hours")
        ax.set_ylim((0,1.15*ax.get_ylim()[2]))
        ax.set_yscale("symlog", linthresh=ythresh, linscale=ythscale, subs=2:9)
        ax.minorticks_on()
        ax.grid(which="major", linestyle="-", linewidth=0.1, color="lightgray")
        ax.tick_params(width=0.1, which="minor")
        ax.legend(loc=loc, labels=L[disp_vars.-1])
        return pl
    end
    #ax.yaxis.set_minor_locator(matplotlib.ticker.SymmetricalLogLocator(base=10.0,linthresh=ythresh))
    #ax.grid(which="minor", linestyle=":", linewidth=0.01)

    PyPlot.ioff()

    PyPlot.rc("text", usetex=false)
    PyPlot.rc("font", size=12)

    #rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    #rcParams["font.size"] = 12
    #rcParams["text.usetex"] = false

    fig, axes = subplots(3,2,
         figsize=(210.0*i2mm,270.0*i2mm))
         #constrained_layout=true

    ax = axes[1,1]
    mydraw([1,2,5,6].+1, [0,5], loc="upper right")
    ax.set_ylim((0,10))
    ax.set_yscale("linear")
    ax.set_ylabel("particles")
    ax.set_title("Binding and fusion")

    ax = axes[1,2]
    mydraw([3,4,7,8].+1, [0,24], loc="upper right")
    ax.set_yscale("linear")
    ax.set_ylabel("numbers")
    ax.set_title("Endocytosis and uncoating")

    ax = axes[2,1]
    mydraw([9,10,12].+1, [0,24], loc="upper left")
    ax.set_yscale("linear")
    ax.set_ylabel("molecules")
    ax.set_title(L"ORF1 translation and gRNA$_{(-)}$ synthesis")

    ax = axes[2,2]
    mydraw([11,13,14].+1, [0,24], loc="lower right") #, ythresh=0.1
    ax.set_ylabel("molecules")
    ax.set_title("Transcription and translation")

    ax = axes[3,1]
    mydraw([15,16,17].+1, [0,24], loc="lower right") #, ythresh=0.1
    ax.set_ylabel("molecules")
    ax.set_title("Translation and nucleocapsid formation")

    ax = axes[3,2]
    mydraw([18,19,20,21].+1, [0,24], loc="upper left")
    #ax.set_yscale("linear")
    ax.set_ylabel("particles")
    ax.set_title("Assembly and release")

    plt.tight_layout()
    PyPlot.savefig(filename)
    plt.close(fig)
end
draw_sol(data, "plot_model_means.pdf")

data = CSV.read("kinetics_medians.csv", DataFrame)
draw_sol(data, "plot_model_medians.pdf")

function draw_sol_full(datas, filename)

    data_means, data_medians, data_q25, data_q75, data_q5, data_q95 = datas

    i2mm = 0.039370077777777776
    function mydraw(disp_vars, tspan, type; loc="upper left",
                ythresh=1.0, ythscale=0.4,
                lstyle="-", colors=nothing, lw=2.0)
        colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
        if type == :mean
              data = data_means
              lstyle="--"
        elseif type == :median
              data = data_medians
        elseif type == :q2575
              data_low = data_q25
              data_high = data_q75
              lw=0
        elseif type == :q0595
              data_low = data_q5
              data_high = data_q95
              lw=0
        end
        if type in [:mean, :median]
              if isnothing(colors)
                    pl = ax.plot(data.t, Array(data[:,disp_vars]),
                                 linestyle=lstyle, linewidth=lw)
              else
                  for (var,col) in zip(disp_vars,colors)
                        pl = ax.plot(data.t, Array(data[:,var]),
                                     linestyle=lstyle, color=col,
                                     linewidth=lw)
                  end
              end
        else
              if isnothing(colors)
                    pl = ax.fill_between(data_low.t,
                                        Array(data_low[:,disp_vars]),
                                        Array(data_high[:,disp_vars]),
                                 alpha=0.5,
                                 linestyle=lstyle, linewidth=lw)
              else
                  for (var,col) in zip(disp_vars,colors)
                        pl = ax.fill_between(data_low.t,
                                             Array(data_low[:,var]),
                                             Array(data_high[:,var]),
                                     alpha=0.5,
                                     linestyle=lstyle, color=col,
                                     linewidth=lw)
                  end
              end
        end
        ax.set_xlim(tspan)
        ax.set_xlabel("t, hours")
        ax.set_ylim((0,1.15*ax.get_ylim()[2]))
        ax.set_yscale("symlog", linthresh=ythresh, linscale=ythscale, subs=2:9)
        ax.minorticks_on()
        ax.grid(which="major", linestyle="-", linewidth=0.1, color="lightgray", alpha=0.5)
        ax.tick_params(width=0.1, which="minor")
        ax.legend(loc=loc, labels=L[disp_vars.-1])
        return pl
    end
    #ax.yaxis.set_minor_locator(matplotlib.ticker.SymmetricalLogLocator(base=10.0,linthresh=ythresh))
    #ax.grid(which="minor", linestyle=":", linewidth=0.01)

    PyPlot.ioff()

    PyPlot.rc("text", usetex=false)
    PyPlot.rc("font", size=12)

    #rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    #rcParams["font.size"] = 12
    #rcParams["text.usetex"] = false

    fig, axes = subplots(3,2,
         figsize=(210.0*i2mm,270.0*i2mm))
         #constrained_layout=true

    ax = axes[1,1]
    mydraw([1,2,5,6].+1, [0,5], :median, loc="upper right")
    mydraw([1,2,5,6].+1, [0,5], :mean, loc="upper right")
    mydraw([1,2,5,6].+1, [0,5], :q2575, loc="upper right")
    ax.set_ylim((0,10))
    ax.set_yscale("linear")
    ax.set_ylabel("particles")
    ax.set_title("Binding and fusion")

    ax = axes[1,2]
    mydraw([3,4,7,8].+1, [0,24], :median, loc="upper right")
    mydraw([3,4,7,8].+1, [0,24], :mean, loc="upper right")
    mydraw([3,4,7,8].+1, [0,24], :q2575, loc="upper right")
    ax.set_yscale("linear")
    ax.set_ylabel("numbers")
    ax.set_title("Endocytosis and uncoating")

    ax = axes[2,1]
    mydraw([9,10,12].+1, [0,24], :median, loc="upper right")
    mydraw([9,10,12].+1, [0,24], :mean, loc="upper right")
    mydraw([9,10,12].+1, [0,24], :q2575, loc="upper right")
    ax.set_yscale("linear")
    ax.set_ylabel("molecules")
    ax.set_title(L"ORF1 translation and gRNA$_{(-)}$ synthesis")

    ax = axes[2,2]
    mydraw([11,13,14].+1, [0,24], :median, loc="lower right") #, ythresh=0.1
    mydraw([11,13,14].+1, [0,24], :mean, loc="lower right") #, ythresh=0.1
    mydraw([11,13,14].+1, [0,24], :q2575, loc="lower right") #, ythresh=0.1
    ax.set_ylabel("molecules")
    ax.set_title("Transcription and translation")

    ax = axes[3,1]
    mydraw([15,16,17].+1, [0,24], :median, loc="lower right") #, ythresh=0.1
    mydraw([15,16,17].+1, [0,24], :mean, loc="lower right") #, ythresh=0.1
    mydraw([15,16,17].+1, [0,24], :q2575, loc="lower right") #, ythresh=0.1
    ax.set_ylabel("molecules")
    ax.set_title("Translation and nucleocapsid formation")

    ax = axes[3,2]
    mydraw([18,19,20,21].+1, [0,24], :median, loc="upper left")
    mydraw([18,19,20,21].+1, [0,24], :mean, loc="upper left")
    mydraw([18,19,20,21].+1, [0,24], :q2575, loc="upper left")
    #ax.set_yscale("linear")
    ax.set_ylabel("particles")
    ax.set_title("Assembly and release")

    plt.tight_layout()
    PyPlot.savefig(filename)
    plt.close(fig)
end

data_means = CSV.read("kinetics_means.csv", DataFrame)
data_medians = CSV.read("kinetics_medians.csv", DataFrame)
data_q25 = CSV.read("kinetics_q25.csv", DataFrame)
data_q75 = CSV.read("kinetics_q75.csv", DataFrame)
data_q5 = CSV.read("kinetics_q5.csv", DataFrame)
data_q95 = CSV.read("kinetics_q95.csv", DataFrame)

datas = [data_means, data_medians, data_q25, data_q75, data_q5, data_q95]

draw_sol_full(datas, "plot_model_full.pdf")
