using DifferentialEquations
using CSV
using PyPlot
PyPlot.ioff()
using LaTeXStrings

function SARS_CoV2_replication(du,u,p,t)
    V_free, V_bound, V_endosome, gRNA_plus,
    DIP_free, DIP_bound, DIP_endosome, DIP_gRNA_plus,
    NSP,
    gRNA_minus, gRNA, DIP_gRNA_minus, DIP_gRNA,
    N, SP, N_gRNA, DIP_N_gRNA,
    V_assembled, V_released, DIP_assembled, DIP_released = u

    (k_bind, d_V, k_diss, k_fuse, k_uncoat, d_endosome,
     d_gRNA, k_transl, f_ORF1, d_NSP, k_tr_minus, K_NSP,
     d_gRNA_minus, k_tr_plus, k_complex, K_N,
     f_N, f_SP, d_N, d_SP,
     n_SP, n_N, K_V_rel, k_assembl,
     d_N_gRNA, k_release, d_assembled,
     d_V_DIP, d_endosome_DIP, d_gRNA_DIP,
     k_trans_minus, k_trans_plus,
     k_trans_minus_DIP, k_trans_plus_DIP,
     k_tr_minus_DIP, d_gRNA_minus_DIP,
     k_tr_plus_DIP, k_complex_DIP,
     n_N_DIP, k_assembl_DIP, n_SP_DIP, K_V_rel_DIP,
     d_N_gRNA_DIP, k_release_DIP, d_assembled_DIP) = p

    #dV_free
    du[1] = -k_bind*V_free - d_V*V_free + k_diss*V_bound
    #dV_bound
    du[2] = k_bind*V_free - (k_fuse + k_diss + d_V)*V_bound
    #dV_endosome
    du[3] = k_fuse*V_bound - (k_uncoat + d_endosome)*V_endosome
    #dgRNA_plus
    du[4] = k_uncoat*V_endosome - d_gRNA*gRNA_plus

    #dDIP_free
    du[5] = -k_bind*DIP_free - d_V_DIP*DIP_free + k_diss*DIP_bound
    #dDIP_bound
    du[6] = k_bind*DIP_free - (k_fuse + k_diss + d_V_DIP)*DIP_bound
    #dDIP_endosome
    du[7] = k_fuse*DIP_bound - (k_uncoat + d_endosome_DIP)*DIP_endosome
    #dDIP_gRNA_plus
    du[8] = k_uncoat*DIP_endosome - d_gRNA_DIP*DIP_gRNA_plus

    #dNSP
    du[9] = k_transl*f_ORF1*gRNA_plus - d_NSP*NSP -
           (k_trans_minus*gRNA_plus + k_trans_plus*gRNA_minus +
            k_trans_minus_DIP*DIP_gRNA_plus + k_trans_plus_DIP*DIP_gRNA_minus)*NSP

    f_RdRp = NSP/(NSP + K_NSP)
    f_complex = N/(N + K_N)
    f_assemb = SP/(SP + K_V_rel*n_SP)
    f_assemb_DIP = SP/(SP + K_V_rel_DIP*n_SP_DIP)

    #dgRNA_minus
    du[10] = k_tr_minus*gRNA_plus*f_RdRp - d_gRNA_minus*gRNA_minus
    #dgRNA
    du[11] = k_tr_plus*gRNA_minus*f_RdRp - (k_complex*f_complex + d_gRNA)*gRNA

    #dDIP_gRNA_minus
    du[12] = k_tr_minus_DIP*DIP_gRNA_plus*f_RdRp - d_gRNA_minus_DIP*DIP_gRNA_minus
    #dDIP_gRNA
    du[13] = k_tr_plus_DIP*DIP_gRNA_minus*f_RdRp - (k_complex_DIP*f_complex + d_gRNA_DIP)*DIP_gRNA

    #dN
    du[14] = k_transl*f_N*gRNA - k_complex*n_N*f_complex*gRNA - k_complex_DIP*n_N_DIP*f_complex*DIP_gRNA - d_N*N

    #dSP
    du[15] = k_transl*f_SP*gRNA - k_assembl*n_SP*f_assemb*N_gRNA - k_assembl_DIP*n_SP_DIP*f_assemb_DIP*DIP_N_gRNA - d_SP*SP

    #dN_gRNA
    du[16] = k_complex*f_complex*gRNA - (k_assembl*f_assemb + d_N_gRNA)*N_gRNA
    #dDIP_N_gRNA
    du[17] = k_complex_DIP*f_complex*DIP_gRNA - (k_assembl_DIP*f_assemb_DIP + d_N_gRNA_DIP)*DIP_N_gRNA

    #dV_assembled
    du[18] = k_assembl*f_assemb*N_gRNA - (k_release + d_assembled)*V_assembled
    #dV_released
    du[19] = k_release*V_assembled - d_V*V_released

    #dDIP_assembled
    du[20] = k_assembl_DIP*f_assemb_DIP*DIP_N_gRNA - (k_release_DIP + d_assembled_DIP)*DIP_assembled
    #dDIP_released
    du[21] = k_release_DIP*DIP_assembled - d_V_DIP*DIP_released
end

# parameters
pt = CSV.File("parameters.txt", header=false, delim=' ', ignorerepeated=true)
pnames = String.(pt.Column1)
p = Float64.(pt.Column2)

println("Obtaining the numerical solution of the model")

# initial condition
u0 = repeat([0.0],21)
V0 = 10.0
u0[1] = V0
u0[5] = V0

# solve the model
t_final = 24.0
tspan = (0.0, t_final)
prob = ODEProblem(SARS_CoV2_replication,u0,tspan,p)
sol = solve(prob,KenCarp5(),reltol=1e-16,abstol=1e-16)

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

function draw_sol(sol, filename)
    i2mm = 0.039370077777777776
    function mydraw(disp_vars, tspan; loc="upper left", ythresh=1.0, ythscale=0.4,
                                      lstyle="-", colors=nothing)
        if isnothing(colors)
            pl = ax.plot(sol.t, Array(sol')[:,disp_vars], linestyle=lstyle, linewidth=2.0)
        else
            for (var,col) in zip(disp_vars,colors)
                pl = ax.plot(sol.t, Array(sol')[:,var], linestyle=lstyle, color=col, linewidth=2.0)
            end
        end
        ax.set_xlim(tspan)
        ax.set_xlabel("t, hours")
        ax.set_ylim((0,1.15*ax.get_ylim()[2]))
        ax.set_yscale("symlog", linthresh=ythresh, linscale=ythscale, subs=2:9)
        ax.minorticks_on()
        ax.grid(which="major", linestyle="-", linewidth=0.1, color="lightgray")
        ax.tick_params(width=0.1, which="minor")
        ax.legend(loc=loc, labels=L[disp_vars])
        return pl
    end
    ##ax.yaxis.set_minor_locator(matplotlib.ticker.SymmetricalLogLocator(base=10.0,linthresh=ythresh))
    #ax.grid(which="minor", linestyle=":", linewidth=0.01)

    PyPlot.ioff()

    PyPlot.rc("text", usetex=false)
    PyPlot.rc("font", size=12)

    fig, axes = subplots(3,2,
         figsize=(210.0*i2mm,270.0*i2mm))
         #constrained_layout=true

    ax = axes[1,1]
    mydraw([1,2,5,6], [0,5], loc="upper right")
    ax.set_ylim((0,10))
    ax.set_yscale("linear")
    ax.set_ylabel("particles")
    ax.set_title("Binding and fusion")

    ax = axes[1,2]
    mydraw([3,4,7,8], [0,24], loc="upper right")
    ax.set_yscale("linear")
    ax.set_ylabel("numbers")
    ax.set_title("Endocytosis and uncoating")

    ax = axes[2,1]
    mydraw([9,10,12], [0,24], loc="upper left")
    ax.set_yscale("linear")
    ax.set_ylabel("molecules")
    ax.set_title(L"ORF1 translation and gRNA$_{(-)}$ synthesis")

    ax = axes[2,2]
    mydraw([11,13,14], [0,24], loc="lower right") #, ythresh=0.1
    ax.set_ylabel("molecules")
    ax.set_title("Transcription and translation")

    ax = axes[3,1]
    mydraw([15,16,17], [0,24], loc="lower right") #, ythresh=0.1
    ax.set_ylabel("molecules")
    ax.set_title("Translation and nucleocapsid formation")

    ax = axes[3,2]
    mydraw([18,19,20,21], [0,24], loc="upper left")
    #ax.set_yscale("linear")
    ax.set_ylabel("particles")
    ax.set_title("Assembly and release")

    plt.tight_layout()
    PyPlot.savefig(filename)
    plt.close(fig)
end
draw_sol(sol, "plot_model_det.pdf")

Vs = 1:20
DIPs = 0:20

Vwt = []
Vdip = []

for (V0, DIP0) in Base.product(Vs, DIPs)
    u0_ = Int64.(copy(u0))
    u0_[1] = V0
    u0_[5] = DIP0

    prob_ = remake(prob, u0=u0_)
    sol_ = solve(prob_,KenCarp5(),reltol=1e-16,abstol=1e-16)

    push!(Vwt, sol_(24, idxs=19))
    push!(Vdip, sol_(24, idxs=21))
end

Vv = [V0 for (V0, DIP0) in Base.product(Vs, DIPs)]
Dd = [DIP0 for (V0, DIP0) in Base.product(Vs, DIPs)]
CSV.write("result_det.csv", DataFrame(V0=Vv[:], DIP0=Dd[:], V_wt=Vwt, V_dip=Vdip))

