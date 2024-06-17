using DifferentialEquations
using CSV

jumps = []
vartojumps = Vector{Vector{Int64}}()
jumptovars = Vector{Vector{Int64}}()

push!(vartojumps, [1,2])
push!(vartojumps, [3,4,5])
push!(vartojumps, [6,7])
push!(vartojumps, [8,17,19,23])
push!(vartojumps, [9,10])
push!(vartojumps, [11,12,13])
push!(vartojumps, [14,15])
push!(vartojumps, [16,21,28])
push!(vartojumps, [18,19,20,21,22,23,25,28,30])
push!(vartojumps, [20,24,25])
push!(vartojumps, [26,27,33,35])
push!(vartojumps, [22,29,30])
push!(vartojumps, [31,32])
push!(vartojumps, [27,32,34])
push!(vartojumps, [36,37,38])
push!(vartojumps, [37,39])
push!(vartojumps, [38,40])
push!(vartojumps, [41,42])
push!(vartojumps, [43])
push!(vartojumps, [44,45])
push!(vartojumps, [46])

#    #dV_free
#    du[1] = -k_bind*V_free - d_V*V_free + k_diss*V_bound
function rate1(u,p,t)
    #k_bind*V_free
    p[1]*u[1]
end
function affect1!(integrator)
    integrator.u[1] -= 1
    integrator.u[2] += 1
    nothing
end
push!(jumps, ConstantRateJump(rate1, affect1!))
push!(jumptovars, [1,2])

function rate2(u,p,t)
    #d_V*V_free
    p[2]*u[1]
end
function affect2!(integrator)
    integrator.u[1] -= 1
    nothing
end
push!(jumps, ConstantRateJump(rate2, affect2!))
push!(jumptovars, [1])

function rate3(u,p,t)
    #k_diss*V_bound
    p[3]*u[2]
end
function affect3!(integrator)
    integrator.u[1] += 1
    integrator.u[2] -= 1
    nothing
end
push!(jumps, ConstantRateJump(rate3, affect3!))
push!(jumptovars, [1,2])

#    #dV_bound
#    du[2] = k_bind*V_free - (k_fuse + k_diss + d_V)*V_bound
function rate4(u,p,t)
    #k_fuse*V_bound
    p[4]*u[2]
end
function affect4!(integrator)
    integrator.u[2] -= 1
    integrator.u[3] += 1
    nothing
end
push!(jumps, ConstantRateJump(rate4, affect4!))
push!(jumptovars, [2,3])

function rate5(u,p,t)
    #d_V*V_bound
    p[2]*u[2]
end
function affect5!(integrator)
    integrator.u[2] -= 1
    nothing
end
push!(jumps, ConstantRateJump(rate5, affect5!))
push!(jumptovars, [2])

#    #dV_endosome
#    du[3] = k_fuse*V_bound - (k_uncoat + d_endosome)*V_endosome
function rate6(u,p,t)
    #k_uncoat*V_endosome
    p[5]*u[3]
end
function affect6!(integrator)
    integrator.u[3] -= 1
    integrator.u[4] += 1
    nothing
end
push!(jumps, ConstantRateJump(rate6, affect6!))
push!(jumptovars, [3,4])

function rate7(u,p,t)
    #d_endosome*V_endosome
    p[6]*u[3]
end
function affect7!(integrator)
    integrator.u[3] -= 1
    nothing
end
push!(jumps, ConstantRateJump(rate7, affect7!))
push!(jumptovars, [3])

#    #dgRNA_plus
#    du[4] = k_uncoat*V_endosome - d_gRNA*gRNA_plus
function rate8(u,p,t)
    #d_gRNA*gRNA_plus
    p[7]*u[4]
end
function affect8!(integrator)
    integrator.u[4] -= 1
    nothing
end
push!(jumps, ConstantRateJump(rate8, affect8!))
push!(jumptovars, [4])

#    #dDIP_free
#    du[5] = -k_bind*DIP_free - d_V_DIP*DIP_free + k_diss*DIP_bound
function rate9(u,p,t)
    #k_bind*DIP_free
    p[1]*u[5]
end
function affect9!(integrator)
    integrator.u[5] -= 1
    integrator.u[6] += 1
    nothing
end
push!(jumps, ConstantRateJump(rate9, affect9!))
push!(jumptovars, [5,6])

function rate10(u,p,t)
    #d_V_DIP*DIP_free
    p[28]*u[5]
end
function affect10!(integrator)
    integrator.u[5] -= 1
    nothing
end
push!(jumps, ConstantRateJump(rate10, affect10!))
push!(jumptovars, [5])

function rate11(u,p,t)
    #k_diss*DIP_bound
    p[3]*u[6]
end
function affect11!(integrator)
    integrator.u[5] += 1
    integrator.u[6] -= 1
    nothing
end
push!(jumps, ConstantRateJump(rate11, affect11!))
push!(jumptovars, [5,6])

#    #dDIP_bound
#    du[6] = k_bind*DIP_free - (k_fuse + k_diss + d_V_DIP)*DIP_bound
function rate12(u,p,t)
    #k_fuse*DIP_bound
    p[4]*u[6]
end
function affect12!(integrator)
    integrator.u[6] -= 1
    integrator.u[7] += 1
    nothing
end
push!(jumps, ConstantRateJump(rate12, affect12!))
push!(jumptovars, [6,7])

function rate13(u,p,t)
    #d_V_DIP*DIP_bound
    p[28]*u[6]
end
function affect13!(integrator)
    integrator.u[6] -= 1
    nothing
end
push!(jumps, ConstantRateJump(rate13, affect13!))
push!(jumptovars, [6])

#    #dDIP_endosome
#    du[7] = k_fuse*DIP_bound - (k_uncoat + d_endosome_DIP)*DIP_endosome
function rate14(u,p,t)
    #k_uncoat*DIP_endosome
    p[5]*u[7]
end
function affect14!(integrator)
    integrator.u[7] -= 1
    integrator.u[8] += 1
    nothing
end
push!(jumps, ConstantRateJump(rate14, affect14!))
push!(jumptovars, [7,8])

function rate15(u,p,t)
    #d_endosome_DIP*DIP_endosome
    p[29]*u[7]
end
function affect15!(integrator)
    integrator.u[7] -= 1
    nothing
end
push!(jumps, ConstantRateJump(rate15, affect15!))
push!(jumptovars, [7])

#    #dDIP_gRNA_plus
#    du[8] = k_uncoat*DIP_endosome - d_gRNA_DIP*DIP_gRNA_plus
function rate16(u,p,t)
    #d_gRNA_DIP*DIP_gRNA_plus
    p[30]*u[8]
end
function affect16!(integrator)
    integrator.u[8] -= 1
    nothing
end
push!(jumps, ConstantRateJump(rate16, affect16!))
push!(jumptovars, [8])

#    du[9] = k_transl*f_ORF1*gRNA_plus - d_NSP*NSP -
#           (k_trans_minus*gRNA_plus + k_trans_plus*gRNA_minus +
#        k_trans_minus_DIP*DIP_gRNA_plus + k_trans_plus_DIP*DIP_gRNA_minus)*NSP
function rate17(u,p,t)
    #k_transl*f_ORF1*gRNA_plus
    p[8]*p[9]*u[4]
end
function affect17!(integrator)
    integrator.u[9] += 1
    nothing
end
push!(jumps, ConstantRateJump(rate17, affect17!))
push!(jumptovars, [9])

function rate18(u,p,t)
    #d_NSP*NSP
    p[10]*u[9]
end
function affect18!(integrator)
    integrator.u[9] -= 1
    nothing
end
push!(jumps, ConstantRateJump(rate18, affect18!))
push!(jumptovars, [9])

function rate19(u,p,t)
    #k_trans_minus*gRNA_plus*NSP
    p[31]*u[4]*u[9]
end
function affect19!(integrator)
    integrator.u[9] -= 1
    nothing
end
push!(jumps, ConstantRateJump(rate19, affect19!))
push!(jumptovars, [9])

function rate20(u,p,t)
    #k_trans_plus*gRNA_minus*NSP
    p[32]*u[10]*u[9]
end
function affect20!(integrator)
    integrator.u[9] -= 1
    nothing
end
push!(jumps, ConstantRateJump(rate20, affect20!))
push!(jumptovars, [9])

function rate21(u,p,t)
    #k_trans_minus_DIP*DIP_gRNA_plus*NSP
    p[33]*u[8]*u[9]
end
function affect21!(integrator)
    integrator.u[9] -= 1
    nothing
end
push!(jumps, ConstantRateJump(rate21, affect21!))
push!(jumptovars, [9])

function rate22(u,p,t)
    #k_trans_plus_DIP*DIP_gRNA_minus*NSP
    p[34]*u[12]*u[9]
end
function affect22!(integrator)
    integrator.u[9] -= 1
    nothing
end
push!(jumps, ConstantRateJump(rate22, affect22!))
push!(jumptovars, [9])

#    f_RdRp = NSP/(NSP + K_NSP)
#    f_complex = N/(N + K_N)
#    f_assemb = SP/(SP + K_V_rel*n_SP)
#    f_assemb_DIP = SP/(SP + K_V_rel_DIP*n_SP_DIP)
#
#    #dgRNA_minus
#    du[10] = k_tr_minus*gRNA_plus*f_RdRp - d_gRNA_minus*gRNA_minus
function rate23(u,p,t)
    #k_tr_minus*gRNA_plus*f_RdRp
    p[11]*u[4]*u[9]/(u[9]+p[12])
end
function affect23!(integrator)
    integrator.u[10] += 1
    nothing
end
push!(jumps, ConstantRateJump(rate23, affect23!))
push!(jumptovars, [10])

function rate24(u,p,t)
    #d_gRNA_minus*gRNA_minus
    p[13]*u[10]
end
function affect24!(integrator)
    integrator.u[10] -= 1
    nothing
end
push!(jumps, ConstantRateJump(rate24, affect24!))
push!(jumptovars, [10])

#    #dgRNA
#    du[11] = k_tr_plus*gRNA_minus*f_RdRp - (k_complex*f_complex + d_gRNA)*gRNA
function rate25(u,p,t)
    #k_tr_plus*gRNA_minus*f_RdRp
    p[14]*u[10]*u[9]/(u[9]+p[12])
end
function affect25!(integrator)
    integrator.u[11] += 1
    nothing
end
push!(jumps, ConstantRateJump(rate25, affect25!))
push!(jumptovars, [11])

function rate26(u,p,t)
    #d_gRNA*gRNA
    p[7]*u[11]
end
function affect26!(integrator)
    integrator.u[11] -= 1
    nothing
end
push!(jumps, ConstantRateJump(rate26, affect26!))
push!(jumptovars, [11])

function rate27(u,p,t)
    #k_complex*f_complex*gRNA
    p[15]*u[14]/(u[14]+p[16])*u[11]
end
function affect27!(integrator)
    integrator.u[11] -= 1
    integrator.u[14] -= integrator.p[22] # -= n_N
    integrator.u[16] += 1
    nothing
end
push!(jumps, ConstantRateJump(rate27, affect27!))
push!(jumptovars, [11,14,16])

#    #dDIP_gRNA_minus
# du[12] = k_tr_minus_DIP*DIP_gRNA_plus*f_RdRp - d_gRNA_minus_DIP*DIP_gRNA_minus
function rate28(u,p,t)
    #k_tr_minus_DIP*DIP_gRNA_plus*f_RdRp
    p[35]*u[8]*u[9]/(u[9]+p[12])
end
function affect28!(integrator)
    integrator.u[12] += 1
    nothing
end
push!(jumps, ConstantRateJump(rate28, affect28!))
push!(jumptovars, [12])

function rate29(u,p,t)
    #d_gRNA_minus_DIP*DIP_gRNA_minus
    p[36]*u[12]
end
function affect29!(integrator)
    integrator.u[12] -= 1
    nothing
end
push!(jumps, ConstantRateJump(rate29, affect29!))
push!(jumptovars, [12])

#    #dDIP_gRNA
#    du[13] = k_tr_plus_DIP*DIP_gRNA_minus*f_RdRp -
#            (k_complex_DIP*f_complex + d_gRNA_DIP)*DIP_gRNA
function rate30(u,p,t)
    #k_tr_plus_DIP*DIP_gRNA_minus*f_RdRp
    p[37]*u[12]*u[9]/(u[9]+p[12])
end
function affect30!(integrator)
    integrator.u[13] += 1
    nothing
end
push!(jumps, ConstantRateJump(rate30, affect30!))
push!(jumptovars, [13])

function rate31(u,p,t)
    #d_gRNA_DIP*DIP_gRNA
    p[30]*u[13]
end
function affect31!(integrator)
    integrator.u[13] -= 1
    nothing
end
push!(jumps, ConstantRateJump(rate31, affect31!))
push!(jumptovars, [13])

function rate32(u,p,t)
    #k_complex_DIP*f_complex*DIP_gRNA
    p[38]*u[14]/(u[14]+p[16])*u[13]
end
function affect32!(integrator)
    integrator.u[13] -= 1
    integrator.u[14] -= integrator.p[39] # -= n_N_DIP
    integrator.u[17] += 1
    nothing
end
push!(jumps, ConstantRateJump(rate32, affect32!))
push!(jumptovars, [13,14,17])

#    #dN
#    du[14] = k_transl*f_N*gRNA - k_complex*n_N*f_complex*gRNA -
#             k_complex_DIP*n_N_DIP*f_complex*DIP_gRNA - d_N*N
function rate33(u,p,t)
    #k_transl*f_N*gRNA
    p[8]*p[17]*u[11]
end
function affect33!(integrator)
    integrator.u[14] += 1
    nothing
end
push!(jumps, ConstantRateJump(rate33, affect33!))
push!(jumptovars, [14])

function rate34(u,p,t)
    #d_N*N
    p[19]*u[14]
end
function affect34!(integrator)
    integrator.u[14] -= 1
    nothing
end
push!(jumps, ConstantRateJump(rate34, affect34!))
push!(jumptovars, [14])

#    #dSP
#    du[15] = k_transl*f_SP*gRNA - k_assembl*n_SP*f_assemb*N_gRNA -
#             k_assembl_DIP*n_SP_DIP*f_assemb_DIP*DIP_N_gRNA - d_SP*SP
function rate35(u,p,t)
    #k_transl*f_SP*gRNA
    p[8]*p[18]*u[11]
end
function affect35!(integrator)
    integrator.u[15] += 1
    nothing
end
push!(jumps, ConstantRateJump(rate35, affect35!))
push!(jumptovars, [15])

function rate36(u,p,t)
    #d_SP*SP
    p[20]*u[15]
end
function affect36!(integrator)
    integrator.u[15] -= 1
    nothing
end
push!(jumps, ConstantRateJump(rate36, affect36!))
push!(jumptovars, [15])

function rate37(u,p,t)
    #k_assembl*f_assemb*N_gRNA
    p[24]*u[15]/(u[15]+p[23]*p[21])*u[16]
end
function affect37!(integrator)
    integrator.u[15] -= integrator.p[21] # -= n_SP
    integrator.u[16] -= 1
    integrator.u[18] += 1
    nothing
end
push!(jumps, ConstantRateJump(rate37, affect37!))
push!(jumptovars, [15,16,18])

function rate38(u,p,t)
    #k_assembl_DIP*f_assemb_DIP*DIP_N_gRNA
    p[40]*u[15]/(u[15]+p[42]*p[41])*u[17]
end
function affect38!(integrator)
    integrator.u[15] -= integrator.p[41] # -= n_SP_DIP
    integrator.u[17] -= 1
    integrator.u[20] += 1
    nothing
end
push!(jumps, ConstantRateJump(rate38, affect38!))
push!(jumptovars, [15,17,20])

#    #dN_gRNA
#    du[16] = k_complex*f_complex*gRNA - (k_assembl*f_assemb + d_N_gRNA)*N_gRNA
function rate39(u,p,t)
    #d_N_gRNA*N_gRNA
    p[25]*u[16]
end
function affect39!(integrator)
    integrator.u[16] -= 1
    nothing
end
push!(jumps, ConstantRateJump(rate39, affect39!))
push!(jumptovars, [16])

#    #dDIP_N_gRNA
#    du[17] = k_complex_DIP*f_complex*DIP_gRNA -
#            (k_assembl_DIP*f_assemb_DIP + d_N_gRNA_DIP)*DIP_N_gRNA
function rate40(u,p,t)
    #d_N_gRNA_DIP*DIP_N_gRNA
    p[43]*u[17]
end
function affect40!(integrator)
    integrator.u[17] -= 1
    nothing
end
push!(jumps, ConstantRateJump(rate40, affect40!))
push!(jumptovars, [17])

#    #dV_assembled
#    du[18] = k_assembl*f_assemb*N_gRNA - (k_release + d_assembled)*V_assembled
#    #dV_released
#    du[19] = k_release*V_assembled - d_V*V_released
function rate41(u,p,t)
    #k_release*V_assembled
    p[26]*u[18]
end
function affect41!(integrator)
    integrator.u[18] -= 1
    integrator.u[19] += 1
    nothing
end
push!(jumps, ConstantRateJump(rate41, affect41!))
push!(jumptovars, [18,19])

function rate42(u,p,t)
    #d_assembled*V_assembled
    p[27]*u[18]
end
function affect42!(integrator)
    integrator.u[18] -= 1
    nothing
end
push!(jumps, ConstantRateJump(rate42, affect42!))
push!(jumptovars, [18])

function rate43(u,p,t)
    #d_V*V_released
    p[2]*u[19]
end
function affect43!(integrator)
    integrator.u[19] -= 1
    nothing
end
push!(jumps, ConstantRateJump(rate43, affect43!))
push!(jumptovars, [19])

#    #dDIP_assembled
#    du[20] = k_assembl_DIP*f_assemb_DIP*DIP_N_gRNA -
#            (k_release_DIP + d_assembled_DIP)*DIP_assembled
#    #dDIP_released
#    du[21] = k_release_DIP*DIP_assembled - d_V_DIP*DIP_released
function rate44(u,p,t)
    #k_release_DIP*DIP_assembled
    p[44]*u[20]
end
function affect44!(integrator)
    integrator.u[20] -= 1
    integrator.u[21] += 1
    nothing
end
push!(jumps, ConstantRateJump(rate44, affect44!))
push!(jumptovars, [20,21])

function rate45(u,p,t)
    #d_assembled_DIP*DIP_assembled
    p[45]*u[20]
end
function affect45!(integrator)
    integrator.u[20] -= 1
    nothing
end
push!(jumps, ConstantRateJump(rate45, affect45!))
push!(jumptovars, [20])

function rate46(u,p,t)
    #d_V_DIP*DIP_released
    p[28]*u[21]
end
function affect46!(integrator)
    integrator.u[21] -= 1
    nothing
end
push!(jumps, ConstantRateJump(rate46, affect46!))
push!(jumptovars, [21])

depgraph = [vcat(vartojumps[vars]...) for vars in jumptovars]

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

discrete_prob = DiscreteProblem(Int64.(u0), tspan, p)
jump_prob = JumpProblem(discrete_prob,
                        RSSACR(), #RSSA(), #SortingDirect(), #Direct(),
                        jumps...,
                        #dep_graph=depgraph,
                        vartojumps_map=vartojumps,
                        jumptovars_map=jumptovars,
                        save_positions=(false,false))

@time sol = solve(jump_prob, SSAStepper(), saveat=1.0) # precompile

