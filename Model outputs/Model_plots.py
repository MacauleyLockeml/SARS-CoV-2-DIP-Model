 # -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 11:18:19 2022

@author: Macauley Locke
@contact: mmmwl@leeds.ac.uk,
@date: 18/03/2022

The following code is used to generate model output figures for a fixed set 
"""


import numpy as np
from scipy import integrate
import pandas as pd
import csv
import matplotlib.pyplot as plt
#These lines are for running the code simultaneously on ARC

def eqs(X,t,f_orf,f_N,f_SP,n_NW,k_bind,d_VW,k_diss,k_fuse,k_uncoat,d_endosomeW,d_gRNAW,d_VD,d_endosomeD,d_gRNAD,k_transl,d_NSP,k_trWT_neg,K_NSP,d_gRNAWT_neg,k_trWT_pos,k_WTcomplex,k_trD_neg,d_gRNAD_neg,k_trD_pos,k_Dcomplex,K_N,d_N,d_SP,n_SPWT,	n_SPD,	n_ND,K_VWT_rel,k_assemblWT,K_VD_rel,k_assemblD,d_NgRNAWT,k_releaseWT,d_assembledWT,d_N_gRNAD,k_releaseD,d_assembledD,k_transWT_neg,k_transWT_pos,k_transD_neg,k_transD_pos):
        
     #Start empty array to save for output
    V=[0]*21
        
        #free WT virus X[0]
    V[0]=-k_bind*X[0]-d_VW*X[0]+k_diss*X[1]
        
        #bound WT virus X[1]
    V[1]=k_bind*X[0]-(k_fuse+k_diss+d_VW)*X[1]
        
        #WT virus an endosome X[2]
    V[2]=k_fuse*X[1]-(k_uncoat+d_endosomeW)*X[2]
        
        #WT gRNA(+) X[3]
    V[3]=k_uncoat*X[2]-d_gRNAW*X[3]
        
        #free DIP X[4]
    V[4]=-k_bind*X[4]-d_VD*X[4]+k_diss*X[5]
        
        #bound DIP X[5]
    V[5]=k_bind*X[4]-(k_fuse+k_diss+d_VD)*X[5]
        
        #DIP in an endosome X[6]
    V[6]=k_fuse*X[5]-(k_uncoat+d_endosomeD)*X[6]
        
        #DIP gRNA(+) X[7]
    V[7]=k_uncoat*X[6]-d_gRNAD*X[7]
    
        #NSP X[8]
    V[8]=k_transl*f_orf*X[3]-d_NSP*X[8]-(k_transWT_neg*X[3]+k_transWT_pos*X[9]+k_transD_neg*X[7]+k_transD_pos*X[11])*X[8]
    
        #WT gRNA(-) X[9]
    V[9]=(k_trWT_neg*X[3])*(X[8]/(X[8]+K_NSP))-d_gRNAWT_neg*X[9]
        
        #WT gRNA X[10]
    V[10]=k_trWT_pos*X[9]*(X[8]/(X[8]+K_NSP))-(k_WTcomplex*(X[13]/(X[13]+K_N))+d_gRNAW)*X[10]
        
        #DIP gRNA(-) X[11]
    V[11]=k_trD_neg*X[7]*(X[8]/(X[8]+K_NSP))-d_gRNAD_neg*X[11]
        
        #DIP gRNA X[12]
    V[12]=k_trD_pos*X[11]*(X[8]/(X[8]+K_NSP))-(k_Dcomplex*(X[13]/(X[13]+K_N))+d_gRNAD)*X[12]
        
        #N protein X[13]
    V[13]=k_transl*f_N*X[10]-k_WTcomplex*n_NW*(X[13]/(X[13]+K_N))*X[10]-k_Dcomplex*n_ND*(X[13]/(X[13]+K_N))*X[12]-d_N*X[13]
        
        #SP proteins X[14]
    V[14]=k_transl*f_SP*X[10]-k_assemblWT*n_SPWT*(X[14]/(X[14]+K_VWT_rel*n_SPWT))*X[15]-d_SP*X[14]-k_assemblD*n_SPD*(X[14]/(X[14]+K_VD_rel*n_SPD))*X[18]
        
        #WT N-gRNA X[15]
    V[15]=k_WTcomplex*(X[13]/(X[13]+K_N))*X[10]-(k_assemblWT*(X[14]/(X[14]+K_VWT_rel*n_SPWT))+d_NgRNAWT)*X[15]
        
        #WT assembled X[16]
    V[16]=k_assemblWT*(X[14]/(X[14]+K_VWT_rel*n_SPWT))*X[15]-(k_releaseWT+d_assembledWT)*X[16]
        
        #WT released X[17]
    V[17]=k_releaseWT*X[16]-d_VW*X[17]
        
        #DIP n-gRNA X[18]
    V[18]=k_Dcomplex*(X[13]/(X[13]+K_N))*X[12]-(k_assemblD*(X[14]/(X[14]+K_VD_rel*n_SPD))+d_N_gRNAD)*X[18]
        
        #DIP assembly X[19]
    V[19]=k_assemblD*(X[14]/(X[14]+K_VD_rel*n_SPD))*X[18]-(k_releaseD+d_assembledD)*X[19]
        
        #DIP released
    V[20]=k_releaseD*X[19]-d_VD*X[20]
        
    return V
    

if __name__ == '__main__':
    t = np.linspace(0.0, 48, 4900) #Timecourse for integration

    samples = df=pd.read_csv('Param_values_original.csv', sep=',',header=None,index_col=(0))
    bounds=df.values.reshape(1,41)[0]
    k_bind=12
    d_VW=0.12
    k_diss=0.61
    k_fuse=0.5
    k_uncoat=0.5
    d_endosomeW=0.06
    d_gRNAW=0.2
    d_VD=0.481
    d_endosomeD=0.00329
    d_gRNAD=0.218
    k_transl=45360
    d_NSP=0.069
    k_trWT_neg=3
    K_NSP=100
    d_gRNAWT_neg=0.1
    k_trWT_pos=1000
    k_WTcomplex=0.4
    k_trD_neg=34
    d_gRNAD_neg=0.218
    k_trD_pos=2540
    k_Dcomplex=0.140
    K_N=5*10**6
    d_N=0.023
    d_SP=0.044
    n_SPWT=2000
    n_SPD=112
    n_ND=53
    K_VWT_rel=1000
    k_assemblWT=1
    K_VD_rel=380
    k_assemblD=0.38
    d_NgRNAWT=0.2
    k_releaseWT=8
    d_assembledWT=0.06
    d_N_gRNAD=0.268
    k_releaseD=105
    d_assembledD=0.00489
    k_transWT_neg=0.0000539
    k_transWT_pos=0.00617
    k_transD_neg=0.0000472
    k_transD_pos= 0.00861
    f_orf=1/21000
    f_N=1/1200
    f_SP=1/10000
    n_NW=456
    
    X0=np.zeros(21)
    X0[0]=10
    X0[4]=10
    
    X, infodict = integrate.odeint(eqs, X0, t, args=(f_orf,f_N,f_SP,n_NW,k_bind,d_VW,k_diss,k_fuse,k_uncoat,d_endosomeW,d_gRNAW,d_VD,d_endosomeD,d_gRNAD,k_transl,d_NSP,k_trWT_neg,K_NSP,d_gRNAWT_neg,k_trWT_pos,k_WTcomplex,k_trD_neg,d_gRNAD_neg,k_trD_pos,k_Dcomplex,K_N,d_N,d_SP,n_SPWT,	n_SPD,	n_ND,K_VWT_rel,k_assemblWT,K_VD_rel,k_assemblD,d_NgRNAWT,k_releaseWT,d_assembledWT,d_N_gRNAD,k_releaseD,d_assembledD,k_transWT_neg,k_transWT_pos,k_transD_neg,k_transD_pos), full_output = 1)
    
    
    _labels=['$[V_{free}^{wt}]$','$[V_{bound}^{wt}]$','$[V_{endosome}^{wt}]$','$[gRNA_{(+)}^{wt}]$','$[V_{free}^{dip}]$','$[V_{bound}^{dip}]$','$[V_{endosome}^{dip}]$','$[gRNA_{(+)}^{dip}]$','$[NSP]$','$[gRNA_{(-)}^{wt}]$','$[gRNA^{wt}]$','$[gRNA_{(-)}^{dip}]$','$[gRNA^{dip}]$','$[N]$','$[SP]$','$[N-gRNA^{wt}]$','$[V_{assembled}^{wt}]$','$[V_{released}^{wt}]$','$[N-gRNA^{dip}]$','$[V_{assembled}^{dip}]$','$[V_{released}^{dip}]$']
    _WT_index=[0,1,2,3,8,9,10,13,14,15,16,17]
    _WT_low=[0,1,2,3]
    _DIP_index=[4,5,6,7,11,12,18,19,20]
    
    def eqs_ori(X,t,f_orf,f_N,f_SP,n_NW,k_bind,d_VW,k_diss,k_fuse,k_uncoat,d_endosomeW,d_gRNAW,d_VD,d_endosomeD,d_gRNAD,k_transl,d_NSP,k_trWT_neg,K_NSP,d_gRNAWT_neg,k_trWT_pos,k_WTcomplex,k_trD_neg,d_gRNAD_neg,k_trD_pos,k_Dcomplex,K_N,d_N,d_SP,n_SPWT,	n_SPD,	n_ND,K_VWT_rel,k_assemblWT,K_VD_rel,k_assemblD,d_NgRNAWT,k_releaseWT,d_assembledWT,d_N_gRNAD,k_releaseD,d_assembledD,k_transWT_neg,k_transWT_pos,k_transD_neg,k_transD_pos):
        
            #Start empty array to save for output
        V=[0]*21
            
            #free WT virus X[0]
        V[0]=-k_bind*X[0]-d_VW*X[0]+k_diss*X[1]
            
            #bound WT virus X[1]
        V[1]=k_bind*X[0]-(k_fuse+k_diss+d_VW)*X[1]
            
            #WT virus an endosome X[2]
        V[2]=k_fuse*X[1]-(k_uncoat+d_endosomeW)*X[2]
            
            #WT gRNA(+) X[3]
        V[3]=k_uncoat*X[2]-d_gRNAW*X[3]
            
            #free DIP X[4]
        V[4]=0#-k_bind*X[4]-d_VD*X[4]+k_diss*X[5]
            
            #bound DIP X[5]
        V[5]=0#k_bind*X[4]-(k_fuse+k_diss+d_VD)*X[5]
            
            #DIP in an endosome X[6]
        V[6]=0#k_fuse*X[5]-(k_uncoat+d_endosomeD)*X[6]
            
            #DIP gRNA(+) X[7]
        V[7]=0#k_uncoat*X[6]-d_gRNAD*X[7]
        
            #NSP X[8]
        V[8]=k_transl*f_orf*X[3]-d_NSP*X[8]#-(k_transWT_neg*X[3]+k_transWT_pos*X[9]+k_transD_neg*X[7]+k_transD_pos*X[11])*X[8]
        
            #WT gRNA(-) X[9]
        V[9]=(k_trWT_neg*X[3])*(X[8]/(X[8]+K_NSP))-d_gRNAWT_neg*X[9]
            
            #WT gRNA X[10]
        V[10]=k_trWT_pos*X[9]*(X[8]/(X[8]+K_NSP))-(k_WTcomplex*(X[13]/(X[13]+K_N))+d_gRNAW)*X[10]
            
            #DIP gRNA(-) X[11]
        V[11]=0#k_trD_neg*X[7]*(X[8]/(X[8]+K_NSP))-d_gRNAD_neg*X[11]
            
            #DIP gRNA X[12]
        V[12]=0#k_trD_pos*X[11]*(X[8]/(X[8]+K_NSP))-(k_Dcomplex*(X[13]/(X[13]+K_N))+d_gRNAD)*X[12]
            
            #N protein X[13]
        V[13]=k_transl*f_N*X[10]-k_WTcomplex*n_NW*(X[13]/(X[13]+K_N))*X[10]-k_Dcomplex*n_ND*(X[13]/(X[13]+K_N))*X[12]-d_N*X[13]
            
            #SP proteins X[14]
        V[14]=k_transl*f_SP*X[10]-k_assemblWT*n_SPWT*(X[14]/(X[14]+K_VWT_rel*n_SPWT))*X[15]-d_SP*X[14]#-k_assemblD*n_SPD*(X[14]/(X[14]+K_VD_rel*n_SPD))*X[18]
            
            #WT N-gRNA X[15]
        V[15]=k_WTcomplex*(X[13]/(X[13]+K_N))*X[10]-(k_assemblWT*(X[14]/(X[14]+K_VWT_rel*n_SPWT))+d_NgRNAWT)*X[15]
            
            #WT assembled X[16]
        V[16]=k_assemblWT*(X[14]/(X[14]+K_VWT_rel*n_SPWT))*X[15]-(k_releaseWT+d_assembledWT)*X[16]
            
            #WT released X[17]
        V[17]=k_releaseWT*X[16]-d_VW*X[17]
            
            #DIP n-gRNA X[18]
        V[18]=0#k_Dcomplex*(X[13]/(X[13]+K_N))*X[12]-(k_assemblD*(X[14]/(X[14]+K_VD_rel*n_SPD))+d_N_gRNAD)*X[18]
            
            #DIP assembly X[19]
        V[19]=0#k_assemblD*(X[14]/(X[14]+K_VD_rel*n_SPD))*X[18]-(k_releaseD+d_assembledD)*X[19]
            
            #DIP released
        V[20]=0#k_releaseD*X[19]-d_VW*X[20]
            
        return V

    k_bind=12
    d_VW=0.12
    k_diss=0.61
    k_fuse=0.5
    k_uncoat=0.5
    d_endosomeW=0.06
    d_gRNAW=0.2
    d_VD=0
    d_endosomeD=0
    d_gRNAD=0
    k_transl=45360
    d_NSP=0.069
    k_trWT_neg=3
    K_NSP=100
    d_gRNAWT_neg=0.1
    k_trWT_pos=1000
    k_WTcomplex=0.4
    k_trD_neg=0
    d_gRNAD_neg=0
    k_trD_pos=0
    k_Dcomplex=0
    K_N=5*10**6
    d_N=0.023
    d_SP=0.044
    n_SPWT=2000
    n_SPD=0
    n_ND=0
    K_VWT_rel=1000
    k_assemblWT=1
    K_VD_rel=0
    k_assemblD=0
    d_NgRNAWT=0.2
    k_releaseWT=8
    d_assembledWT=0.06
    d_N_gRNAD=0
    k_releaseD=0
    d_assembledD=0
    k_transWT_neg=0
    k_transWT_pos=0
    k_transD_neg=0
    k_transD_pos= 0
    f_orf=1/21000
    f_N=1/1200
    f_SP=1/10000
    n_NW=456
    
    X0=np.zeros(21)
    X0[0]=10
    
    t2 = np.linspace(0.0, 48, 4900)

        
    Orig, infodict = integrate.odeint(eqs_ori, X0, t2, args=(f_orf,f_N,f_SP,n_NW,k_bind,d_VW,k_diss,k_fuse,k_uncoat,d_endosomeW,d_gRNAW,d_VD,d_endosomeD,d_gRNAD,k_transl,d_NSP,k_trWT_neg,K_NSP,d_gRNAWT_neg,k_trWT_pos,k_WTcomplex,k_trD_neg,d_gRNAD_neg,k_trD_pos,k_Dcomplex,K_N,d_N,d_SP,n_SPWT,	n_SPD,	n_ND,K_VWT_rel,k_assemblWT,K_VD_rel,k_assemblD,d_NgRNAWT,k_releaseWT,d_assembledWT,d_N_gRNAD,k_releaseD,d_assembledD,k_transWT_neg,k_transWT_pos,k_transD_neg,k_transD_pos), full_output = 1)
    _labels2=['$[V_{free}]$','$[V_{bound}]$','$[V_{endosome}]$','$[gRNA_{(+)}]$','$[NSP]$','$[gRNA_{(-)}]$','$[gRNA]$','$[N]$','$[SP]$','$[N-gRNA]$','$[V_{assembled}]$','$[V_{released}]$']
        

    figl,axl=plt.subplots(figsize=(6,6),dpi=300)
    axl.plot(t, X[:,0], label=_labels[0])
    axl.plot(t, X[:,1], label=_labels[1])
    axl.plot(t, X[:,4], label=_labels[4])
    axl.plot(t, X[:,5], label=_labels[5])
    axl.plot(t, Orig[:,0], label=_labels2[0],linestyle="dashed")
    axl.plot(t, Orig[:,1], label=_labels2[1],linestyle="dashed")
    #axl.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
    axl.set(title='Binding and fusion')
    axl.grid('on')
    plt.xlim(xmin=0,xmax=10)
    plt.xlabel('t, hours')
    plt.ylabel('particles')
    plt.legend(ncol=3,loc='upper center', bbox_to_anchor=(0.5, -0.10),fancybox=True, shadow=True,prop={'size': 14})
    figl.savefig("Binding and fusion.png",bbox_inches="tight")
    
    fig2,axl=plt.subplots(figsize=(6,6),dpi=300)
    axl.plot(t, X[:,2], label=_labels[2])
    axl.plot(t, X[:,3], label=_labels[3])
    axl.plot(t, X[:,6], label=_labels[6])
    axl.plot(t, X[:,7], label=_labels[7])
    axl.plot(t, Orig[:,2], label=_labels2[2],linestyle="dashed")
    axl.plot(t, Orig[:,3], label=_labels2[3],linestyle="dashed")
    #axl.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
    axl.set(title='Endocytosis and uncoating')
    axl.grid('on')
    plt.xlim(xmin=0,xmax=25)
    plt.xlabel('t, hours')
    plt.ylabel('numbers')
    plt.legend(ncol=3,loc='upper center', bbox_to_anchor=(0.5, -0.10),fancybox=True, shadow=True,prop={'size': 14})
    fig2.savefig("Endocytosis and uncoating.png",bbox_inches="tight")
    
    fig3,axl=plt.subplots(figsize=(6,6),dpi=300)
    axl.plot(t, X[:,8], label=_labels[8])
    axl.plot(t, X[:,9], label=_labels[9])
    axl.plot(t, X[:,11], label=_labels[11])
    axl.plot(t, Orig[:,4], label=_labels2[4],linestyle="dashed")
    axl.plot(t, Orig[:,5], label=_labels2[5],linestyle="dashed")
    #axl.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
    axl.set(title='ORF1 translation and $gRNA_(-)$ Synthesis')
    axl.grid('on')
    plt.xlim(xmin=0,xmax=25)
    plt.ylim(ymin=0, ymax=50)
    plt.xlabel('t, hours')
    plt.ylabel('molecules')
    plt.legend(ncol=3,loc='upper center', bbox_to_anchor=(0.5, -0.10),fancybox=True, shadow=True,prop={'size': 14})
    fig3.savefig("ORF1 translation and gRNA_(-) Synthesis.png",bbox_inches="tight")
     
    
    fig4,axl=plt.subplots(figsize=(6,6),dpi=300)
    axl.plot(t, X[:,13], label=_labels[13])
    axl.plot(t, X[:,10], label=_labels[10])
    axl.plot(t, X[:,12], label=_labels[12])
    axl.plot(t, Orig[:,6], label=_labels2[6],linestyle="dashed")
    axl.plot(t, Orig[:,7], label=_labels2[7],linestyle="dashed")
    #axl.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
    axl.set(title='Transcription and translation')
    axl.grid('on')
    plt.yscale('symlog')
    plt.xlabel('t, hours')
    plt.ylabel('molecules')
    plt.legend(ncol=3,loc='upper center', bbox_to_anchor=(0.5, -0.10),fancybox=True, shadow=True,prop={'size': 14})
    fig4.savefig("Transcription and translation.png",bbox_inches="tight")
    
    fig5,axl=plt.subplots(figsize=(6,6),dpi=300)
    axl.plot(t, X[:,14], label=_labels[14])
    axl.plot(t, X[:,15], label=_labels[15])
    axl.plot(t, X[:,18], label=_labels[18])
    axl.plot(t, Orig[:,8], label=_labels2[8],linestyle="dashed")
    axl.plot(t, Orig[:,9], label=_labels2[9],linestyle="dashed")
    #axl.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
    axl.set(title='Translation and Nucleocaspid formation')
    axl.grid('on')
    plt.yscale('symlog')
    plt.xlabel('t, hours')
    plt.ylabel('molecules')
    plt.legend(ncol=3,loc='upper center', bbox_to_anchor=(0.5, -0.10),fancybox=True, shadow=True,prop={'size': 14})
    fig5.savefig("Translation and Nucleocaspid formation.png",bbox_inches="tight")
    
    fig6,axl=plt.subplots(figsize=(6,6),dpi=300)
    axl.plot(t, X[:,16], label=_labels[16])
    axl.plot(t, X[:,17], label=_labels[17])
    axl.plot(t, X[:,19], label=_labels[19])
    axl.plot(t, X[:,20], label=_labels[20])
    axl.plot(t, Orig[:,10], label=_labels2[10],linestyle="dashed")
    axl.plot(t, Orig[:,11], label=_labels2[11],linestyle="dashed")
    #axl.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
    axl.set(title='Assembly and release')
    axl.grid('on')
    plt.yscale('symlog')
    plt.xlabel('t, hours')
    plt.ylabel('particles')
    plt.legend(ncol=3,loc='upper center', bbox_to_anchor=(0.5, -0.10),fancybox=True, shadow=True,prop={'size': 14})
    fig6.savefig("Assembly and release.png",bbox_inches="tight")
    
    t2 = np.linspace(4, 48, 45)
    
    fig7,ax2=plt.subplots(figsize=(6,6),dpi=300)
    ax2.plot(t2, (X[4:49,17]-X[4:49,20])/X[4:49,20], label=_labels[17])
    #axl.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
    ax2.set(title='Assembly and release')
    ax2.grid('on')
    plt.yscale('symlog')
    plt.xlabel('t, hours')
    plt.legend(loc='right', fancybox=True, shadow=True, title='Parameter',prop={'size': 11})
    fig7.savefig("Assembly and release ratio.png",bbox_inches="tight")
    
    
    
    
    
    
    fig8,axs=plt.subplots(3,2,figsize=(10,10),dpi=300)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=1.24, wspace=0.6, hspace=0.5)
    plt.subplot(3,2,1)
    plt.plot(t, X[:,0], label=_labels[0])
    plt.plot(t, X[:,1], label=_labels[1])
    plt.plot(t, X[:,4], label=_labels[4])
    plt.plot(t, X[:,5], label=_labels[5])
    plt.plot(t, Orig[:,0], label=_labels2[0],linestyle="dashed")
    plt.plot(t, Orig[:,1], label=_labels2[1],linestyle="dashed")
    #axl.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
    axs[0,0].set(title='Binding and fusion')
    #axs[0,0].grid('on')
    plt.xlim(xmin=0,xmax=10)
    plt.xlabel('t, hours')
    plt.ylabel('particles')
    plt.legend(ncol=3,loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True)

    
    plt.subplot(3,2,2)
    plt.plot(t, X[:,2], label=_labels[2])
    plt.plot(t, X[:,3], label=_labels[3])
    plt.plot(t, X[:,6], label=_labels[6])
    plt.plot(t, X[:,7], label=_labels[7])
    plt.plot(t, Orig[:,2], label=_labels2[2],linestyle="dashed")
    plt.plot(t, Orig[:,3], label=_labels2[3],linestyle="dashed")
    #axl.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
    axs[0,1].set(title='Endocytosis and uncoating')
    #axs[0,1].grid('on')
    plt.xlim(xmin=0,xmax=25)
    plt.xlabel('t, hours')
    plt.ylabel('numbers')
    plt.legend(ncol=3,loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True)
    
    plt.subplot(3,2,3)
    plt.plot(t, X[:,8], label=_labels[8])
    plt.plot(t, X[:,9], label=_labels[9])
    plt.plot(t, X[:,11], label=_labels[11])
    plt.plot(t, Orig[:,8], label=_labels2[4],linestyle="dashed")
    plt.plot(t, Orig[:,9], label=_labels2[5],linestyle="dashed")
    #axl.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
    axs[1,0].set(title='ORF1 translation and $gRNA_{(-)}$ Synthesis')
    #axs[0,2].grid('on')
    plt.xlim(xmin=0,xmax=25)
    plt.ylim(ymin=0, ymax=50)
    plt.xlabel('t, hours')
    plt.ylabel('molecules')
    plt.legend(ncol=3,loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True)

     
    
    plt.subplot(3,2,4)
    
    plt.plot(t, X[:,10], label=_labels[10])
    plt.plot(t, X[:,12], label=_labels[12])
    plt.plot(t, X[:,13], label=_labels[13])
    plt.plot(t, Orig[:,10], label=_labels2[6],linestyle="dashed")
    plt.plot(t, Orig[:,13], label=_labels2[7],linestyle="dashed")
    #axl.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
    axs[1,1].set(title='Transcription and translation')
    #axs[1,0].grid('on')
    plt.yscale('symlog')
    plt.xlabel('t, hours')
    plt.ylabel('molecules')
    plt.legend(ncol=3,loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True)
    
    
    plt.subplot(3,2,5)
    plt.plot(t, X[:,14], label=_labels[14])
    plt.plot(t, X[:,15], label=_labels[15])
    plt.plot(t, X[:,18], label=_labels[18])
    plt.plot(t, Orig[:,14], label=_labels2[8],linestyle="dashed")
    plt.plot(t, Orig[:,15], label=_labels2[9],linestyle="dashed")
    #axl.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
    axs[2,0].set(title='Translation and Nucleocaspid formation')
    #axs[1,1].grid('on')
    plt.yscale('symlog')
    plt.xlabel('t, hours')
    plt.ylabel('molecules')
    plt.legend(ncol=3,loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True)
    
    plt.subplot(3,2,6)
    plt.plot(t, X[:,16], label=_labels[16])
    plt.plot(t, X[:,17], label=_labels[17])
    plt.plot(t, X[:,19], label=_labels[19])
    plt.plot(t, X[:,20], label=_labels[20])
    plt.plot(t, Orig[:,16], label=_labels2[10],linestyle="dashed")
    plt.plot(t, Orig[:,17], label=_labels2[11],linestyle="dashed")
    #axl.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
    axs[2,1].set(title='Assembly and release')
   # axs[1,2].grid('on')
    plt.yscale('symlog')
    plt.xlabel('t, hours')
    plt.ylabel('particles')
    plt.legend(ncol=3,loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True)
    fig8.savefig("outputs.png",bbox_inches="tight")

            