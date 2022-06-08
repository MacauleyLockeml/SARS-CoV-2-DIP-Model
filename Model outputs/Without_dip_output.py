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
    V[4]=0#-k_bind*X[4]-d_VD*X[4]+k_diss*X[5]
        
        #bound DIP X[5]
    V[5]=0#k_bind*X[4]-(k_fuse+k_diss+d_VD)*X[5]
        
        #DIP in an endosome X[6]
    V[6]=0#k_fuse*X[5]-(k_uncoat+d_endosomeD)*X[6]
        
        #DIP gRNA(+) X[7]
    V[7]=0#k_uncoat*X[6]-d_gRNAD*X[7]
    
        #NSP X[8]
    V[8]=k_transl*f_orf*X[3]-d_NSP*X[8]-(k_transWT_neg*V[3]+k_transWT_pos*X[9]+k_transD_neg*X[7]+k_transD_pos*X[11])*X[8]
    
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
    

if __name__ == '__main__':
    t = np.linspace(0.0, 48, 4900) #Timecourse for integration

    samples = df=pd.read_csv('Param_values_original.csv', sep=',',header=None,index_col=(0))
    bounds=df.values.reshape(1,41)[0]
    k_bind,d_VW,k_diss,k_fuse,k_uncoat,d_endosomeW,d_gRNAW,d_VD,d_endosomeD,d_gRNAD,k_transl,d_NSP,k_trWT_neg,K_NSP,d_gRNAWT_neg,k_trWT_pos,k_WTcomplex,k_trD_neg,d_gRNAD_neg,k_trD_pos,k_Dcomplex,K_N,d_N,d_SP,n_SPWT,	n_SPD,	n_ND,K_VWT_rel,k_assemblWT,K_VD_rel,k_assemblD,d_NgRNAWT,k_releaseWT,d_assembledWT,d_N_gRNAD,k_releaseD,d_assembledD,k_transWT_neg,k_transWT_pos,k_transD_neg,k_transD_pos= bounds
    f_orf=1/21000
    f_N=1/1200
    f_SP=1/10000
    n_NW=456
    
    X0=np.zeros(21)
    X0[0]=10
    
    X, infodict = integrate.odeint(eqs, X0, t, args=(f_orf,f_N,f_SP,n_NW,k_bind,d_VW,k_diss,k_fuse,k_uncoat,d_endosomeW,d_gRNAW,d_VD,d_endosomeD,d_gRNAD,k_transl,d_NSP,k_trWT_neg,K_NSP,d_gRNAWT_neg,k_trWT_pos,k_WTcomplex,k_trD_neg,d_gRNAD_neg,k_trD_pos,k_Dcomplex,K_N,d_N,d_SP,n_SPWT,	n_SPD,	n_ND,K_VWT_rel,k_assemblWT,K_VD_rel,k_assemblD,d_NgRNAWT,k_releaseWT,d_assembledWT,d_N_gRNAD,k_releaseD,d_assembledD,k_transWT_neg,k_transWT_pos,k_transD_neg,k_transD_pos), full_output = 1)
    
    
    _labels=['$[V_{free}^{WT}]$','$[V_{bound}^{WT}]$','$[V_{endosome}^{WT}]$','$[gRNA_{(+)}^{WT}]$','$[V_{free}^{dip}]$','$[V_{bound}^{dip}]$','$[V_{endosome}^{dip}]$','$[gRNA_{(+)}^{dip}]$','$[NSP]$','$[gRNA_{(-)}^{WT}]$','$[gRNA^{WT}]$','$[gRNA_{(-)}^{dip}]$','$[gRNA^{dip}]$','$[N]$','$[SP]$','$[N-gRNA^{WT}]$','$[V_{assembled}^{WT}]$','$[V_{released}^{WT}]$','$[N-gRNA^{dip}]$','$[V_{assembled}^{dip}]$','$[V_{released}^{dip}]$']
    _WT_index=[0,1,2,3,8,9,10,13,14,15,16,17]
    _WT_low=[0,1,2,3]
    _DIP_index=[4,5,6,7,11,12,18,19,20]
    

    figl,axl=plt.subplots(figsize=(6,6),dpi=300)
    axl.plot(t, X[:,0], label=_labels[0])
    axl.plot(t, X[:,1], label=_labels[1])
    #axl.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
    axl.set(title='Binding and fusion')
    axl.grid('on')
    plt.xlim(xmin=0,xmax=5)
    plt.legend(loc='upper right', fancybox=True, shadow=True, title='Parameter',prop={'size': 11})
    figl.savefig("Binding and fusion.png")
    
    figl,axl=plt.subplots(figsize=(6,6),dpi=300)
    axl.plot(t, X[:,2], label=_labels[2])
    axl.plot(t, X[:,3], label=_labels[3])
    #axl.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
    axl.set(title='Endocytosis and uncoating')
    axl.grid('on')
    plt.xlim(xmin=0,xmax=25)
    plt.legend(loc='upper right', fancybox=True, shadow=True, title='Parameter',prop={'size': 11})
    figl.savefig("Endocytosis and uncoating.png")
    
    figl,axl=plt.subplots(figsize=(6,6),dpi=300)
    axl.plot(t, X[:,8], label=_labels[8])
    axl.plot(t, X[:,9], label=_labels[9])
    #axl.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
    axl.set(title='ORF1 translation and $gRNA_(-)$ Synthesis')
    axl.grid('on')
    plt.xlim(xmin=0,xmax=25)
    plt.ylim(ymin=0, ymax=50)
    plt.legend(loc='upper right', fancybox=True, shadow=True, title='Parameter',prop={'size': 11})
    figl.savefig("ORF1 translation and gRNA_(-) Synthesis.png")
     
    
    figl,axl=plt.subplots(figsize=(6,6),dpi=300)
    axl.plot(t, X[:,10], label=_labels[10])
    axl.plot(t, X[:,13], label=_labels[13])
    #axl.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
    axl.set(title='Transcription and translation')
    axl.grid('on')
    plt.yscale('symlog')
    plt.legend(loc='lower right', fancybox=True, shadow=True, title='Parameter',prop={'size': 11})
    figl.savefig("Transcription and translation.png")
    
    figl,axl=plt.subplots(figsize=(6,6),dpi=300)
    axl.plot(t, X[:,14], label=_labels[14])
    axl.plot(t, X[:,15], label=_labels[15])
    #axl.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
    axl.set(title='Translation and Nucleocaspid formation')
    axl.grid('on')
    plt.yscale('symlog')
    plt.legend(loc='lower right', fancybox=True, shadow=True, title='Parameter',prop={'size': 11})
    figl.savefig("Translation and Nucleocaspid formation.png")
    
    figl,axl=plt.subplots(figsize=(6,6),dpi=300)
    axl.plot(t, X[:,16], label=_labels[16])
    axl.plot(t, X[:,17], label=_labels[17])
    #axl.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
    axl.set(title='Assembly and release')
    axl.grid('on')
    plt.yscale('symlog')
    plt.legend(loc='right', fancybox=True, shadow=True, title='Parameter',prop={'size': 11})
    figl.savefig("Assembly and release.png")

            