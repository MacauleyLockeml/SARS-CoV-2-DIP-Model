# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 10:59:03 2022

@author: macau
"""

import os, glob
import pandas as pd
import numpy as np
from scipy import integrate
import csv
import matplotlib.pyplot as plt

origpath=os.path.dirname(__file__)

os.chdir(os.path.dirname(__file__) + r"\TEST2")

extension = 'csv'
all_filenames = [i for i in glob.glob('*.{}'.format(extension))]



combined_csv = pd.concat([pd.read_csv(f,header=None,index_col=None) for f in all_filenames ])

Bayes=combined_csv.to_numpy()
_samples=10**6#9*10**6
_accepted=0.001
_percent=int(_accepted*_samples)

sorttemp=Bayes[Bayes[:,45].argsort()]    
sort=sorttemp[0:_percent,]

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


    #samples = df=pd.read_csv('Param_values_original.csv', sep=',',header=None,index_col=(0))
    #bounds=df.values.reshape(1,41)[0]
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

t2 = np.linspace(0.0, 48, 490)
t= np.linspace(0.0, 48, 49) 
    
Orig, infodict = integrate.odeint(eqs_ori, X0, t2, args=(f_orf,f_N,f_SP,n_NW,k_bind,d_VW,k_diss,k_fuse,k_uncoat,d_endosomeW,d_gRNAW,d_VD,d_endosomeD,d_gRNAD,k_transl,d_NSP,k_trWT_neg,K_NSP,d_gRNAWT_neg,k_trWT_pos,k_WTcomplex,k_trD_neg,d_gRNAD_neg,k_trD_pos,k_Dcomplex,K_N,d_N,d_SP,n_SPWT,	n_SPD,	n_ND,K_VWT_rel,k_assemblWT,K_VD_rel,k_assemblD,d_NgRNAWT,k_releaseWT,d_assembledWT,d_N_gRNAD,k_releaseD,d_assembledD,k_transWT_neg,k_transWT_pos,k_transD_neg,k_transD_pos), full_output = 1)



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

median=[]

for i in range(46):
    median.append(np.median(sort[:,i]))

f_orf,f_N,f_SP,n_NW,k_bind,d_VW,k_diss,k_fuse,k_uncoat,d_endosomeW,d_gRNAW,d_VD,d_endosomeD,d_gRNAD,k_transl,d_NSP,k_trWT_neg,K_NSP,d_gRNAWT_neg,k_trWT_pos,k_WTcomplex,k_trD_neg,d_gRNAD_neg,k_trD_pos,k_Dcomplex,K_N,d_N,d_SP,n_SPWT,	n_SPD,	n_ND,K_VWT_rel,k_assemblWT,K_VD_rel,k_assemblD,d_NgRNAWT,k_releaseWT,d_assembledWT,d_N_gRNAD,k_releaseD,d_assembledD,k_transWT_neg,k_transWT_pos,k_transD_neg,k_transD_pos,distance=median






X0=np.zeros(21)
X0[0]=10
X0[4]=10

NEW, infodict = integrate.odeint(eqs, X0, t2, args=(f_orf,f_N,f_SP,n_NW,k_bind,d_VW,k_diss,k_fuse,k_uncoat,d_endosomeW,d_gRNAW,d_VD,d_endosomeD,d_gRNAD,k_transl,d_NSP,k_trWT_neg,K_NSP,d_gRNAWT_neg,k_trWT_pos,k_WTcomplex,k_trD_neg,d_gRNAD_neg,k_trD_pos,k_Dcomplex,K_N,d_N,d_SP,n_SPWT,	n_SPD,	n_ND,K_VWT_rel,k_assemblWT,K_VD_rel,k_assemblD,d_NgRNAWT,k_releaseWT,d_assembledWT,d_N_gRNAD,k_releaseD,d_assembledD,k_transWT_neg,k_transWT_pos,k_transD_neg,k_transD_pos), full_output = 1)

fig = plt.figure(figsize=(9,5),dpi=300)
ax = fig.add_subplot(111)
ax.plot(t2,Orig[:,17], label='without DIPS-$[V_{released}]$', color='orange')
#ax.fill_between(t, high_final, low_final, alpha=0.3, label='95% CI', color='orange')
ax.plot(t2, NEW[:,17], label='with DIPS-$[V_{released}^{wt}]$', color='darkorchid')
ax.plot(t2, NEW[:,20], label='with DIPS-$[V_{released}^{dip}]$', color='royalblue')
#ax.plot(t, (np.log(med_final2)-0), label='Median', color='darkorchid')
#ax.fill_between(t, high_final2, low_final2, alpha=0.3, label='95% CI', color='darkorchid')
#ax.grid('on')
plt.ylabel("fold change")
plt.yscale('symlog')
plt.legend(ncol=3,loc='upper center', bbox_to_anchor=(0.5, -0.10),fancybox=True, shadow=True,prop={'size': 14})
fig.savefig('V_released', bbox_inches='tight')
print(Orig[240,17],NEW[240,17])
print(Orig[480,17],NEW[480,17])



labels=['f_orf','f_N','f_SP','n_NW','k_bind','d_VW','k_diss','k_fuse','k_uncoat','d_endosomeW','d_gRNAW','$d_{V}^{dip}$','$d_{endosome}^{dip}$','$d_{gRNA}^{dip}$','k_transl','d_NSP','k_trWT_neg','K_NSP','d_gRNAWT_neg','k_trWT_pos','k_WTcomplex','$k_{tr(-)}^{dip}$','$d_{gRNA(-)}^{dip}$','$k_{tr(+)}^{dip}$','$k_{complex}^{dip}$','K_N','d_N','d_SP','n_SPWT','$n_{SP}^{dip}$','$n_{N}^{dip}$','K_VWT_rel','k_assemblWT','$K_{V_{rel}}^{dip}$','$k_{assembl}^{dip}$','d_NgRNAWT','k_releaseWT','d_assembledWT','$d_{N-gRNA}^{dip}$','$k_{release}^{dip}$','$d_{assembled}^{dip}$','$k_{trans(-)}^{wt}$','$k_{trans(+)}^{wt}$','$k_{trans(-)}^{dip}$','$k_{trans(+)}^{dip}$']
with open('summary.csv',"w",newline='') as csvfile:
        tofile = csv.writer(csvfile)
        tofile.writerow(["Parameter","median","mean","lower Confidence","Upper confidence","Smallest Distace"])
        for j in range(len(labels)):
            
            tofile.writerow([labels[j],np.median(sort[:,j]),np.mean(sort[:,j]),np.percentile(sort[:,j],2.5),np.percentile(sort[:,j],97.5),sort[0,j]])
                
csvfile.close()



    # d_VD=pow(10,np.random.uniform(-1.2,0.55))
    # d_endosomeD=pow(10,np.random.uniform(-4,-0.93))
    # d_gRNAD=pow(10,np.random.uniform(-1.16,-0.162))
    # k_trD_neg=pow(10,np.random.uniform(0,3))
    # d_gRNAD_neg=pow(10,np.random.uniform(-1.30,0))
    # k_trD_pos=pow(10,np.random.uniform(2.79,4.14))
    # k_Dcomplex=pow(10,np.random.uniform(-1.69,0))
    # n_SPD=pow(10,np.random.uniform(1,3.1))
    # n_ND=pow(10,np.random.uniform(1,2.35))
    # K_VD_rel=pow(10,np.random.uniform(1,4.31))
    # k_assemblD=pow(10,np.random.uniform(-2,1.31))
    # d_N_gRNAD=pow(10,np.random.uniform(-1.16,0))
    # k_releaseD=pow(10,np.random.uniform(0.9,3.15))
    # d_assembledD=pow(10,np.random.uniform(-4,-0.62))
    # k_transWT_neg=pow(10,np.random.uniform(-5,-3.7))
    # k_transWT_pos=pow(10,np.random.uniform(-2.21,-1.86))
    # k_transD_neg=pow(10,np.random.uniform(-5.69,-3))
    # k_transD_pos= pow(10,np.random.uniform(-2.9,-1.17))
    #plt.hist(sort[:,2], bins=np.logspace(np.log10(0.0001),np.log10(100), 20),color="darkorchid",alpha=0.4)
    #plt.hist(np.linspace(pow(10,-4),pow(10,2),_percent),alpha=0.5,bins=20)
    

_ignore=[0,1,2,3,4,5,6,7,8,9,10,14,15,16,17,18,19,20,25,26,27,28,31,32,35,36,37]

search_range=[[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],
                         [0,0],[0,0],[-1.2,0.55],[-4,-0.93],[-1.16,-0.162],[0,0],
                         [0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,3],[-1.3,0],[2.79,4.14],[-1.69,0],
                         [0,0],[0,0],[0,0],[0,0],[1,3.1],[1,2.35],[0,0],[0,0],[1,4.31],[-2,1.31],
                         [0,0],[0,0],[0,0],[-1.16,0],[0.9,3.14],[-4,-0.62],[-5,-3.7],[-2.22,-1.86]
                         ,[-5.69,-3],[-2.9,-1.17]]

for j in range(45):
    if j not in _ignore:
        fig=plt.figure(figsize=(9,5),dpi=300)
        plt.hist(sort[:,j], bins=np.linspace(pow(10,search_range[j][0]),pow(10,search_range[j][1]),10),color="darkorchid",alpha=1)
        plt.xlabel(labels[j])
        plt.axvline(np.median(sort[:,j]), color='black', linestyle='dashed', linewidth=2)
        plt.hist(np.linspace(pow(10,search_range[j][0]),pow(10,search_range[j][1]),_percent),alpha=0.5,bins=10)