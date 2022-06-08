from SALib.analyze import sobol
import numpy as np
import time
import pandas as pd
#These lines are for having run the code simultaneously on ARC
model = np.empty((0,490))
#
for i in range(40):
    data = np.loadtxt(open("Model_" + str(i+1) + ".csv", "rb"),delimiter=",",skiprows=0)
    model = np.vstack((model,data))

#model=np.loadtxt(open("Model.csv", "rb"), delimiter=",", skiprows=0)
#Restate the problem in this program, the same as in the sampling program
df=pd.read_csv('Param_values.csv', sep=',',header=None,index_col=(0))
bounds=df.values.tolist()

names=df.index.values

print(len(names))

problem = {
    'num_vars': len(names), 
    'names': names ,
    'bounds': bounds
}
print(len(model))
#Define empty datasets for first and total order Sobol indices. Could also
#include second order here. .
S1_df = np.empty((0,len(names)))
S1_conf_df = np.empty((0,len(names)))
ST_df = np.empty((0,len(names)))
ST_conf_df = np.empty((0,len(names)))


for i in range(49): #Length of your timecourse after cut
   Si_sob = sobol.analyze(problem, model[:,i], calc_second_order=True)

   S1 = Si_sob['S1']
   S1_df = np.vstack((S1_df,S1))

   S1_conf = Si_sob['S1_conf']
   S1_conf_df = np.vstack((S1_conf_df,S1_conf))
   
   S2=Si_sob['S2']

   ST = Si_sob['ST']
   ST_df = np.vstack((ST_df,ST))

   ST_conf = Si_sob['ST_conf']
   ST_conf_df = np.vstack((ST_conf_df,ST_conf))

#Saves to text files
np.savetxt("S1.txt", S1_df, fmt='%1.6f', delimiter=' ', newline='\n')
np.savetxt("S1_conf.txt", S1_conf_df, fmt='%1.6f', delimiter=' ', newline='\n')

np.savetxt("S2.txt", S2, fmt='%1.6f', delimiter=' ', newline='\n')

np.savetxt("ST.txt", ST_df, fmt='%1.6f', delimiter=' ', newline='\n')
np.savetxt("ST_conf.txt", ST_conf_df, fmt='%1.6f', delimiter=' ', newline='\n')

