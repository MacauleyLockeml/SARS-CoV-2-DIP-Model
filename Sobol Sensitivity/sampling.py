from SALib.sample import saltelli
import csv
import numpy as np
import pandas as pd
# Define the model inputs, sampling logarithmically through the ranges



df=pd.read_csv('Param_values.csv', sep=',',header=None,index_col=(0))
bounds=df.values.tolist()

names=df.index.values

problem = {
    'num_vars': len(names), 
    'names': names ,
    'bounds': bounds
}

# Generate samples
N=10000
D=len(names) #no. of parameters
param_values = saltelli.sample(problem, N,calc_second_order=True) #N*(D+2) 
#The following code separates the sample into chunks for running simultaneously on ARC

def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]

chunk = list(chunks(param_values,21000))

for i in range(40):
  chunki = chunk[i]
  np.savetxt("sample_" + str(i+1) + ".txt", chunki, fmt='%1.6f', delimiter=' ', newline='\n')

#This line saves the full set of parameter values as a text file
#np.savetxt("sample.txt", param_values, fmt='%1.6f', delimiter=' ', newline='\n')