import numpy as np
import matplotlib.pyplot as plt

#Reads in all of the data
S1 = np.loadtxt("S1.txt")
S1_conf = np.loadtxt("S1_conf.txt")
S1_upper = S1 + S1_conf
S1_lower = S1 - S1_conf

St = np.loadtxt("ST.txt")
St_conf = np.loadtxt("ST_conf.txt")
St_upper = St + St_conf
St_lower = St - St_conf

#Timecourse to plot over
t = np.linspace(0,48,49)

#Labels for the parameters
labels_p = ['$k_{bind}$','$d_V^{WT}$','$k_{diss}$','$k_{fuse}$','$k_{uncoat}$','$d_{endosome}^{WT}$','$d_{gRNA}^{WT}$','$d_{V}^{dip}$','$d_{endosome}^{dip}$','$d_{gRNA}^{dip}$','$k_{transl}$','$d_{NSP}$','$k_{tr(-)}^{WT}$','$K_{NSP}$','$d_{gRNA(-)}^{WT}$','$k_{tr(+)}^{WT}$','$k_{complex}^{WT}$','$k_{tr(-)}^{dip}$','$d_{gRNA(-)}^{dip}$','$k_{tr(+)}^{dip}$','$k_{complex}^{dip}$','$K_{N}$','$d_{N}$','$d_{SP}$','$n_{SP}^{WT}$',	'$n_{SP}^{dip}$','$n_{N}^{dip}$','$K_{V_{rel}}^{WT}$','$k_{assembl}^{WT}$','$K_{V_{rel}}^{dip}$','$k_{assembl}^{dip}$','$d_{N-gRNA}^{WT}$','$k_{release}^{WT}$','$d_{assembled}^{WT}$','$d_{N-gRNA}^{dip}$','$k_{release}^{dip}$','$d_{assembled}^{dip}$','$k_{trans(-)}^{WT}$','$k_{trans(+)}^{WT}$','$k_{tr(-)}^{dip}$','$k_{tr(+)}^{dip}$']

_common_index=[0,2,3,4,10,11,13,21,22,23]
_WT_index=[1,5,6,12,14,15,16,24,27,28,31,32,33]
_DIP_index=[7,8,9,17,18,19,20,25,26,29,30,34,35,36]

#Plotting the first order sobol indices with confidence intervals
f = plt.figure(figsize=(9,7),dpi=300)
ax = f.add_subplot(111)
count=0
for i in _common_index:
    if count<10:
        ax.plot(t, S1[:,i], label=labels_p[i],linestyle='dashed')
        count+1
        ax.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
    else:
        ax.plot(t, S1[:,i], label=labels_p[i],linestyle='dashed')
        
        ax.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
ax.grid('on')
handles, labels = ax.get_legend_handles_labels()
lgd = ax.legend(handles, labels, ncol=6,loc='upper center', bbox_to_anchor=(0.5, -0.10),fancybox=True, shadow=True, title='Parameter',prop={'size': 14})
plt.title('First order Sobol')
plt.ylabel('Sobol index')
plt.xlabel('Time (hours)')
plt.show()
f.savefig('S1_all_Common', bbox_extra_artists=(lgd,), bbox_inches='tight')


fw = plt.figure(figsize=(9,7),dpi=300)
axw = fw.add_subplot(111)
count=0
for j in _WT_index:
    if count<10:
        axw.plot(t, S1[:,j], label=labels_p[j],linestyle='dashed')
        count+1
        #ax.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
    else:
        axw.plot(t, S1[:,j], label=labels_p[j],linestyle='dashed')
        
        #ax.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
axw.grid('on')
handles, labels = axw.get_legend_handles_labels()
lgd = axw.legend(handles, labels, ncol=6,loc='upper center', bbox_to_anchor=(0.5, -0.10),fancybox=True, shadow=True, title='Parameter',prop={'size': 14})
plt.title('First order Sobol')
plt.ylabel('Sobol index')
plt.xlabel('Time (hours)')
plt.show()
f.savefig('S1_all_WT', bbox_extra_artists=(lgd,), bbox_inches='tight')



fd = plt.figure(figsize=(9,7),dpi=300)
axd = fd.add_subplot(111)
count=0
for j in _DIP_index:
    if count<10:
        axd.plot(t, S1[:,j], label=labels_p[j],linestyle='dashed')
        count+1
        #ax.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
    else:
        axd.plot(t, S1[:,j], label=labels_p[j],linestyle='dashed')
        
        #ax.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
axd.grid('on')
handles, labels = axd.get_legend_handles_labels()
lgd = axd.legend(handles, labels, ncol=6,loc='upper center', bbox_to_anchor=(0.5, -0.10),fancybox=True, shadow=True, title='Parameter',prop={'size': 14})
plt.title('First order Sobol')
plt.ylabel('Sobol index')
plt.xlabel('Time (hours)')
plt.show()
f.savefig('S1_all', bbox_extra_artists=(lgd,), bbox_inches='tight')


#Plotting the total order sobol indices with confidence intervals
f2 = plt.figure(figsize=(9,7), dpi=120) #dpi increases the quality and size of the image
ax2 = f2.add_subplot(111)
for i in range(len(labels)):
    if i<10:
        ax2.plot(t, St[:,i], label=labels[i])
        #ax2.fill_between(t, St_lower[:,i], St_upper[:,i], alpha=0.3)
    else:
        ax2.plot(t, St[:,i], label=labels[i],linestyle='dashed')
        #ax2.fill_between(t, St_lower[:,i], St_upper[:,i], alpha=0.3)
        
ax2.grid('on')
handles2, labels2 = ax2.get_legend_handles_labels()
lgd = ax2.legend(handles2, labels2,ncol=6,loc='upper center', bbox_to_anchor=(0.5, -0.10),fancybox=True, shadow=True, title='Parameter',prop={'size': 14})
#plt.title('Total order Sobol indicies for timecourse of activated TBK-1 for model 1')
plt.ylabel('Sobol index',fontsize=14)
plt.yticks(fontsize=13)
plt.xlabel('Time (hours)',fontsize=14)
plt.xticks(fontsize=13)
plt.show()
f2.savefig('ST_all_model_1', bbox_extra_artists=(lgd,), bbox_inches='tight')

# #Individual subplots for each of the parameters. For the first order indices
# f3, ax3 = plt.subplots(3,2,figsize=(12,8))
# plt.subplots_adjust(left=None, bottom=None, right=None, top=1.2, wspace=None, hspace=0.4)

# plt.subplot(3,2,1)
# plt.plot(t, S1[:,10], label=labels[10])
# plt.fill_between(t, S1_lower[:,10], S1_upper[:,10], alpha=0.3)
# plt.grid()
# plt.xlabel('Time (hours)')
# plt.ylabel('Sobol index')
# plt.title('Parameter: $n_{B}$')

# plt.subplot(3,2,2)
# plt.plot(t, S1[:,4], label=labels[4])
# plt.fill_between(t, S1_lower[:,4], S1_upper[:,4], alpha=0.3)
# plt.grid()
# plt.xlabel('Time (hours)')
# plt.ylabel('Sobol index')
# plt.title('Parameter: $k_{B}$')

# plt.subplot(3,2,3)
# plt.plot(t, S1[:,6], label=labels[6])
# plt.fill_between(t, S1_lower[:,6], S1_upper[:,6], alpha=0.3)
# plt.grid()
# plt.xlabel('Time (hours)')
# plt.ylabel('Sobol index')
# plt.title('Parameter: $\kappa_{V}$')

# plt.subplot(3,2,4)
# plt.plot(t, S1[:,8], label=labels[8])
# plt.fill_between(t, S1_lower[:,8], S1_upper[:,8], alpha=0.3)
# plt.grid()
# plt.xlabel('Time (hours)')
# plt.ylabel('Sobol index')
# plt.title('Parameter: $n_{R}$')

# plt.subplot(3,2,5)
# plt.plot(t, S1[:,7], label=labels[7])
# plt.fill_between(t, S1_lower[:,7], S1_upper[:,7], alpha=0.3)
# plt.grid()
# plt.xlabel('Time (hours)')
# plt.ylabel('Sobol index')
# plt.title('Parameter: $n_{D}$')

# plt.subplot(3,2,6)
# plt.plot(t, S1[:,9], label=labels[9])
# plt.fill_between(t, S1_lower[:,9], S1_upper[:,9], alpha=0.3)
# plt.grid()
# plt.xlabel('Time (hours)')
# plt.ylabel('Sobol index')
# plt.title('Parameter: $n_{V}$')

# plt.suptitle('Timecourses 1-6 of the first order Sobol indices from most to least important \n at the end of the timecourse', y=1.3, fontsize=15)
# plt.savefig('S1_ind_top.png', bbox_inches='tight')
# plt.show()



# f4, ax4=plt.subplots(3,2,figsize=(12,8))
# plt.subplots_adjust(left=None, bottom=None, right=None, top=1.2, wspace=None, hspace=0.4)

# plt.subplot(3,2,1)
# plt.plot(t, S1[:,0], label=labels[0])
# plt.fill_between(t, S1_lower[:,0], S1_upper[:,0], alpha=0.3)
# plt.grid()
# plt.xlabel('Time (hours)')
# plt.ylabel('Sobol index')
# plt.title('Parameter: $k_{R}$')

# plt.subplot(3,2,2)
# plt.plot(t, S1[:,1], label=labels[1])
# plt.fill_between(t, S1_lower[:,1], S1_upper[:,1], alpha=0.3)
# plt.grid()
# plt.xlabel('Time (hours)')
# plt.ylabel('Sobol index')
# plt.title('Parameter: $q_{R}$')

# plt.subplot(3,2,3)
# plt.plot(t, S1[:,2], label=labels[2])
# plt.fill_between(t, S1_lower[:,2], S1_upper[:,2], alpha=0.3)
# plt.grid()
# plt.xlabel('Time (hours)')
# plt.ylabel('Sobol index')
# plt.title('Parameter: $k_{V}$')

# plt.subplot(3,2,4)
# plt.plot(t, S1[:,3], label=labels[3])
# plt.fill_between(t, S1_lower[:,3], S1_upper[:,3], alpha=0.3)
# plt.grid()
# plt.xlabel('Time (hours)')
# plt.ylabel('Sobol index')
# plt.title('Parameter: $q_{V}$')

# plt.subplot(3,2,5)
# plt.plot(t, S1[:,5], label=labels[5])
# plt.fill_between(t, S1_lower[:,5], S1_upper[:,5], alpha=0.3)
# plt.grid()
# plt.xlabel('Time (hours)')
# plt.ylabel('Sobol index')
# plt.title('Parameter: $q_{B}$')

# # plt.subplot(4,2,6)
# # plt.plot(t, S1[:,11], label=labels[11])
# # plt.fill_between(t, S1_lower[:,11], S1_upper[:,11], alpha=0.3)
# # plt.grid()
# # plt.xlabel('Time (hours)')
# # plt.ylabel('Sobol index')
# # plt.title('Parameter: $n_{1}0$')

# # plt.subplot(4,2,7)
# # plt.plot(t, S1[:,12], label=labels[12])
# # plt.fill_between(t, S1_lower[:,12], S1_upper[:,12], alpha=0.3)
# # plt.grid()
# # plt.xlabel('Time (hours)')
# # plt.ylabel('Sobol index')
# # plt.title('Parameter: $n_{2}0$')

# # plt.subplot(4,2,8)
# # plt.plot(t, S1[:,13], label=labels[13])
# # plt.fill_between(t, S1_lower[:,13], S1_upper[:,13], alpha=0.3)
# # plt.grid()
# # plt.xlabel('Time (hours)')
# # plt.ylabel('Sobol index')
# # plt.title('Parameter: $n_{3}0$')

# plt.suptitle('Timecourses of the first order Sobol indices for remaing parameters \n at the end of the timecourse', y=1.3, fontsize=15)
# plt.savefig('S1_ind.png', bbox_inches='tight')
# plt.show()

# #Individual subplots for each of the parameters. For the total order indices
# f5, ax5 = plt.subplots(3,2,figsize=(12,8))
# plt.subplots_adjust(left=None, bottom=None, right=None, top=1.2, wspace=None, hspace=0.4)

# plt.subplot(3,2,1)
# plt.plot(t, St[:,10], label=labels[10], color='green')
# plt.fill_between(t, St_lower[:,10], St_upper[:,10], alpha=0.3, color='green')
# plt.grid()
# plt.xlabel('Time (hours)')
# plt.ylabel('Sobol index')
# plt.title('Parameter: $n_{B}$')

# plt.subplot(3,2,2)
# plt.plot(t, St[:,4], label=labels[4], color='green')
# plt.fill_between(t, St_lower[:,4], St_upper[:,4], alpha=0.3, color='green')
# plt.grid()
# plt.xlabel('Time (hours)')
# plt.ylabel('Sobol index')
# plt.title('Parameter: $k_{b}$')

# plt.subplot(3,2,3)
# plt.plot(t, St[:,6], label=labels[6], color='green')
# plt.fill_between(t, St_lower[:,6], St_upper[:,6], alpha=0.3, color='green')
# plt.grid()
# plt.xlabel('Time (hours)')
# plt.ylabel('Sobol index')
# plt.title('Parameter: $\kappa_{v}$')

# plt.subplot(3,2,4)
# plt.plot(t, St[:,7], label=labels[7], color='green')
# plt.fill_between(t, St_lower[:,7], St_upper[:,7], alpha=0.3, color='green')
# plt.grid()
# plt.xlabel('Time (hours)')
# plt.ylabel('Sobol index')
# plt.title('Parameter: $n_{D}$')

# plt.subplot(3,2,5)
# plt.plot(t, St[:,5], label=labels[5], color='green')
# plt.fill_between(t, St_lower[:,5], St_upper[:,5], alpha=0.3, color='green')
# plt.grid()
# plt.xlabel('Time (hours)')
# plt.ylabel('Sobol index')
# plt.title('Parameter: $q_{b}$')

# plt.subplot(3,2,6)
# plt.plot(t, St[:,8], label=labels[8], color='green')
# plt.fill_between(t, St_lower[:,8], St_upper[:,8], alpha=0.3, color='green')
# plt.grid()
# plt.xlabel('Time (hours)')
# plt.ylabel('Sobol index')
# plt.title('Parameter: $n_{R}$')

# plt.suptitle('Total order Sobol indicies for the parameters in the model to the timecourse of $n_{3}$, \n 6 most important', y=1.3, fontsize=15)
# plt.savefig('ST_ind_top.png', bbox_inches='tight')
# plt.show()


# f6, ax6=plt.subplots(3,2,figsize=(12,8))
# plt.subplots_adjust(left=None, bottom=None, right=None, top=1.2, wspace=None, hspace=0.4)

# plt.subplot(3,2,1)
# plt.plot(t, St[:,0], label=labels[0], color='green')
# plt.fill_between(t, St_lower[:,0], St_upper[:,0], alpha=0.3, color='green')
# plt.grid()
# plt.xlabel('Time (hours)')
# plt.ylabel('Sobol index')
# plt.title('Parameter: $k_{r}$')

# plt.subplot(3,2,2)
# plt.plot(t, St[:,1], label=labels[1], color='green')
# plt.fill_between(t, St_lower[:,1], St_upper[:,1], alpha=0.3, color='green')
# plt.grid()
# plt.xlabel('Time (hours)')
# plt.ylabel('Sobol index')
# plt.title('Parameter: $q_{r}$')

# plt.subplot(3,2,3)
# plt.plot(t, St[:,2], label=labels[2], color='green')
# plt.fill_between(t, St_lower[:,2], St_upper[:,2], alpha=0.3, color='green')
# plt.grid()
# plt.xlabel('Time (hours)')
# plt.ylabel('Sobol index')
# plt.title('Parameter: $k_{v}$')

# plt.subplot(3,2,4)
# plt.plot(t, St[:,3], label=labels[3], color='green')
# plt.fill_between(t, St_lower[:,3], St_upper[:,3], alpha=0.3, color='green')
# plt.grid()
# plt.xlabel('Time (hours)')
# plt.ylabel('Sobol index')
# plt.title('Parameter: $q_{v}$')

# plt.subplot(3,2,5)
# plt.plot(t, St[:,9], label=labels[9], color='green')
# plt.fill_between(t, St_lower[:,9], St_upper[:,9], alpha=0.3, color='green')
# plt.grid()
# plt.xlabel('Time (hours)')
# plt.ylabel('Sobol index')
# plt.title('Parameter: $q_{b}$')

# # plt.subplot(4,2,6)
# # plt.plot(t, St[:,11], label=labels[11], color='green')
# # plt.fill_between(t, St_lower[:,11], St_upper[:,11], alpha=0.3, color='green')
# # plt.grid()
# # plt.xlabel('Time (hours)')
# # plt.ylabel('Sobol index')
# # plt.title('Parameter: $n_{1}0$')

# # plt.subplot(4,2,7)
# # plt.plot(t, St[:,12], label=labels[12], color='green')
# # plt.fill_between(t, St_lower[:,12], St_upper[:,12], alpha=0.3, color='green')
# # plt.grid()
# # plt.xlabel('Time (hours)')
# # plt.ylabel('Sobol index')
# # plt.title('Parameter: $n_{2}0$')

# # plt.subplot(4,2,8)
# # plt.plot(t, St[:,13], label=labels[13], color='green')
# # plt.fill_between(t, St_lower[:,13], St_upper[:,13], alpha=0.3, color='green')
# # plt.grid()
# # plt.xlabel('Time (hours)')
# # plt.ylabel('Sobol index')
# # plt.title('Parameter: $n_{3}0$')

# plt.suptitle('Total order Sobol indicies for the parameters in the model to the timecourse of $n_{3}$, \n remaining parameters', y=1.3, fontsize=15)
# plt.savefig('St_ind.png', bbox_inches='tight')
# plt.show()