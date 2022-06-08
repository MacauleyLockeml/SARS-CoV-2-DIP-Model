import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms
import seaborn as sns

#Reads in all of the data
S1 = np.loadtxt("S1.txt")
S1_conf = np.loadtxt("S1_conf.txt")
S1_upper = S1 + S1_conf
S1_lower = S1 - S1_conf


St = np.loadtxt("ST.txt")
St_conf = np.loadtxt("ST_conf.txt")
St_upper = St + St_conf
St_lower = St - St_conf


S1W = np.loadtxt("S1W.txt")
S1_confW = np.loadtxt("S1_confW.txt")
S1_upperW = S1W + S1_confW
S1_lowerW = S1W - S1_confW


StW = np.loadtxt("STW.txt")
St_confW = np.loadtxt("ST_confW.txt")
St_upperW = StW + St_confW
St_lowerW = StW - St_confW


#Timecourse to plot over
t = np.linspace(0,48,49)

#Labels for the parameters
labels_p = ['$k_{bind}$','$d_V^{WT}$','$k_{diss}$','$k_{fuse}$','$k_{uncoat}$','$d_{endosome}^{WT}$','$d_{gRNA}^{WT}$','$d_{V}^{dip}$','$d_{endosome}^{dip}$','$d_{gRNA}^{dip}$','$k_{transl}$','$d_{NSP}$','$k_{tr(-)}^{WT}$','$K_{NSP}$','$d_{gRNA(-)}^{WT}$','$k_{tr(+)}^{WT}$','$k_{complex}^{WT}$','$k_{tr(-)}^{dip}$','$d_{gRNA(-)}^{dip}$','$k_{tr(+)}^{dip}$','$k_{complex}^{dip}$','$K_{N}$','$d_{N}$','$d_{SP}$','$n_{SP}^{WT}$',	'$n_{SP}^{dip}$','$n_{N}^{dip}$','$K_{V_{rel}}^{WT}$','$k_{assembl}^{WT}$','$K_{V_{rel}}^{dip}$','$k_{assembl}^{dip}$','$d_{N-gRNA}^{WT}$','$k_{release}^{WT}$','$d_{assembled}^{WT}$','$d_{N-gRNA}^{dip}$','$k_{release}^{dip}$','$d_{assembled}^{dip}$','$k_{trans(-)}^{WT}$','$k_{trans(+)}^{WT}$','$k_{tr(-)}^{dip}$','$k_{tr(+)}^{dip}$']

_common_index=[0,2,3,4,10,11,13,21,22,23]
_WT_index=[1,5,6,12,14,15,16,24,27,28,31,32,33]
_DIP_index=[7,8,9,17,18,19,20,25,26,29,30,34,35,36]

# #Plotting the first order sobol indices with confidence intervals
# f = plt.figure(figsize=(9,7),dpi=300)
# ax = f.add_subplot(111)
# count=0
# for i in _common_index:
#     if count<10:
#         ax.plot(t, S1[:,i], label=labels_p[i],linestyle='dashed')
#         count+1
#         #ax.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
#     else:
#         ax.plot(t, S1[:,i], label=labels_p[i],linestyle='dashed')
        
#         #ax.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
# ax.grid('on')
# handles, labels = ax.get_legend_handles_labels()
# lgd = ax.legend(handles, labels, ncol=6,loc='upper center', bbox_to_anchor=(0.5, -0.10),fancybox=True, shadow=True, title='Parameter',prop={'size': 14})
# plt.title('First order Sobol')
# plt.ylabel('Sobol index')
# plt.xlabel('Time (hours)')
# plt.show()
# f.savefig('S1_all_Common', bbox_extra_artists=(lgd,), bbox_inches='tight')


# fw = plt.figure(figsize=(9,7),dpi=300)
# axw = fw.add_subplot(111)
# count=0
# for j in _WT_index:
#     if count<10:
#         axw.plot(t, S1[:,j], label=labels_p[j],linestyle='dashed')
#         count+1
#         #ax.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
#     else:
#         axw.plot(t, S1[:,j], label=labels_p[j],linestyle='dashed')
        
#         #ax.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
# axw.grid('on')
# handles, labels = axw.get_legend_handles_labels()
# lgd = axw.legend(handles, labels, ncol=6,loc='upper center', bbox_to_anchor=(0.5, -0.10),fancybox=True, shadow=True, title='Parameter',prop={'size': 14})
# plt.title('First order Sobol')
# plt.ylabel('Sobol index')
# plt.xlabel('Time (hours)')
# plt.show()
# f.savefig('S1_all_WT', bbox_extra_artists=(lgd,), bbox_inches='tight')



# fd = plt.figure(figsize=(9,7),dpi=300)
# axd = fd.add_subplot(111)
# count=0
# for j in _DIP_index:
#     if count<10:
#         axd.plot(t, S1[:,j], label=labels_p[j],linestyle='dashed')
#         count+1
#         #ax.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
#     else:
#         axd.plot(t, S1[:,j], label=labels_p[j],linestyle='dashed')
        
#         #ax.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
# axd.grid('on')
# handles, labels = axd.get_legend_handles_labels()
# lgd = axd.legend(handles, labels, ncol=6,loc='upper center', bbox_to_anchor=(0.5, -0.10),fancybox=True, shadow=True, title='Parameter',prop={'size': 14})
# plt.title('First order Sobol')
# plt.ylabel('Sobol index')
# plt.xlabel('Time (hours)')
# plt.show()
# f.savefig('S1_all', bbox_extra_artists=(lgd,), bbox_inches='tight')


# #Plotting the total order sobol indices with confidence intervals
# f2 = plt.figure(figsize=(9,7), dpi=120) #dpi increases the quality and size of the image
# ax2 = f2.add_subplot(111)
# for i in range(len(labels)):
#     if i<10:
#         ax2.plot(t, St[:,i], label=labels[i])
#         #ax2.fill_between(t, St_lower[:,i], St_upper[:,i], alpha=0.3)
#     else:
#         ax2.plot(t, St[:,i], label=labels[i],linestyle='dashed')
#         #ax2.fill_between(t, St_lower[:,i], St_upper[:,i], alpha=0.3)
        
# ax2.grid('on')
# handles2, labels2 = ax2.get_legend_handles_labels()
# lgd = ax2.legend(handles2, labels2,ncol=6,loc='upper center', bbox_to_anchor=(0.5, -0.10),fancybox=True, shadow=True, title='Parameter',prop={'size': 14})
# #plt.title('Total order Sobol indicies for timecourse of activated TBK-1 for model 1')
# plt.ylabel('Sobol index',fontsize=14)
# plt.yticks(fontsize=13)
# plt.xlabel('Time (hours)',fontsize=14)
# plt.xticks(fontsize=13)
# plt.show()
# f2.savefig('ST_all_model_1', bbox_extra_artists=(lgd,), bbox_inches='tight')


# _index=ind = np.arange(len(labels_p))
# width=0.5

# fig, ax = plt.subplots(figsize=(6,20),dpi=300)
# ax.barh( np.arange(len(labels_p)),St[24,:], label="DIP" )
# ax.barh(np.arange(len(labels_p)),StW[24,:],label=("WT"),color="orange" )
# ax.set_yticks(ind)
# ax.set_yticklabels(labels_p,fontsize=20)
# plt.xticks(fontsize=20,rotation=90)
# plt.yticks(fontsize=20)
# ax.set_ylabel('Total Sobol index',fontsize=20)
# ax.set_xlabel('Parameters',fontsize=20)
# ax.legend()


figall, axes = plt.subplots(ncols=2, sharey=True, figsize=[8, 12],dpi=300)
DIPp=axes[0].barh(np.arange(len(labels_p)),St[48,:],xerr=St_conf[48,:], label="Total order DIP", align='center', color='darkorange', zorder=10)
axes[0].barh(np.arange(len(labels_p)),S1[48,:],xerr=S1_conf[48,:], label="First order DIP", align='center', color='mediumpurple', zorder=10)
axes[0].set(title='$N-gRNA^{dip}$')
axes[0].set_xlabel("Sobol Index")
axes[0].set_xlim(xmax=1)

axes[0].legend(ncol=3,loc='upper center', bbox_to_anchor=(0.5, -0.10),fancybox=True, shadow=True,prop={'size': 14})
axes[1].barh(np.arange(len(labels_p)),StW[48,:], xerr=St_confW[48,:], label="Total order WT", align='center', color='royalblue', zorder=10)
axes[1].barh(np.arange(len(labels_p)),S1W[48,:],xerr=S1_confW[48,:],label=("First order WT"), align='center', color='orangered', zorder=10)
axes[1].set(title='$N-gRNA^{WT}$')
axes[1].set_xlim(xmax=1)
axes[1].set_xlabel("Sobol Index")

axes[0].set(yticks=np.arange(len(labels_p)), yticklabels=labels_p)
axes[0].yaxis.tick_right()
# axes[0].set_xlim(xmax=1.25)
axes[0].invert_xaxis()


for ax in axes:
    ax.margins(0.03)
    ax.grid(True)

figall.tight_layout()
figall.subplots_adjust(wspace=0.4)

plt.setp(axes[0].yaxis.get_majorticklabels(), ha='center')

# Create offset transform by some points in x direction
dx = 36 / 72.
dy = 0 / 72.
offset = matplotlib.transforms.ScaledTranslation(dx, dy, figall.dpi_scale_trans)
# apply offset transform to all y ticklabels.
for label in axes[0].yaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)
    
axes[0].legend(ncol=3,loc='upper center', bbox_to_anchor=(1.25, -0.05),fancybox=True, shadow=True,prop={'size': 14})
axes[1].legend(ncol=3,loc='upper center', bbox_to_anchor=(-0.15, -0.10),fancybox=True, shadow=True,prop={'size': 14})
    
figall.savefig('Sobol_all', bbox_inches='tight')





sns.set(style="darkgrid")
sns.set_color_codes("muted")


