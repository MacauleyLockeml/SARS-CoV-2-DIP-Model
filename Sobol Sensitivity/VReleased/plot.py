import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms
import seaborn as sns

#Reads in all of the data
S1 = np.loadtxt("S1.txt")
S1_conf = np.loadtxt("S1_conf.txt")
S1_upper = S1 + S1_conf
S1_lower = S1 - S1_conf
# S1=S1[0::10]
# S1_upper=S1_upper[0::10]
# S1_lower=S1_lower[0::10]
# S1_conf=S1_conf[0::10]

St = np.loadtxt("ST.txt")
St_conf = np.loadtxt("ST_conf.txt")
St_upper = St + St_conf
St_lower = St - St_conf
# St=St[0::10]
# St_upper=St_upper[0::10]
# St_lower=St_lower[0::10]
# St_conf=St_conf[0::10]

S1W = np.loadtxt("S1W.txt")
S1_confW = np.loadtxt("S1_confW.txt")
S1_upperW = S1W + S1_confW
S1_lowerW = S1W - S1_confW


StW = np.loadtxt("STW.txt")
St_confW = np.loadtxt("ST_confW.txt")
St_upperW = StW + St_confW
St_lowerW = StW - St_confW


S2=np.loadtxt("S2.txt")

#Timecourse to plot over
t = np.linspace(0,48,49)

#Labels for the parameters
labels_p = ['$k_{bind}$','$d_V^{wt}$','$k_{diss}$','$k_{fuse}$','$k_{uncoat}$','$d_{endosome}^{wt}$','$d_{gRNA}^{wt}$','$d_{V}^{dip}$','$d_{endosome}^{dip}$','$d_{gRNA}^{dip}$','$k_{transl}$','$d_{NSP}$','$k_{tr(-)}^{wt}$','$K_{NSP}$','$d_{gRNA(-)}^{wt}$','$k_{tr(+)}^{wt}$','$k_{complex}^{wt}$','$k_{tr(-)}^{dip}$','$d_{gRNA(-)}^{dip}$','$k_{tr(+)}^{dip}$','$k_{complex}^{dip}$','$K_{N}$','$d_{N}$','$d_{SP}$','$n_{SP}^{wt}$',	'$n_{SP}^{dip}$','$n_{N}^{dip}$','$K_{V_{rel}}^{wt}$','$k_{assembl}^{wt}$','$K_{V_{rel}}^{dip}$','$k_{assembl}^{dip}$','$d_{N-gRNA}^{wt}$','$k_{release}^{wt}$','$d_{assembled}^{wt}$','$d_{N-gRNA}^{dip}$','$k_{release}^{dip}$','$d_{assembled}^{dip}$','$k_{trans(-)}^{wt}$','$k_{trans(+)}^{wt}$','$k_{trans(-)}^{dip}$','$k_{trans(+)}^{dip}$']

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
axes[0].set(title='$[V_{released}^{dip}]$')
axes[0].set_xlabel("Sobol Index")
axes[0].set_xlim(xmax=1)

axes[0].legend(ncol=3,loc='upper center', bbox_to_anchor=(0.5, -0.10),fancybox=True, shadow=True,prop={'size': 14})
axes[1].barh(np.arange(len(labels_p)),StW[48,:], xerr=St_confW[48,:], label="Total order WT", align='center', color='royalblue', zorder=10)
axes[1].barh(np.arange(len(labels_p)),S1W[48,:],xerr=S1_confW[48,:],label=("First order WT"), align='center', color='orangered', zorder=10)
axes[1].set(title='$[V_{released}^{wt}]$')
axes[1].set_xlim(xmax=1)
axes[1].set_xlabel("Sobol Index")

axes[0].set(yticks=np.arange(len(labels_p)), yticklabels=labels_p)
axes[0].yaxis.tick_right()
# axes[0].set_xlim(xmax=1.25)
axes[0].invert_xaxis()


for ax in axes:
    ax.margins(0.03)
    ax.grid(False)

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
# sns.set(style="darkgrid")
# sns.set_color_codes("muted")
figall.savefig('Sobol_all_Vreleased', bbox_inches='tight')


# figs1, axes1 = plt.subplots(ncols=2, sharey=True, figsize=[8, 12],dpi=300)
# DIP1=axes1[0].barh(np.arange(len(labels_p)),S1[48,:],xerr=S1_conf[48,:], label="DIP", align='center', color='#43e653', zorder=10)

# axes1[0].set(title='$V_{released}^{dip}$')
# axes1[0].set_xlabel("First Order Sobol Index")


# axes1[1].barh(np.arange(len(labels_p)),S1W[48,:],xerr=S1_confW[48,:],label=("WT"), align='center', color='#ed1c3c', zorder=10)
# axes1[1].set(title='$V_{released}^{WT}$')
# axes1[1].set_xlim(xmin=-0.5,xmax=0.5)
# axes1[1].set_xlabel("First Order Sobol Index")

# axes1[0].set(yticks=np.arange(len(labels_p)), yticklabels=labels_p)
# axes1[0].yaxis.tick_right()
# axes1[0].set_xlim(xmin=-0.5,xmax=0.5)
# axes1[0].invert_xaxis()


# for ax in axes1:
#     ax.margins(0.03)
#     ax.grid(True)

# figs1.tight_layout()
# figs1.subplots_adjust(wspace=0.4)

# plt.setp(axes[0].yaxis.get_majorticklabels(), ha='center')

# # Create offset transform by some points in x direction
# dx = 36 / 72.
# dy = 0 / 72.
# offset = matplotlib.transforms.ScaledTranslation(dx, dy, figs1.dpi_scale_trans)
# # apply offset transform to all y ticklabels.
# for label in axes1[0].yaxis.get_majorticklabels():
#     label.set_transform(label.get_transform() + offset)
    
# figs1.savefig('S1_all_vs_released', bbox_inches='tight')

# print(np.sum(S1[48,:]))

# sns.set(style="darkgrid")
# sns.set_color_codes("muted")




# for i in range(len(labels_p)):
#     figl,axl=plt.subplots(figsize=(6,6),dpi=300)
#     axl.plot(t, S1[:,i], label=labels_p[i],linestyle='dashed')
#     #axl.fill_between(t, S1_lower[:,i], S1_upper[:,i], alpha=0.3)
#     axl.set(title=labels_p[i])
#     figl.savefig("DIP-{y}.png".format(y=labels_p[i]))
    
    
# for i in range(len(labels_p)):
#     figl,axl=plt.subplots(figsize=(6,6),dpi=300)
#     axl.plot(t, S1W[:,i], label=labels_p[i],linestyle='dashed')
#     axl.fill_between(t, S1_lowerW[:,i], S1_upperW[:,i], alpha=0.3)
#     axl.set(title=labels_p[i])
#     figl.savefig("WT-{y}.png".format(y=labels_p[i]))

# for i in range(len(labels_p)):
#     figl,axl=plt.subplots(figsize=(6,6),dpi=300)
#     axl.plot(t, St[:,i], label=labels_p[i],linestyle='dashed')
#     axl.fill_between(t, St_lower[:,i], St_upper[:,i], alpha=0.3)
#     axl.set(title=labels_p[i])
#     figl.savefig("DIP-{y}.png".format(y=labels_p[i]))

# def heatmap(data, row_labels, col_labels, ax=None,
#             cbar_kw={}, cbarlabel="", **kwargs):
#     """
#     Create a heatmap from a numpy array and two lists of labels.

#     Parameters
#     ----------
#     data
#         A 2D numpy array of shape (M, N).
#     row_labels
#         A list or array of length M with the labels for the rows.
#     col_labels
#         A list or array of length N with the labels for the columns.
#     ax
#         A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
#         not provided, use current axes or create a new one.  Optional.
#     cbar_kw
#         A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
#     cbarlabel
#         The label for the colorbar.  Optional.
#     **kwargs
#         All other arguments are forwarded to `imshow`.
#     """

#     if not ax:
#         ax = plt.gca()

#     # Plot the heatmap
#     im = ax.imshow(data, **kwargs)

#     # Create colorbar
#     cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
#     cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

#     # Show all ticks and label them with the respective list entries.
#     ax.set_xticks(np.arange(data.shape[1]), labels=col_labels)
#     ax.set_yticks(np.arange(data.shape[0]), labels=row_labels)

#     # Let the horizontal axes labeling appear on top.
#     ax.tick_params(top=True, bottom=False,
#                     labeltop=True, labelbottom=False)

#     # Rotate the tick labels and set their alignment.
#     plt.setp(ax.get_xticklabels(), rotation=-90, ha="right",
#               rotation_mode="anchor")

#     # Turn spines off and create white grid.
#     ax.spines[:].set_visible(False)

#     ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
#     ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
#     ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
#     ax.tick_params(which="minor", bottom=False, left=False)

#     return im, cbar
# def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
#                       textcolors=("black", "white"),
#                       threshold=None, **textkw):
#     """
#     A function to annotate a heatmap.

#     Parameters
#     ----------
#     im
#         The AxesImage to be labeled.
#     data
#         Data used to annotate.  If None, the image's data is used.  Optional.
#     valfmt
#         The format of the annotations inside the heatmap.  This should either
#         use the string format method, e.g. "$ {x:.2f}", or be a
#         `matplotlib.ticker.Formatter`.  Optional.
#     textcolors
#         A pair of colors.  The first is used for values below a threshold,
#         the second for those above.  Optional.
#     threshold
#         Value in data units according to which the colors from textcolors are
#         applied.  If None (the default) uses the middle of the colormap as
#         separation.  Optional.
#     **kwargs
#         All other arguments are forwarded to each call to `text` used to create
#         the text labels.
#     """

#     if not isinstance(data, (list, np.ndarray)):
#         data = im.get_array()

#     # Normalize the threshold to the images color range.
#     if threshold is not None:
#         threshold = im.norm(threshold)
#     else:
#         threshold = im.norm(data.max())/2.

#     # Set default alignment to center, but allow it to be
#     # overwritten by textkw.
#     kw = dict(horizontalalignment="center",
#               verticalalignment="center")
#     kw.update(textkw)

#     # Get the formatter in case a string is supplied
#     if isinstance(valfmt, str):
#         valfmt = plt.ticker.StrMethodFormatter(valfmt)

#     # Loop over the data and create a `Text` for each "pixel".
#     # Change the text's color depending on the data.
#     texts = []
#     for i in range(data.shape[0]):
#         for j in range(data.shape[1]):
#             kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
#             text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
#             texts.append(text)

#     return texts


# fig, ax = plt.subplots(figsize=(14,13),dpi=300)

# im, cbar = heatmap(S2, labels_p, labels_p, ax=ax,
#                     cmap="RdYlBu", cbarlabel="Sobol Index")
# texts = annotate_heatmap(im, valfmt="{x:.1f} t")
# ax.grid(False)
# fig.tight_layout()
# fig.savefig('S2', bbox_extra_artists=(lgd,), bbox_inches='tight')
