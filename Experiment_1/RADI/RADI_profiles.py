import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import AutoMinorLocator
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np

filepath = "C:/Users/hanna/Documents/GitHub/Ar_dissolution_HvdM/Experiment_1/RADI"

#%% 

#Append data of average pH's per 'set'
A_set = []
for i in range(0, 7):
    A_set.append(pd.read_excel(filepath + "/RADI_results.xlsx",
                          header=0, sheet_name='Cuvette A'))

B_set = []
for i in range(0, 7):
    B_set.append(pd.read_excel(filepath + "/RADI_results.xlsx",
                          header=0, sheet_name='Cuvette B'))

#%% Creating plot
fig, (ax) = plt.subplots(3, 2,
                               dpi=300, 
                               figsize=(5, 12))

#List colours for plot
colors = ['xkcd:bright lavender',
            'xkcd:soft blue',
            'xkcd:topaz',
            'xkcd:lightish green',
            'xkcd:yellowish',
            'xkcd:tangerine',
            'xkcd:watermelon']

#List labels for plot
labels =    ['1hr',
            '2hrs',
            '4hrs',
            '7hrs',
            '25hrs',
            '50hrs',
            '72hrs']

times = [1, 2, 4, 7, 25, 50, 70]

#Plotting each profile
for i in range(0, 7):
    #CUVA
    ax[0,0].plot('pH.{}'.format(times[i]), 'Depth (cm)', data = A_set[i], 
             c = colors[i], label=labels[ i], alpha=0.9, lw=1.5, zorder=10)
    ax[1,0].plot('Omega_C.{}'.format(times[i]), 'Depth (cm)', data = A_set[i], 
             c = colors[i], label=labels[i], alpha=0.9, lw=1.5, zorder=10)   
    ax[2,0].plot('Calcite_prod.{}'.format(times[i]), 'Depth (cm)', data = A_set[i],
                 c=colors[i], label=labels[i], alpha=0.9, lw=1.5, zorder=10)
    ax[2,0].plot('Ar_prod.{}'.format(times[i]), 'Depth (cm)', data = A_set[i],
                 c=colors[i], label=labels[i], alpha=0.9, lw=1.5, linestyle=':', zorder=10)
    #CUVB
    ax[0,1].plot('pH.{}'.format(times[i]), 'Depth (cm)', data = B_set[i], 
             c = colors[i], label=labels[i], alpha=0.9, lw=1.5, zorder=10)
    ax[1,1].plot('Omega_C.{}'.format(times[i]), 'Depth (cm)', data = B_set[i], 
             c = colors[i], label=labels[i], alpha=0.9, lw=1.5, zorder=10)
    ax[2,1].plot('Calcite_prod.{}'.format(times[i]), 'Depth (cm)', data = B_set[i],
                 c=colors[i], label=labels[i], alpha=0.9, lw=1.5, zorder=10)
    ax[2,1].plot('Ar_prod.{}'.format(times[i]), 'Depth (cm)', data = B_set[i],
                 c=colors[i], label=labels[i], alpha=0.9, lw=1.5, linestyle=':', zorder=10)   
#Add baselines
    ax[0,0].vlines(x=7.4, ymin=-3.25, ymax=1.2, color='grey', linestyle='--', alpha=0.8)
    ax[0,1].vlines(x=7.4, ymin=-3.25, ymax=1.2, color='grey', linestyle='--', alpha=0.8)
    ax[1,0].vlines(x=0.31, ymin=-3.25, ymax=1.2, color='grey', linestyle='--', alpha=0.8)
    ax[1,1].vlines(x=0.31, ymin=-3.25, ymax=1.2, color='grey', linestyle='--', alpha=0.8)

#%%
#Adding sediment and pteropod shading    
ax[0,0].axhspan(0, -2.9, color='xkcd:sand', alpha=0.15, lw=0)
ax[0,0].axhline(y=0, color='black', linewidth=0.5, linestyle='--', alpha=0.4)
ax[0,0].axhline(y=0.15, color='black', linewidth=0.5, linestyle=':', alpha=0.4)
ax[0,0].axhspan(0.15, 0, color='grey', alpha=0.3, lw=0)

ax[1,0].axhspan(0, -2.9, color='xkcd:sand', alpha=0.15, lw=0)
ax[1,0].axhline(y=0, color='black', linewidth=0.5, linestyle='--', alpha=0.4)
ax[1,0].axhline(y=0.15, color='black', linewidth=0.5, linestyle=':', alpha=0.4)
ax[1,0].axhspan(0.15, 0, color='grey', alpha=0.3, lw=0)

ax[2,0].axhspan(0, -2.9, color='xkcd:sand', alpha=0.15, lw=0)
ax[2,0].axhline(y=0, color='black', linewidth=0.5, linestyle='--', alpha=0.4)
ax[2,0].axhline(y=0.15, color='black', linewidth=0.5, linestyle=':', alpha=0.4)
ax[2,0].axhspan(0.15, 0, color='grey', alpha=0.3, lw=0)

ax[0,1].axhspan(0, -2.9, color='xkcd:sand', alpha=0.15, lw=0)
ax[0,1].axhline(y=0, color='black', linewidth=0.5, linestyle='--', alpha=0.4)

ax[1,1].axhspan(0, -2.9, color='xkcd:sand', alpha=0.15, lw=0)
ax[1,1].axhline(y=0, color='black', linewidth=0.5, linestyle='--', alpha=0.4)

ax[2,1].axhspan(0, -2.9, color='xkcd:sand', alpha=0.15, lw=0)
ax[2,1].axhline(y=0, color='black', linewidth=0.5, linestyle='--', alpha=0.4)
#%%
#Obtain handles and labels from ax1
handles, labels = ax[0,0].get_legend_handles_labels()

#Handles is a list, so append pteropod, SWI and calcite sand patches
patch0 = mpatches.Patch(color='white', alpha=0.3, edgecolor='white')
patch1 = mpatches.Patch(color='grey', label='Pteropods', alpha=0.3)
patch2 = mpatches.Patch(color='xkcd:sand', label='Calcite sand', alpha=0.15)
line1 =  Line2D([0], [0], color='black', alpha=0.4, linewidth=0.5, linestyle='--', label='SWI')
line2 =  Line2D([0], [0], color='grey', alpha=0.8, linewidth=1.5, linestyle='--', label='Baseline')

line3 = Line2D([0], [0], color='xkcd:watermelon', linewidth=1.5, linestyle='-', label='Calcite')
line4 = Line2D([0], [0], color='xkcd:watermelon', linewidth=1.5, linestyle=':', label='Aragonite')

for i in [line2, patch0, patch1, patch2, line1]:
    handles.append(i)  

#Plot legend
leg = fig.legend(handles=handles, 
            bbox_to_anchor=(1.18, 0.6), 
            fontsize='x-small', 
            ncol=1,
            title='Hours elapsed:',
            title_fontsize='x-small')

leg2 = fig.legend(handles=[line3, line4],
            bbox_to_anchor=(1.18, 0.22), 
            fontsize='x-small', 
            ncol=1,
            title='CaCO$_{3}$ polymorph:',
            title_fontsize='x-small')

leg._legend_box.align = "left"
leg2._legend_box.align = "left"

#%%
#Axes CUVA
ax[0,0].xaxis.set_major_locator(ticker.MultipleLocator(0.25))
ax[0,0].xaxis.set_minor_locator(AutoMinorLocator(5))
ax[0,0].yaxis.set_minor_locator(AutoMinorLocator(5))
ax[0,0].grid(alpha=0.3, which='both')
ax[0,0].set_xlim(7.2, 8)
ax[0,0].set_ylim(-2.9, 1.15)

ax[1,0].xaxis.set_major_locator(ticker.MultipleLocator(0.25))
ax[1,0].xaxis.set_minor_locator(AutoMinorLocator(5))
ax[1,0].yaxis.set_minor_locator(AutoMinorLocator(5))
ax[1,0].grid(alpha=0.3, which='both')
ax[1,0].set_xlim(0.15, 1.3)
ax[1,0].set_ylim(-2.9, 1.15)

ax[2,0].xaxis.set_major_locator(ticker.MultipleLocator(5))
ax[2,0].xaxis.set_minor_locator(AutoMinorLocator(5))
ax[2,0].yaxis.set_minor_locator(AutoMinorLocator(5))
ax[2,0].grid(alpha=0.3, which='both')
ax[2,0].set_xlim(-12, 8)
ax[2,0].set_ylim(-2.9, 1.15)

#Axes CUVB
ax[0,1].xaxis.set_major_locator(ticker.MultipleLocator(0.25))
ax[0,1].xaxis.set_minor_locator(AutoMinorLocator(5))
ax[0,1].yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax[0,1].grid(alpha=0.3, which='both')
ax[0,1].set_xlim(7.2, 8)
ax[0,1].set_ylim(-2.9, 1.15)

ax[1,1].xaxis.set_major_locator(ticker.MultipleLocator(0.25))
ax[1,1].xaxis.set_minor_locator(AutoMinorLocator(5))
ax[1,1].yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax[1,1].grid(alpha=0.3, which='both')
ax[1,1].set_xlim(0.15, 1.3)
ax[1,1].set_ylim(-2.9, 1.15)

ax[2,1].xaxis.set_major_locator(ticker.MultipleLocator(2.5))
ax[2,1].xaxis.set_minor_locator(AutoMinorLocator(5))
ax[2,1].yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax[2,1].grid(alpha=0.3, which='both')
ax[2,1].set_xlim(-5, 1)
ax[2,1].set_ylim(-2.9, 1.15)

#Labels
#fig.suptitle("pH microprofiles", fontsize=16)
ax[0,0].title.set_text('Cuvette A')
ax[0,1].title.set_text('Cuvette B')

ax[0,0].set_xlabel("$pH_{T}$")
ax[0,1].set_xlabel("$pH_{T}$")
ax[0,0].set_ylabel('Depth (cm)')
ax[0,1].axes.yaxis.set_ticklabels([])

ax[1,0].set_xlabel("$Ω_{ca}$")
ax[1,1].set_xlabel("$Ω_{ca}$")
ax[1,0].set_ylabel('Depth (cm)')
ax[1,1].axes.yaxis.set_ticklabels([])

ax[2,0].set_xlabel("CaCO$_3$ production \n($mol$ $m^{-3}$ $solid$ $yr^{-1}$)")
ax[2,1].set_xlabel("CaCO$_3$ production \n($mol$ $m^{-3}$ $solid$ $yr^{-1}$)")
ax[2,0].set_ylabel('Depth (cm)')
ax[2,1].axes.yaxis.set_ticklabels([])

#%%
plt.tight_layout()

plt.subplots_adjust(wspace=None, hspace=None)
plt.savefig("figures/RADI_profiles.png", bbox_inches='tight')
plt.show()

