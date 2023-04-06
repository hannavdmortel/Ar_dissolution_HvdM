import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker

filepath = "C:/Users/hanna/Documents/GitHub/Ar_dissolution_HvdM/Experiment_2/RADI"

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
fig, (ax1, ax2) = plt.subplots(2, 1,
                               dpi=300, 
                               figsize=(4, 5))

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
    ax1.scatter('_Omega_C.{}'.format(times[i]), 'Calcite_prod.{}'.format(times[i]), data = A_set[i], 
             c = colors[i], label=labels[i], alpha=0.8, s=15, edgecolor='none')

    #CUVB
    ax2.scatter('_Omega_C.{}'.format(times[i]), 'Calcite_prod.{}'.format(times[i]), data = B_set[i], 
             c = colors[i], label=labels[i], alpha=0.8, s=15, edgecolor='none')

#%%
#Obtain handles and labels from ax1
handles, labels = ax1.get_legend_handles_labels()

#Plot legend
leg = fig.legend(handles=handles, 
            bbox_to_anchor=(1.2, 0.65), 
            fontsize='x-small', 
            ncol=1,
            title='Hours elapsed:',
            title_fontsize='x-small')

leg._legend_box.align = "left"

#%%
#Axes CUVA
ax1.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.025))
ax1.yaxis.set_minor_locator(ticker.MultipleLocator(1))
ax1.yaxis.set_minor_locator(ticker.MultipleLocator(0.25))
ax1.grid(alpha=0.3, which='both')

#Axes CUVB
ax2.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
ax2.xaxis.set_minor_locator(ticker.MultipleLocator(0.025))
ax2.yaxis.set_minor_locator(ticker.MultipleLocator(1))
ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.25))
ax2.grid(alpha=0.3, which='both')

#Labels
#fig.suptitle("pH microprofiles", fontsize=16)
ax1.title.set_text('Cuvette A')
ax2.title.set_text('Cuvette B')

ax1.set_xlabel("1-$Ω_{ca}$")
ax1.set_ylabel('Calcite production \n($mol$ $m^{-3}$ $solid$ $yr^{-1}$)')
# ax1.axes.yaxis.set_ticklabels([])

ax2.set_xlabel("1-$Ω_{ca}$")
ax2.set_ylabel('Calcite production \n($mol$ $m^{-3}$ $solid$ $yr^{-1}$)')
# ax2.axes.yaxis.set_ticklabels([])


#%%
plt.tight_layout()

plt.subplots_adjust(wspace=None, hspace=None)
plt.savefig("figures/RADI_rates.png", bbox_inches='tight')
plt.show()

