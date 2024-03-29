import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import AutoMinorLocator
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np
from scipy.interpolate import interp1d

filepath = "C:/Users/hanna/Documents/GitHub/Ar_dissolution_HvdM/Experiment_2/Lab"


#%% First for cuvette A:

#Append data of average pH's per 'set'
A_set = []
for i in range(0, 7):
    A_set.append(pd.read_excel(filepath + "/CUV_A.xlsx",
                          header=0, sheet_name='AVG{}'.format(i+1)))

    #Logical for standard deviation < 0.02
    L_stddev = (A_set[i]['Std.dev1']<0.02)&(A_set[i]['Std.dev2']<0.02)
    A_set[i] = A_set[i][L_stddev]
    
    #Create columns with average of remaining pH and std dev
    A_set[i]['Average pH'] = A_set[i][['pH1', 'pH2']].mean(axis=1)
    A_set[i]['[H+]'] = (10**9)*(10**(-A_set[i]['Average pH']))
    
    #For TRIS and std dev microelectrode
    A_set[i]['Total std dev'] = np.sqrt(((A_set[i][['Std.dev1', 'Std.dev2']].mean(axis=1))**2)+0.006**2)
    A_set[i]['H+ std dev'] = -A_set[i]['Total std dev']*A_set[i]['[H+]']*np.log(10)
#%% Then for cuvette B:

#Append data of average pH's per 'set'
B_set = []
for i in range(0, 5):
    B_set.append(pd.read_excel(filepath + "/CUV_B.xlsx",
                          header=0, sheet_name='AVG{}'.format(i+1)))

    #Logical for standard deviation < 0.02
    L_stddev = (B_set[i]['Std.dev1']<0.02)&(B_set[i]['Std.dev2']<0.02)
    B_set[i] = B_set[i][L_stddev]
    
    #Create columns with average of remaining pH and std dev
    B_set[i]['Average pH'] = B_set[i][['pH1', 'pH2']].mean(axis=1)
    B_set[i]['[H+]'] = (10**9)*(10**(-B_set[i]['Average pH']))

    #For TRIS and microelectrode std dev
    B_set[i]['Total std dev'] = np.sqrt(((B_set[i][['Std.dev1', 'Std.dev2']].mean(axis=1))**2)+0.006**2)
    B_set[i]['H+ std dev'] = -B_set[i]['Total std dev']*B_set[i]['[H+]']*np.log(10)

#%% Creating plot
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2,
                               dpi=300, 
                               figsize=(5, 10))

#List colours for plot
colors_A = ['xkcd:bright lavender',
            'xkcd:soft blue',
            'xkcd:topaz',
            'xkcd:lightish green',
            'xkcd:yellowish',
            'xkcd:tangerine',
            'xkcd:watermelon']

colors_B = ['xkcd:bright lavender',
            'xkcd:soft blue',
            'xkcd:topaz',
            'xkcd:yellowish',
            'xkcd:tangerine']

#List labels for plot
labels_A = ['1hr',
            '2hrs',
            '4hrs',
            '7hrs',
            '25hrs',
            '50hrs',
            '72hrs']


#Plotting each profile as a line plot
#CUV A
for i in range(0, 7):
    #Alter depth to line up SWI up both cuvettes
    A_set[i]['Depth corrected'] = A_set[i]['Depth (cm)']+0.1
    #Put 0 reading at SWI, and negative going into sediment
    A_set[i]['Depth corrected'] = (A_set[i]['Depth corrected']-1.15)
    
    #Plot average pH
    ax1.plot('Average pH', 'Depth corrected', data = A_set[i], 
              c = colors_A[i], label=labels_A[i], alpha=0.9, lw=1.5)
    #Plot [H+]
    ax3.plot('[H+]', 'Depth corrected', data = A_set[i], 
             c = colors_A[i], label=labels_A[i], alpha=0.9, lw=1.5)
    
    #Do error shading in little steps to avoid incorrect shading
    A_set[i]['Depth +0.01'] = (A_set[i]['Depth corrected']+0.01)
    A_set[i]['Depth +0.02'] = (A_set[i]['Depth corrected']+0.02)
    A_set[i]['Depth +0.03'] = (A_set[i]['Depth corrected']+0.03)
    A_set[i]['Depth +0.04'] = (A_set[i]['Depth corrected']+0.04)
    A_set[i]['Depth +0.05'] = (A_set[i]['Depth corrected']+0.05)
    A_set[i]['Depth +0.06'] = (A_set[i]['Depth corrected']+0.06)
    A_set[i]['Depth +0.07'] = (A_set[i]['Depth corrected']+0.07)
    A_set[i]['Depth +0.08'] = (A_set[i]['Depth corrected']+0.08)
    A_set[i]['Depth +0.09'] = (A_set[i]['Depth corrected']+0.09)
    A_set[i]['Depth +0.1'] = (A_set[i]['Depth corrected']+0.1)
    
    A_set[i]['Depth -0.01'] = (A_set[i]['Depth corrected']-0.01)
    A_set[i]['Depth -0.02'] = (A_set[i]['Depth corrected']-0.02)
    A_set[i]['Depth -0.03'] = (A_set[i]['Depth corrected']-0.03)
    A_set[i]['Depth -0.04'] = (A_set[i]['Depth corrected']-0.04)
    A_set[i]['Depth -0.05'] = (A_set[i]['Depth corrected']-0.05)
    A_set[i]['Depth -0.06'] = (A_set[i]['Depth corrected']-0.06)
    A_set[i]['Depth -0.07'] = (A_set[i]['Depth corrected']-0.07)
    A_set[i]['Depth -0.08'] = (A_set[i]['Depth corrected']-0.08)
    A_set[i]['Depth -0.09'] = (A_set[i]['Depth corrected']-0.09)
    A_set[i]['Depth -0.1'] = (A_set[i]['Depth corrected']-0.1)
    
    depths_A_up = [A_set[i]['Depth corrected'],
                  A_set[i]['Depth +0.01'], 
                  A_set[i]['Depth +0.02'],
                  A_set[i]['Depth +0.03'],
                  A_set[i]['Depth +0.04'],
                  A_set[i]['Depth +0.05'],
                  A_set[i]['Depth +0.06'],
                  A_set[i]['Depth +0.07'],
                  A_set[i]['Depth +0.08'],
                  A_set[i]['Depth +0.09'],
                  A_set[i]['Depth +0.1']]
    
    for j in range(len(depths_A_up)-1):
        ax1.fill_between(x = A_set[i]['Average pH'],
                      y1 = depths_A_up[j], 
                      y2 = depths_A_up[j+1], 
                      color=colors_A[i], edgecolor='none', alpha=0.4)
        ax3.fill_between(x = A_set[i]['[H+]'],
                      y1 = depths_A_up[j], 
                      y2 = depths_A_up[j+1], 
                      color=colors_A[i], edgecolor='none', alpha=0.4)
        
    depths_A_down = [A_set[i]['Depth corrected'],
                  A_set[i]['Depth -0.01'], 
                  A_set[i]['Depth -0.02'],
                  A_set[i]['Depth -0.03'],
                  A_set[i]['Depth -0.04'],
                  A_set[i]['Depth -0.05'],
                  A_set[i]['Depth -0.06'],
                  A_set[i]['Depth -0.07'],
                  A_set[i]['Depth -0.08'],
                  A_set[i]['Depth -0.09'],
                  A_set[i]['Depth -0.1']]
    
    for j in range(len(depths_A_down)-1):
        ax1.fill_between(x = A_set[i]['Average pH'],
                      y1 = depths_A_down[j], 
                      y2 = depths_A_down[j+1], 
                      color=colors_A[i], edgecolor='none', alpha=0.4)
        ax3.fill_between(x = A_set[i]['[H+]'],
                      y1 = depths_A_down[j], 
                      y2 = depths_A_down[j+1], 
                      color=colors_A[i], edgecolor='none', alpha=0.4)
    
    #Plot x-uncertainty: standard deviations pH
    ax1.fill_betweenx(y = A_set[i]['Depth corrected'], 
                      x1 = A_set[i]['Average pH']+(A_set[i]['Total std dev']), 
                      x2 = A_set[i]['Average pH']-(A_set[i]['Total std dev']),
                      color=colors_A[i], edgecolor='none', alpha=0.4)
    #Plot x-uncertainty: standard deviations [H+]
    ax3.fill_betweenx(y = A_set[i]['Depth corrected'], 
                      x1 = A_set[i]['[H+]']+(A_set[i]['H+ std dev']), 
                      x2 = A_set[i]['[H+]']-(A_set[i]['H+ std dev']),
                      color=colors_A[i], edgecolor='none', alpha=0.4)
    

#CUV B
for i in range(0, 5):
    #Put 0 reading at SWI
    B_set[i]['Depth corrected'] = (B_set[i]['Depth (cm)']-1.15)
    
    #Plot average pH
    ax2.plot('Average pH', 'Depth corrected', data = B_set[i], 
              c = colors_B[i], alpha=0.9, lw=1.5)
    #Plot average [H+]
    ax4.plot('[H+]', 'Depth corrected', data = B_set[i], 
             c = colors_B[i], alpha=0.9, lw=1.5)
    
    #Do error shading in little steps to avoid incorrect shading
    B_set[i]['Depth +0.01'] = (B_set[i]['Depth corrected']+0.01)
    B_set[i]['Depth +0.02'] = (B_set[i]['Depth corrected']+0.02)
    B_set[i]['Depth +0.03'] = (B_set[i]['Depth corrected']+0.03)
    B_set[i]['Depth +0.04'] = (B_set[i]['Depth corrected']+0.04)
    B_set[i]['Depth +0.05'] = (B_set[i]['Depth corrected']+0.05)
    B_set[i]['Depth +0.06'] = (B_set[i]['Depth corrected']+0.06)
    B_set[i]['Depth +0.07'] = (B_set[i]['Depth corrected']+0.07)
    B_set[i]['Depth +0.08'] = (B_set[i]['Depth corrected']+0.08)
    B_set[i]['Depth +0.09'] = (B_set[i]['Depth corrected']+0.09)
    B_set[i]['Depth +0.1'] = (B_set[i]['Depth corrected']+0.1)
    
    B_set[i]['Depth -0.01'] = (B_set[i]['Depth corrected']-0.01)
    B_set[i]['Depth -0.02'] = (B_set[i]['Depth corrected']-0.02)
    B_set[i]['Depth -0.03'] = (B_set[i]['Depth corrected']-0.03)
    B_set[i]['Depth -0.04'] = (B_set[i]['Depth corrected']-0.04)
    B_set[i]['Depth -0.05'] = (B_set[i]['Depth corrected']-0.05)
    B_set[i]['Depth -0.06'] = (B_set[i]['Depth corrected']-0.06)
    B_set[i]['Depth -0.07'] = (B_set[i]['Depth corrected']-0.07)
    B_set[i]['Depth -0.08'] = (B_set[i]['Depth corrected']-0.08)
    B_set[i]['Depth -0.09'] = (B_set[i]['Depth corrected']-0.09)
    B_set[i]['Depth -0.1'] = (B_set[i]['Depth corrected']-0.1)
    
    depths_B_up = [B_set[i]['Depth corrected'],
                  B_set[i]['Depth +0.01'], 
                  B_set[i]['Depth +0.02'],
                  B_set[i]['Depth +0.03'],
                  B_set[i]['Depth +0.04'],
                  B_set[i]['Depth +0.05'],
                  B_set[i]['Depth +0.06'],
                  B_set[i]['Depth +0.07'],
                  B_set[i]['Depth +0.08'],
                  B_set[i]['Depth +0.09'],
                  B_set[i]['Depth +0.1']]
    
    for j in range(len(depths_B_up)-1):
        ax2.fill_between(x = B_set[i]['Average pH'],
                      y1 = depths_B_up[j], 
                      y2 = depths_B_up[j+1], 
                      color=colors_A[i], edgecolor='none', alpha=0.4)
        ax4.fill_between(x = B_set[i]['[H+]'],
                      y1 = depths_B_up[j], 
                      y2 = depths_B_up[j+1], 
                      color=colors_A[i], edgecolor='none', alpha=0.4)
        
    depths_B_down = [B_set[i]['Depth corrected'],
                  B_set[i]['Depth -0.01'], 
                  B_set[i]['Depth -0.02'],
                  B_set[i]['Depth -0.03'],
                  B_set[i]['Depth -0.04'],
                  B_set[i]['Depth -0.05'],
                  B_set[i]['Depth -0.06'],
                  B_set[i]['Depth -0.07'],
                  B_set[i]['Depth -0.08'],
                  B_set[i]['Depth -0.09'],
                  B_set[i]['Depth -0.1']]
    
    for j in range(len(depths_B_down)-1):
        ax2.fill_between(x = B_set[i]['Average pH'],
                      y1 = depths_B_down[j], 
                      y2 = depths_B_down[j+1], 
                      color=colors_A[i], edgecolor='none', alpha=0.4)
        ax4.fill_between(x = B_set[i]['[H+]'],
                      y1 = depths_B_down[j], 
                      y2 = depths_B_down[j+1], 
                      color=colors_A[i], edgecolor='none', alpha=0.4)
    
    #Plot x-uncertainty: standard deviations pH
    ax2.fill_betweenx(y = B_set[i]['Depth corrected'], 
                      x1 = B_set[i]['Average pH']- B_set[i]['Total std dev'], 
                      x2 = B_set[i]['Average pH']+ B_set[i]['Total std dev'],
                      color=colors_B[i], edgecolor='none', alpha=0.4)
    #Plot x-uncertainty: standard deviations [H+]
    ax4.fill_betweenx(y = B_set[i]['Depth corrected'], 
                      x1 = B_set[i]['[H+]']-(B_set[i]['H+ std dev']), 
                      x2 = B_set[i]['[H+]']+(B_set[i]['H+ std dev']),
                      color=colors_B[i], edgecolor='none', alpha=0.4)

    
#Adding sediment and pteropod shading    
for ax in [ax1, ax2, ax3, ax4]:
    ax.axhspan(0, 2.9, color='xkcd:sand', alpha=0.15, lw=0)
    ax.axhline(y=0, color='black', linewidth=0.5, linestyle='--', alpha=0.4)

for ax in [ax1, ax3]:
    ax.axhline(y=-0.15, color='black', linewidth=0.5, linestyle=':', alpha=0.4)
    ax.axhspan(-0.15, 0, color='grey', alpha=0.3, lw=0)

#Add baselines pH
for ax in [ax1, ax2]:
    ax.vlines(x=7.404, ymin=3.25, ymax=-1.2, color='grey', linestyle='--', alpha=0.8)

#Add baselines [H+]
for ax in [ax3, ax4]:
    ax.vlines(x=(10**9)*(10**(-7.404)), ymin=3.25, ymax=-1.2, color='grey', linestyle='--', alpha=0.8)
    
#Obtain handles and labels from ax1
handles, labels = ax1.get_legend_handles_labels()

#Handles is a list, so append pteropod, SWI and calcite sand patches
patch0 = mpatches.Patch(color='white', alpha=0.3, edgecolor='white')
patch1 = mpatches.Patch(color='grey', label='Pteropods', alpha=0.3)
patch2 = mpatches.Patch(color='xkcd:sand', label='Calcite sand', alpha=0.15)
line1 =  Line2D([0], [0], color='black', alpha=0.4, linewidth=0.5, linestyle='--', label='SWI')
line2 =  Line2D([0], [0], color='grey', alpha=0.8, linewidth=1.5, linestyle='--', label='Baseline')
line3 =  Line2D([0], [0], color='grey', alpha=0.8, linewidth=1.5, linestyle='--', label='Baseline')

for i in [line2, patch0, patch1, patch2, line1]:
    handles.append(i)

#Plot legend
# leg = fig.legend(handles=handles, 
#             bbox_to_anchor=(1.2, 0.65), 
#             fontsize='x-small', 
#             ncol=1,
#             title='Hours elapsed:',
#             title_fontsize='x-small')

# leg._legend_box.align = "left"

#%% Adding baselines for analysis
# ax1.axvline(x=A_set[0]['Average pH'][0:20].mean(), color=colors_A[0], linewidth=1, linestyle='--', alpha=0.6)
# ax2.axvline(x=B_set[0]['Average pH'][0:20].mean(), color=colors_B[0], linewidth=1, linestyle='--', alpha=0.6)

#ALL baselines
# A_baselines=[]
# B_baselines=[]
# for i in range(0, 6):
#     A_baselines.append(A_set[i]['Average pH'][0:20].mean())
#     ax1.axvline(x=A_baselines[i], color=colors_A[i], linewidth=1, linestyle='-', alpha=0.8)
# for i in range(0, 5):
#     B_baselines.append(B_set[i]['Average pH'][0:20].mean()) 
#     ax2.axvline(x=B_baselines[i], color=colors_B[i], linewidth=1, linestyle='-', alpha=0.8)

#%%
#Axes
ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
ax1.grid(alpha=0.3, which='both')

#pH
for ax in [ax1, ax2]:
    ax.set_xlim(7.3, 8.1)
    ax.set_ylim(2.9, -1.15, )
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.25))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax.grid(alpha=0.3, which='both')
    ax.set_xlabel(r"pH$_{T}$")

#[H+]
for ax in [ax3, ax4]:
    ax.set_xlim((10**9)*(10**(-8.3)), (10**9)*(10**(-7.35)))
    ax.set_ylim(2.9, -1.15)
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax.invert_xaxis()
    ax.grid(alpha=0.3, which='both')
    ax.set_xlabel("[H$^+$] (nM)")

#Labels
fig.suptitle("pH microprofiles", fontsize=16)
ax1.title.set_text('CUV P')
ax2.title.set_text('CUV CTRL')

ax1.set_ylabel('Depth (cm)')
ax2.axes.yaxis.set_ticklabels([])
ax3.set_ylabel('Depth (cm)')
ax4.axes.yaxis.set_ticklabels([])

plt.tight_layout()
plt.savefig("figures/CUV_profiles_4.png", bbox_inches='tight')
plt.show()
