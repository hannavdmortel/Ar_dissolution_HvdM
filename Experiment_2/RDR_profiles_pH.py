import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
import numpy as np

filepath = "C:/Users/hanna/Documents/GitHub/Ar_dissolution_HvdM/Experiment_2"

#%%Import and plot Unisense RDR microprofile data

#For figure
fig, ax = plt.subplots(dpi=300, figsize=(3.5, 5))
colors = ['xkcd:watermelon', 
           'xkcd:tangerine',
           'xkcd:topaz']

# colors1 = ['#fa96a0', '#fd4659', '#fa1633'] #watermelon
# colors2 = ['#fcbc68', '#ff9408', '#ed8702'] #tangerine
# colors3 = ['#bff2ef', '#13bbaf', '#049187'] #topaz

colors_scat = ['#c5c7c5', '#929591', '#606160']

labels = ['Center', 
          '3.5 cm from center',
          '7 cm from center']

#Initial run
run1 = {}
profile_numbers = [1, 2, 3]
for i in profile_numbers:
    run1["{}".format(i)] = pd.read_excel(filepath + "/Data/RDR_profiles_2022-02-03.xlsx", 
                       header=0,
                       sheet_name="Data (Profile experiment {})".format(i),
                       )
    #Change units of depth and rename
    run1["{}".format(i)]["Depth (um)"] = (run1["{}".format(i)]["Depth (um)"]/10000)
    run1["{}".format(i)]["Depth (um)"] = run1["{}".format(i)].rename(columns={"Depth (um)":"Depth (cm)"}, inplace=True)
    L = (run1["{}".format(i)]["Depth (cm)"] < 8.5)
    run1["{}".format(i)]["Depth (cm)"]=run1["{}".format(i)]["Depth (cm)"][L]
    
    #Omit standard dev > 0.2
    L2 = (run1["{}".format(i)]["Std. dev,  (pH unit)Sensor 3 - pH"] < 0.1)
    run1["{}".format(i)] = run1["{}".format(i)][L2]
    
    #Plot data
    run1["{}".format(i)].plot.scatter(x="pH calibrated with TRIS", y="Depth (cm)", ax=ax, 
                   alpha=0.5, s=3, edgecolor='none', color= colors_scat[i-1])

#Plot average
df1 = pd.DataFrame(data = [
    run1['1']['pH calibrated with TRIS'], 
    run1['2']['pH calibrated with TRIS'], 
    run1['3']['pH calibrated with TRIS']])
df1 = pd.DataFrame.transpose(df1)
avg1 = df1.mean(axis='columns')
run1mean = plt.plot(avg1, run1['1']["Depth (cm)"], color=colors[0], linewidth=1.5, label='Mean run 1')

#Second run
run2 = {}
profile_numbers = [1, 2, 3]
for i in profile_numbers:
    run2["{}".format(i)] = pd.read_excel(filepath + "/Data/RDR_profiles_2022-03-02.xlsx", 
                       header=0,
                       sheet_name="Data (Profile experiment {})".format(i),
                       )
    #Change units depth and rename
    run2["{}".format(i)]["Depth (um)"] = run2["{}".format(i)]["Depth (um)"]/10000
    run2["{}".format(i)]["Depth (um)"] = run2["{}".format(i)].rename(columns={"Depth (um)":"Depth (cm)"}, inplace=True)
    
    #Omit standard dev > 0.2
    L = (run2["{}".format(i)]["Std. dev,  (pH unit)Sensor 3 - pH"] < 0.1)
    run2["{}".format(i)] = run2["{}".format(i)][L]
    
    #Plot data
    run2["{}".format(i)].plot.scatter(x="pH calibrated with TRIS", y="Depth (cm)", ax=ax, 
                   alpha= 0.5, s=3, edgecolor='none', color= colors_scat[i-1])

#Plot average
df2 = pd.DataFrame(data = [
    run2['1']['pH calibrated with TRIS'], 
    run2['2']['pH calibrated with TRIS'], 
    run2['3']['pH calibrated with TRIS']])
df2 = pd.DataFrame.transpose(df2)
avg2 = df2.mean(axis='columns')
run2mean = plt.plot(avg2, run2['3']["Depth (cm)"], color=colors[1], linewidth=1.5, label='Mean run 2')

#Final run
run3 = {}
profile_numbers = [1, 2, 3]
for i in profile_numbers:
    run3["{}".format(i)] = pd.read_excel(filepath + "/Data/RDR_profiles_2022-03-18.xlsx", 
                       header=0,
                       sheet_name="Data (Profile experiment {})".format(i),
                       )
    #Change units depth and rename
    run3["{}".format(i)]["Depth (um)"] = run3["{}".format(i)]["Depth (um)"]/10000
    run3["{}".format(i)]["Depth (um)"] = run3["{}".format(i)].rename(columns={"Depth (um)":"Depth (cm)"}, inplace=True)
    run3["{}".format(i)]['Depth (cm)'] = run3["{}".format(i)]['Depth (cm)']+0.35
    #Omit standard dev > 0.2
    L = (run3["{}".format(i)]["Std. dev,  (pH unit)Sensor 3 - pH"] < 0.1)
    run3["{}".format(i)] = run3["{}".format(i)][L]
    
#Fix depths for profiles 2 and 3
run3["2"]['Depth (cm)'] = run3["2"]['Depth (cm)']-0.31
run3["3"]['Depth (cm)'] = run3["3"]['Depth (cm)']-0.47

    #Plot data
for i in profile_numbers:
    run3["{}".format(i)].plot.scatter(x="pH calibrated with TRIS", y="Depth (cm)", ax=ax, 
                   alpha=0.5, s=3, edgecolor='none', color= colors_scat[i-1])
    
#Plot average
df3 = pd.DataFrame(data = [
    run3['1']['pH calibrated with TRIS'], 
    run3['2']['pH calibrated with TRIS'], 
    run3['3']['pH calibrated with TRIS']])
df3 = pd.DataFrame.transpose(df3)
avg3 = df3.mean(axis='columns')
run3mean = plt.plot(avg3, run3['2']["Depth (cm)"], color=colors[2], linewidth=1.5, label='Mean run 3')

#Add baseline
ax.vlines(x=7.404, ymin=0, ymax=8.5, color='grey', linestyle='--', alpha=0.8)

#%% Uncertainty propagation and shading
df1['error'] = np.sqrt(((run1['1']['Std. dev,  (pH unit)Sensor 3 - pH'])**2) + ((run1['3']['Std. dev,  (pH unit)Sensor 3 - pH'])**2) + ((run1['1']['Std. dev,  (pH unit)Sensor 3 - pH'])**2) + (3*0.006**2))
df2['error'] = np.sqrt(((run2['1']['Std. dev,  (pH unit)Sensor 3 - pH'])**2) + ((run2['3']['Std. dev,  (pH unit)Sensor 3 - pH'])**2) + ((run2['1']['Std. dev,  (pH unit)Sensor 3 - pH'])**2) + (3*0.006**2))
df3['error'] = np.sqrt(((run3['1']['Std. dev,  (pH unit)Sensor 3 - pH'])**2) + ((run3['3']['Std. dev,  (pH unit)Sensor 3 - pH'])**2) + ((run3['1']['Std. dev,  (pH unit)Sensor 3 - pH'])**2) + (3*0.006**2))

#Plot x-uncertainty: standard deviations pH
plt.fill_betweenx(y = run1['1']["Depth (cm)"], 
                  x1 = avg1-df1['error'], 
                  x2 = avg1+df1['error'],
                  color=colors[0], edgecolor='none', alpha=0.3)

plt.fill_betweenx(y = run2['3']["Depth (cm)"], 
                  x1 = avg2-df2['error'], 
                  x2 = avg2+df2['error'],
                  color=colors[1], edgecolor='none', alpha=0.3)

plt.fill_betweenx(y = run3['2']["Depth (cm)"], 
                  x1 = avg3-df3['error'], 
                  x2 = avg3+df3['error'],
                  color=colors[2], edgecolor='none', alpha=0.3)

#Plot y-uncertainty: sediment height difference
plt.fill_between(x = avg1,
                  y1 = run1['1']["Depth (cm)"]-0.1, 
                  y2 = run1['1']["Depth (cm)"]+0.1, 
                  color=colors[0], edgecolor='none', alpha=0.3)

plt.fill_between(x = avg2,
                  y1 = run2['3']["Depth (cm)"]-0.1, 
                  y2 = run2['3']["Depth (cm)"]+0.1, 
                  color=colors[1], edgecolor='none', alpha=0.3)

plt.fill_between(x = avg3,
                  y1 = run3['2']["Depth (cm)"]-0.1, 
                  y2 = run3['2']["Depth (cm)"]+0.1, 
                  color=colors[2], edgecolor='none', alpha=0.3)

#%% Formatting

#Labels
# t1 = ax.text(0.09, 0.6, '$Ω_{ca}$ =\n0.22', size=5, transform=ax.transAxes, ha='center')
# t1.set_bbox(dict(facecolor='xkcd:watermelon', alpha=0.3, edgecolor='xkcd:watermelon'))
# t2 = ax.text(0.55, 0.6, '$Ω_{ca}$ =\n0.46', size=5, transform=ax.transAxes, ha='center')
# t2.set_bbox(dict(facecolor='xkcd:tangerine', alpha=0.3, edgecolor='xkcd:tangerine'))
# t3 = ax.text(0.85, 0.6, '$Ω_{ca}$ =\n1.09', size=5, transform=ax.transAxes, ha='center')
# t3.set_bbox(dict(facecolor='xkcd:topaz', alpha=0.3, edgecolor='xkcd:topaz'))

#Axes
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.invert_yaxis()
ax.set_ylim(top=0, bottom=8.5)
#ax.set_xlim(left=7.5, right=8.2)
ax.grid(alpha=0.3, which='both')
#ax.set_xlim(7.3, 8)

#Change y-axis labels to be negative in sediment
a=ax.get_yticks().tolist()
new_yticks = [5, 4, 3, 2, 1, 0, -1, -2, -3, -4]
ax.axes.yaxis.set_ticklabels(new_yticks)

#Labels
#ax.set_title("RDR profiles")
ax.set_xlabel("$pH_{T}$")

#Add sediment shading
ax.axhspan(5, 8.5, color='xkcd:sand', alpha=0.15, lw=0)
ax.axhline(y=5, color='black', linewidth=0.5, linestyle='--', alpha=0.4)

#Legend
labels =    ['Center', 
             '3.5 cm from center',
             '7 cm from center',
             'Calcite sand',
             'Mean run 1',
             'Mean run 2',
             'Mean run 3',
             'SWI']

#legend_dict=dict(zip(labels_patches,cpatches, labels_lines))
patchList = []
for i in [0, 1, 2]:
    patchList.append([mpatches.Patch(facecolor=colors_scat[i], label=labels[i], alpha=0.7)])
patchList.append(mpatches.Patch(facecolor='xkcd:sand', label='Sand', alpha=0.15, lw=0))
for i in [0, 1, 2]:
    patchList.append([Line2D([0], [0], color=colors[i], linewidth=1, label=labels[i+3])])
patchList.append(Line2D([0], [0], color='black', alpha=0.4, linewidth=0.5, linestyle='--', label='SWI'))
   
plt.legend(handles=patchList, 
            labels=labels, 
            ncol=2, 
            fontsize='x-small',
            bbox_to_anchor=(0.95, -0.13),
            handler_map = {list: HandlerTuple(None)})

plt.tight_layout()
plt.savefig("figures/RDR_profiles_pH.png")
plt.show()
