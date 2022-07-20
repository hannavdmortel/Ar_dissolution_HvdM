from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np

#Respiration test in two beakers, one with poison and one without

#%% Data
#O2 in umol/kg, no given standard deviations

depths = [
    1.5,
    0.5,
    -0.5,
    -2.5]

beaker1_ini = [
    262.697423,
    268.5617658,
    257.2519619,
    261.0218965]

beaker2_ini = [
    219.9714972,
    217.4582074,
    215.7826809,
    223.3225502]

beaker1_mid = [
    269.9734191,
    268.2449199,
    257.873925,
    256.1454258]

beaker2_mid = [
    213.7971965,
    213.3650717,
    207.7474494,
    181.3878373]

beaker1_fin = [
    280.7765388,
    280.7765388,
    279.4801644,
    280.344414]

beaker2_fin = [
    234.1070616,
    230.2179384,
    220.7111931,
    197.8085793]

#%%O2 plots
fig, ax = plt.subplots(dpi=300, figsize=(3.5, 3))

ax.scatter(beaker1_ini, depths, color='xkcd:tangerine', s=10)
ax.plot(beaker1_ini, depths, color='xkcd:tangerine', label = 'Day 1 $[O_{2}]$ poisoned')

ax.scatter(beaker1_mid, depths, color='xkcd:tangerine', s=10)
ax.plot(beaker1_mid, depths, color='xkcd:tangerine', linestyle='--', label = 'Day 4 $[O_{2}]$ poisoned')

ax.scatter(beaker1_fin, depths, color='xkcd:tangerine', s=10)
ax.plot(beaker1_fin, depths, color='xkcd:tangerine', label = 'Day 7 $[O_{2}]$ not poisoned')

ax.scatter(beaker2_ini, depths, color='xkcd:topaz', s=10)
ax.plot(beaker2_ini, depths, color='xkcd:topaz', label = 'Day 1 $[O_{2}]$ not poisoned')

ax.scatter(beaker2_mid, depths, color='xkcd:topaz', s=10)
ax.plot(beaker2_mid, depths, color='xkcd:topaz', linestyle=':', label = 'Day 2 $[O_{2}]$ not poisoned')

ax.scatter(beaker2_fin, depths, color='xkcd:topaz', s=10)
ax.plot(beaker2_fin, depths, color='xkcd:topaz', linestyle=':', label = 'Day 3 $[O_{2}]$ not poisoned')

#Uncertainties
#Inferred from other runs and multiple measurements
ax.fill_betweenx(y = np.array(depths), 
                  x1 = np.array(beaker1_ini)+10, 
                  x2 = np.array(beaker1_ini)-10,
                  color='xkcd:tangerine', edgecolor='none', 
                  alpha=0.2)
ax.fill_betweenx(y = np.array(depths), 
                  x1 = np.array(beaker1_mid)+10, 
                  x2 = np.array(beaker1_mid)-10,
                  color='xkcd:tangerine', edgecolor='none', 
                  alpha=0.2)
ax.fill_betweenx(y = np.array(depths), 
                  x1 = np.array(beaker1_fin)+10, 
                  x2 = np.array(beaker1_fin)-10,
                  color='xkcd:tangerine', edgecolor='none', 
                  alpha=0.2)
ax.fill_betweenx(y = np.array(depths), 
                  x1 = np.array(beaker2_ini)+10, 
                  x2 = np.array(beaker2_ini)-10,
                  color='xkcd:topaz', edgecolor='none', 
                  alpha=0.2)
ax.fill_betweenx(y = np.array(depths), 
                  x1 = np.array(beaker2_mid)+10, 
                  x2 = np.array(beaker2_mid)-10,
                  color='xkcd:topaz', edgecolor='none', 
                  alpha=0.2)
ax.fill_betweenx(y = np.array(depths), 
                  x1 = np.array(beaker2_fin)+10, 
                  x2 = np.array(beaker2_fin)-10,
                  color='xkcd:topaz', edgecolor='none', 
                  alpha=0.2)

#Labels and gridlines
ax.set_xlabel('$[O_{2}]$ (μmol/kg)')
ax.set_ylabel('Depth (cm)')

ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.grid(alpha=0.3, which='both')

#SWI
ax.axhspan(-0.2, 0.2, color='grey', alpha=0.3, lw=0, label='SWI')

#Plot legend
leg = fig.legend(bbox_to_anchor=(1.5, 0.65), 
            fontsize='small')
#Save fig
plt.savefig("figures/resp_test.png", bbox_inches='tight')

#%% Plot delta O2

dO2_beaker1 = np.array(beaker1_mid)-np.array(beaker1_ini)
dO2_beaker2 = np.array(beaker2_mid)-np.array(beaker2_ini)

dO2_beaker1b = np.array(beaker1_fin)-np.array(beaker1_ini)
dO2_beaker2b = np.array(beaker2_fin)-np.array(beaker2_ini)

fig2, ax = plt.subplots(dpi=300, figsize=(3.5, 4))

ax.scatter(dO2_beaker1, depths, color='xkcd:tangerine', s=10)
ax.plot(dO2_beaker1, depths, color='xkcd:tangerine', label = 'Poisoned Day 4 vs Day 1')

ax.scatter(dO2_beaker2, depths, color='xkcd:topaz', s=10)
ax.plot(dO2_beaker2, depths, color='xkcd:topaz', label = 'Not poisoned Day 4 vs Day 1')

ax.scatter(dO2_beaker1b, depths, color='xkcd:tangerine', s=10)
ax.plot(dO2_beaker1b, depths, color='xkcd:tangerine', linestyle='--', label = 'Poisoned Day 7 vs Day 1')

ax.scatter(dO2_beaker2b, depths, color='xkcd:topaz', s=10)
ax.plot(dO2_beaker2b, depths, color='xkcd:topaz', linestyle='--', label = 'Not poisoned Day 7 vs Day 1')


#Uncertainties
ax.fill_betweenx(y = np.array(depths), 
                  x1 = dO2_beaker1+10, 
                  x2 = dO2_beaker1-10,
                  color='xkcd:tangerine', edgecolor='none', 
                  alpha=0.2)
ax.fill_betweenx(y = np.array(depths), 
                  x1 = dO2_beaker2+10, 
                  x2 = dO2_beaker2-10,
                  color='xkcd:topaz', edgecolor='none', 
                  alpha=0.2)
ax.fill_betweenx(y = np.array(depths), 
                  x1 = dO2_beaker1b+10, 
                  x2 = dO2_beaker1b-10,
                  color='xkcd:tangerine', edgecolor='none', 
                  alpha=0.2)
ax.fill_betweenx(y = np.array(depths), 
                  x1 = dO2_beaker2b+10, 
                  x2 = dO2_beaker2b-10,
                  color='xkcd:topaz', edgecolor='none', 
                  alpha=0.2)

#SWI
ax.axhspan(-0.2, 0.2, color='grey', alpha=0.3, lw=0, label='SWI')

#Labels and gridlines
ax.set_xlabel('Δ$[O_{2}]$ (μmol/kg)')
ax.set_ylabel('Depth (cm)')

ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.grid(alpha=0.3, which='both')

#Plot legend
leg = fig2.legend(bbox_to_anchor=(0.8, 0), 
            fontsize='small')

plt.savefig("figures/resp_test_dO2.png", bbox_inches='tight')
