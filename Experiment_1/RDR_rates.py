import PyCO2SYS as pyco2
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator
import matplotlib.ticker as ticker

#RDR experiment dissolution rate vs degree of undersaturation

#Sediment
SA              = 17.67*10**-3          #area sediment disk [m2]
m_sed           = 1.210                 #mass wet sediment in disk reactor [kg]
porosity_sed    = 0.38                  #[%]
v_sed           = 760                   #[mL] 
v_sw_sed        = v_sed*porosity_sed    #seawater in sediment

#Seawater
V_sw    = 1480.4        #Volume in disk reactor [mL]
rho_sw  = 1.02*10**-3   #Density [kg/mL]
m_sw    = V_sw*rho_sw   #Mass [kg]

#Residence time seawater
r_inflow = 11.37    #Inflow rate [mL/hr]
RT = V_sw/r_inflow  #Residence time [hrs]

#pH
pH_initial = 7.404662

#Top of water column
# pH_RDR1 = 6.95
# pH_RDR2 = 7.42
# pH_RDR3 = 7.86

#Average water column
pH_RDR1 = 7.188251974361844
pH_RDR2 = 7.536897532785713
pH_RDR3 = 7.580079757175361

#TA
TA_initial  = 681  #from QuAAtro [umol/kg]
#TA_initial    = 886    #from Vindta [umol/kg]

TA_RDR1 = 973 #QuAAtro [umol/kg]
TA_RDR2 = 942 #QuAAtro [umol/kg]
TA_RDR3 = 967 #QuAAtro [umol/kg]

#Delta TA
dTA_0        = TA_RDR3-TA_initial #TOTAL [umol/kg]
dTA_1        = TA_RDR1-TA_initial #[umol/kg]
dTA_2        = TA_RDR2-TA_initial   #[umol/kg]
dTA_3        = TA_RDR3-TA_initial   #[umol/kg]

#Delta[Ca2+]
dCa_0        = dTA_0/2      #[umol/kg]  
dCa_1        = dTA_1/2      #[umol/kg]
dCa_2        = dTA_2/2      #[umol/kg]
dCa_3        = dTA_3/2      #[umol/kg]  

#Dissolution rate
r0_dissolution = (m_sw*dCa_0)/(SA*RT)   #[umol/(m2 hr)]
r1_dissolution = (m_sw*dCa_1)/(SA*RT)   #[umol/(m2 hr)]
r2_dissolution = (m_sw*dCa_2)/(SA*RT)   #[umol/(m2 hr)]
r3_dissolution = (m_sw*dCa_3)/(SA*RT)   #[umol/(m2 hr)]

#%% PyCO2SYS

parms = dict(salinity=28.3, 
             temperature=20.9, 
             opt_k_carbonic=10,
             par1_type = 1, 
             par2_type = 3,)

#Initial
RDR0_ = pyco2.sys(**parms, 
    par1=TA_initial, 
    par2=pH_initial,
    uncertainty_into=[
        "saturation_calcite", 
        "saturation_aragonite", 
        "pCO2"],
    uncertainty_from={
        "par1": 5, #umol/kg
        "par2": 0.005, 
        "temperature": 0.05, #guess
        **pyco2.uncertainty_OEDG18
    }    
    ) 

#RDR1
RDR_one = pyco2.sys(**parms, 
    par1=TA_RDR1, 
    par2=pH_RDR1,
    uncertainty_into=[
        "saturation_calcite", 
        "saturation_aragonite", 
        "pCO2"],
    uncertainty_from={
        "par1": 5.4, #umol/kg
        "par2": 0.0147, #sqrt(0.006 (from tris)**2 + average error water column**2)
        "temperature": 0.05, #guess
        **pyco2.uncertainty_OEDG18
    }    
    ) 

#RDR2
RDR_two = pyco2.sys(**parms, 
    par1=TA_RDR2, 
    par2=pH_RDR2,
    uncertainty_into=[
        "saturation_calcite", 
        "saturation_aragonite", 
        "pCO2"],
    uncertainty_from={
        "par1": 5.4, #umol/kg
        "par2": 0.0150, #sqrt(0.006 (from tris)**2 + average error water column**2)
        "temperature": 0.05, #guess
        **pyco2.uncertainty_OEDG18
    }    
    )

#RDR3
RDR_three = pyco2.sys(**parms, 
    par1=TA_RDR3, 
    par2=pH_RDR3,
    uncertainty_into=[
        "saturation_calcite", 
        "saturation_aragonite", 
        "pCO2"],
    uncertainty_from={
        "par1": 5.4, #umol/kg
        "par2": 0.0144, #sqrt(0.006 (from tris)**2 + average error water column**2)
        "temperature": 0.05, #guess
        **pyco2.uncertainty_OEDG18
    }    
    )


#%% Plot saturation state against dissolution rate
fig, ax1 = plt.subplots(dpi=300, figsize=(4.3, 3))

rates = [
    r1_dissolution, 
    r2_dissolution, 
    r3_dissolution]
sat_states = [
    1-RDR_one['saturation_calcite'],
    1-RDR_two['saturation_calcite'],
    1-RDR_three['saturation_calcite']]

sat_states2 = [
    1-RDR_one['saturation_aragonite'],
    1-RDR_two['saturation_aragonite'],
    1-RDR_three['saturation_aragonite']]

colors = ['xkcd:watermelon', 'xkcd:tangerine', 'xkcd:topaz']

y_err = [np.sqrt((0.05/m_sw)**2+(5.4/dTA_1)**2+2*((0.1*10**-3)/SA)**2+(8/RT)**2)*r1_dissolution,
         np.sqrt((0.05/m_sw)**2+(5.4/dTA_2)**2+2*((0.1*10**-3)/SA)**2+(8/RT)**2)*r2_dissolution,
         np.sqrt((0.05/m_sw)**2+(5.4/dTA_3)**2+2*((0.1*10**-3)/SA)**2+(8/RT)**2)*r3_dissolution]

x_err = [RDR_one['u_saturation_calcite'], 
         RDR_two['u_saturation_calcite'], 
         RDR_three['u_saturation_calcite']]

x_err2 = [RDR_one['u_saturation_aragonite'], 
         RDR_two['u_saturation_aragonite'], 
         RDR_three['u_saturation_aragonite']]

#Plot scatter points
ax2 = ax1.twiny()
ax2.scatter(x=sat_states2, y=rates, c=colors, zorder=2)

#Add error bars
# ax1.errorbar(x=sat_states, y=rates, yerr=y_err, ecolor=colors, c='none', zorder=1)
# ax1.errorbar(x=sat_states, y=rates, xerr=x_err, ecolor=colors,  c='none', zorder=1)

ax2.errorbar(x=sat_states2, y=rates, yerr=y_err, ecolor=colors, c='none', zorder=1)
ax2.errorbar(x=sat_states2, y=rates, xerr=x_err2, ecolor=colors, c='none', zorder=1)

#Fill between
# ytop = [(rates[0]+y_err[0]), (rates[1]+y_err[1]), (rates[2]+y_err[2])]
# ybottom = [(rates[0]-y_err[0]), (rates[1]-y_err[1]), (rates[2]-y_err[2])]
# ax1.fill_between(x=sat_states, y1=ytop, y2=ybottom, color='grey', edgecolor='none', alpha=0.15)

ax1.grid(alpha=0.3, which='both', zorder=0)

#Add calcite dissolution rate law from Naviaux 2019
sat_ca = np.linspace(0.25, 0.8, 20)
ax1.plot(sat_ca, (10**(-10.8)*sat_ca**4.7*3.6e13), 
         color='xkcd:grey', label = 'Calcite dissolution law (Naviaux, 2019)')
#Uncertainty shading:
# ax1.fill_between(x = sat_ca, 
#                  y1= (10**(-10.8-0.4)*sat_ca**(4.7+0.7)*3.6e13),
#                  y2= (10**(-10.8+0.4)*sat_ca**(4.7-0.7)*3.6e13),
#                  color='grey', alpha=0.3, edgecolor='none')

#Add aragonite dissolution rate law from 
sat_ar = np.linspace(0, 0.83, 50)
ax2.plot(sat_ar, (0.013*(sat_ar**1.37)*51403), 
         color='xkcd:grey', label = 'Aragonite dissolution law (Dong, 2019)', linestyle='--')
#Uncertainty shading (only for rate order):
# ax2.fill_between(x = sat_ar, 
#                   y1= (0.013*(sat_ar**(1.37+0.18))*51403),
#                   y2= (0.013*(sat_ar**(1.37-0.18))*51403),
#                   color='red', alpha=0.3, edgecolor='none')

#Labels
labels = ['RDR1', 'RDR2', 'RDR3']
ax1.text(sat_states[0]+0.085, 50, labels[0], fontsize=7.8, ha='right')
ax1.text(sat_states[1]+0.019, 50, labels[1], fontsize=7.8, ha='right')
ax1.text(sat_states[2]+0, 50, labels[2], fontsize=7.8, ha='right')

ax1.set_xlabel(r'$(1-Ω_{ca})$')
ax2.set_xlabel(r'$(1-Ω_{ar})$')
ax1.set_ylabel('Dissolution rate ($μmol$ $m^{-2}$ $hr^{-1}$)', va='center', rotation='vertical')

ax1.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
ax1.xaxis.set_minor_locator(AutoMinorLocator(10))

ax2.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
ax2.xaxis.set_minor_locator(AutoMinorLocator(10))

ax1.yaxis.set_minor_locator(AutoMinorLocator(4))

#Legend
leg = fig.legend( 
            bbox_to_anchor=(0.6, 0.825), 
            fontsize='xx-small', 
            ncol=1)

ax1.set_ylim(0, 500)
ax1.set_xlim(0.225, 0.8)
ax2.set_xlim(0.53, 0.83)
plt.tight_layout()
plt.savefig("Figures/rate_vs_omega.png")

#%%
#Time between measurements
dt_RDR1 = 455.5 #hrs
dt_RDR2 = 647.75 #hrs
dt_RDR3 = 380.42 #hrs


