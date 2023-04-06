#pH and omega

import PyCO2SYS as pyco2
import numpy as np

#Calculate CO2 system with pH and TA
temp = 20.9 #degrees
pco2 = 480 #ppm
sal = 30

#Initial
pH_ini = 7.4
TA_ini = 681.18

#CUVA
pH_CUVA = 7.55 #average water column last run
TA_CUVA = 2770.59

#CUVB
pH_CUVB = 7.59 #average water column last run
TA_CUVB = 1526.86

parms = dict(salinity=sal, temperature=temp, opt_pH_scale=1, opt_k_carbonic=10)

CUVA = pyco2.sys(**parms, par1=pH_CUVA, par2=TA_CUVA,
                    par1_type = 3, par2_type = 1,
                    uncertainty_into=[
                            "saturation_calcite", 
                            "saturation_aragonite", 
                            "pCO2",
                            "k_calcite",
                            "k_aragonite"],
                    uncertainty_from={
                            "par2": 5.4, #umol/kg
                            "par1": 0.05,
                            "temperature": 0.05, #guess
                            **pyco2.uncertainty_OEDG18
                            })
CUVB = pyco2.sys(**parms, par1=pH_CUVB, par2=TA_CUVB,
                    par1_type = 3, par2_type = 1,
                    uncertainty_into=[
                            "saturation_calcite", 
                            "saturation_aragonite", 
                            "pCO2",
                            "k_calcite",
                            "k_aragonite"],
                    uncertainty_from={
                            "par2": 5.4, #umol/kg
                            "par1": 0.02,
                            "temperature": 0.05, #guess
                            **pyco2.uncertainty_OEDG18
                            })
initial = pyco2.sys(**parms, par1=pH_ini, par2=TA_ini,
                    par1_type = 3, par2_type = 1,
                    uncertainty_into=[
                            "saturation_calcite", 
                            "saturation_aragonite", 
                            "pCO2",
                            "k_calcite",
                            "k_aragonite"],
                    uncertainty_from={
                            "par2": 5.4, #umol/kg
                            "par1": 0.02,
                            "temperature": 0.05, #guess
                            **pyco2.uncertainty_OEDG18
                            })

print('1: ', CUVA['saturation_calcite'],
      '\nu:', CUVA['u_saturation_calcite'],
      '\n'+'2:', CUVB['saturation_calcite'],
      '\nu:', CUVB['u_saturation_calcite'],
      '\n'+'3:', initial['saturation_calcite'],
      '\nu:', initial['u_saturation_calcite'])

print('\n\n1: ', CUVA['saturation_aragonite'],
      '\nu:', CUVA['u_saturation_aragonite'],
      '\n'+'2:', CUVB['saturation_aragonite'],
      '\nu:', CUVB['u_saturation_aragonite'],
      '\n'+'3:', initial['saturation_aragonite'],
      '\nu:', initial['u_saturation_aragonite'])

