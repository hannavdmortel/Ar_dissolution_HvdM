#pH and omega

import PyCO2SYS as pyco2
import numpy as np

#Calculate CO2 system with pH and TA
temp = 20.9 #degrees
pco2 = 480 #ppm
sal = 28.3

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
                    par1_type = 3, par2_type = 1)
CUVB = pyco2.sys(**parms, par1=pH_CUVB, par2=TA_CUVB,
                    par1_type = 3, par2_type = 1)
initial = pyco2.sys(**parms, par1=pH_ini, par2=TA_ini,
                    par1_type = 3, par2_type = 1)

print('1: ', CUVA['saturation_calcite'],
      '\n'+'2:', CUVB['saturation_calcite'],
      '\n'+'3:', initial['saturation_calcite'])

