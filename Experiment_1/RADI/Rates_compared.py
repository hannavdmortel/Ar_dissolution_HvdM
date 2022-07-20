import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker 
from matplotlib.ticker import AutoMinorLocator

#Rate expressions from RADIv1 paper (Sulpis et al., 2022) vs our rate expressions

#Calcite and aragonite concentrations
Ca = 10e4 # mol/m3
Ar = 10e4 # mol/m3

#List of saturation states
Ω_ca  = np.linspace(-1, 1, 400)
Ω_ca1 = np.linspace(0.828, 1, 400)
Ω_ca2 = np.linspace(-1, 0.828, 400)

Ω_ar1 = np.linspace(0.835, 1, 400)
Ω_ar2 = np.linspace(-1, 0.835, 400)

##RATE CONSTANTS##
#Dissolution
# 0.828 < Ω_ca < 1 
radi_k_diss_ca1 = 6.3e-3 # /yr
k_diss_ca1 = 1.7e-4/100*365 # /yr
# Ω_ca =< 0.828 
radi_k_diss_ca2 = 20 # /yr
k_diss_ca2 = 0.54/100*365 # /yr

# 0.835 < Ω_ar < 1 
radi_k_diss_ar1 = 3.8e-3 # /yr
k_diss_ar1 = 0.0166/100*365 # /yr
# Ω_ar =< 0.835 
radi_k_diss_ar2 = 4.2e-2 # /yr
k_diss_ar2 = 0.18/100*365 # /yr

#Precipitation
radi_k_prec_ca = 0.4 # mol/m3/a
k_prec_ca = 1.5e5 # mol/m3/a

##RATE EXPRESSIONS##
#Traditional
R_diss_ca = Ca*0.1*365*(1-Ω_ca)**4.5

# 0.828 < Ω_ca < 1 
radi_R_diss_ca1 = Ca*radi_k_diss_ca1*(1-Ω_ca1)**0.11
R_diss_ca1 = Ca*k_diss_ca1*(1-Ω_ca1)**0.11

# Ω_ca =< 0.828 
radi_R_diss_ca2 = Ca*radi_k_diss_ca2*(1-Ω_ca2)**4.7
R_diss_ca2 = Ca*k_diss_ca2*(1-Ω_ca2)**4.7

# 0.835 < Ω_ar < 1 
radi_R_diss_ar1 = Ar*radi_k_diss_ar1*(1-Ω_ar1)**0.13
R_diss_ar1 = Ar*k_diss_ar1*(1-Ω_ar1)**0.13

# Ω_ar =< 0.835 
radi_R_diss_ar2 = Ar*radi_k_diss_ar2*(1-Ω_ar2)**1.46
R_diss_ar2 = Ar*k_diss_ar2*(1-Ω_ar2)**1.46

#Precipitation
radi_R_prec_ca = radi_k_prec_ca*abs(Ω_ca-1)**1.76
R_prec_ca = k_prec_ca*abs(Ω_ca-1)**1.76

##PLOTS##
fig, (ax1, ax2) = plt.subplots(2, 1, dpi=300, figsize=(4, 5))

#Plot rate of CaCO3 dissolution against 1-Ω
ax1.plot(1-Ω_ca, R_diss_ca, c='xkcd:grey', linestyle='--', label='Traditional calcite dissolution \nlaw, e.g. Archer (1991)')

ax1.plot(1-Ω_ca1, radi_R_diss_ca1, c='xkcd:topaz', label='RADI $R_{diss, ca}$')
ax1.plot(1-Ω_ca2, radi_R_diss_ca2, c='xkcd:topaz')

ax1.plot(1-Ω_ca1, R_diss_ca1, c='xkcd:topaz', linestyle=':', label='Our $R_{diss, ca}$')
ax1.plot(1-Ω_ca2, R_diss_ca2, c='xkcd:topaz', linestyle=':')

ax1.plot(1-Ω_ar1, radi_R_diss_ar1, c='xkcd:watermelon', label='RADI $R_{diss, ar}$')
ax1.plot(1-Ω_ar2, radi_R_diss_ar2, c='xkcd:watermelon')

ax1.plot(1-Ω_ar1, R_diss_ar1, c='xkcd:watermelon', linestyle=':', label='Our $R_{diss, ar}$')
ax1.plot(1-Ω_ar2, R_diss_ar2, c='xkcd:watermelon', linestyle=':')

#Plot rate of CaCO3 precipitation against 1-Ω
ax2.plot(Ω_ca-1, radi_R_prec_ca, c='xkcd:lightish green', label='RADI $R_{prec, ca}$')
ax2.plot(Ω_ca-1, R_prec_ca, c='xkcd:lightish green', linestyle=':', label='Our $R_{prec, ca}$')

#Labels and gridlines
ax1.set_xlabel('1-Ω')
ax1.set_ylabel('Rate of $ {CaCO_3}$ dissolution \n($mol$ $m^{-3}$ $a^{-1}$)', fontsize=8)
ax1.set_yscale('log')

ax2.set_ylabel('Rate of $ {CaCO_3}$ precipitation \n($mol$ $m^{-3}$ $a^{-1}$)', fontsize=8)
ax2.set_yscale('log')
ax2.invert_yaxis()

ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
locmin = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10) * .1,
                                      numticks=100)
ax1.yaxis.set_minor_locator(locmin)
ax1.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax2.yaxis.set_minor_locator(locmin)
ax2.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

ax1.set_xlim(-0.05, 0.6)
ax2.set_xlim(-0.6, 0.05)
ax1.set_ylim(1e-5, 1e7)
ax2.set_ylim(1e7, 1e-5)

#Plot legend
leg = fig.legend(bbox_to_anchor=(0.82, -0.01), 
            fontsize='small')

ax1.grid(alpha=0.3, which='major')
ax2.grid(alpha=0.3, which='major')

fig.tight_layout()

plt.savefig("figures/rates_compared.png", bbox_inches='tight')

