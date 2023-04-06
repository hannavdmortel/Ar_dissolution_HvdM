
CO2SYS_data = CO2SYS(dtalkwf/rho_sw*1e6,dtCO2wf/rho_sw*1e6,1,2,S,T,T,P,P,Si_molPerKg(1,1).*1e6,PO4_molPerKg(1,1).*1e6,1,10,1);
dCO3wf = CO2SYS_data(:,22)' .* 10^-6.* rho_sw;          %[mol/kg] CO3 
OmegaCwf = dCawf.*dCO3wf./ (kspc.* rho_sw.^2); %[no unit] calcite saturation state
pHwf = CO2SYS_data(:,37)';          %[mol/kg] CO3 

%%
clear export
export(1,:)=pHwf(1,[1 2 4 7 13 45 90 130]);
export(2:70,:)=pHf(:,[1 2 4 7 13 45 90 130]);

%%
clear export
export(1,:)=OmegaCwf(1,[1 2 4 7 13 45 90 130]);
export(2:70,:)=OmegaCf(:,[1 2 4 7 13 45 90 130]);

%%
clear export
export(1,:)=dtalkwf(1,[1 2 4 7 13 45 90 130]).*1e6./rho_sw;
export(2:70,:)=dtalkf(:,[1 2 4 7 13 45 90 130]).*1e6./rho_sw;

%%
clear export
export(1,:)=dtCO2wf(1,[1 2 4 7 13 45 90 130]).*1e6./rho_sw;
export(2:70,:)=dtCO2f(:,[1 2 4 7 13 45 90 130]).*1e6./rho_sw;

%%
clear export
export(1,:)=dCawf(1,[1 2 4 7 13 45 90 130]).*1e6./rho_sw;
export(2:70,:)=dCaf(:,[1 2 4 7 13 45 90 130]).*1e6./rho_sw;

%%
clear export
export(1,:)=NaN(1,8);
export(2:70,:)=pcalcitef(:,[1 2 4 7 13 45 90 130])./27000*100;

%%
clear export
export(1,:)=NaN(1,8);
export(2:70,:)=paragonitef(:,[1 2 4 7 13 45 90 130])./29300*100;

%%
clear export
export(1,:)=NaN(1,8);
export(2:70,:)=TotR_paragonitef(:,[1 2 4 7 13 45 90 130])./29300*100;

%%
clear export
export(1,:)=NaN(1,8);
export(2:70,:)=TotR_pcalcitef(:,[1 2 4 7 13 45 90 130])./29300*100;
