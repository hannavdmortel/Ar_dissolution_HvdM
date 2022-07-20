% SET UP ENVIRONMENTAL CONDITIONS %
%%% beach test for olivine application %%%

clear all
Station= "Beach"; %give a name to the simulation

%% definition of the spatial domain
z_max=34e-3;     %[m] bottom sediment depth, should be a multiple of z_res
z_res=0.5e-3;     %[m] depth resolution
ndepths = 1 + z_max/z_res;     %[no unit] number of depth layers /!\ no change required /!\
depths = linspace(0, z_max, ndepths); %[m] depth /!\ no change required /!\
z_res = ones(size(depths))*z_res; %[m] depth resolution /!\ no change required /!\

%% definition of the temporal domain
stoptime =1;       %[a] total timespan of the simulation
interval=1/512000;          %[a] time resolution (1/60000 is nine minutes, 1/8760 is one hour; 1/365.2 is a day), 1/512 is good when olivine is in the simulation
t_length=stoptime/interval;      %[no unit] number of time layers /!\ no change required /!\

%% bottom-water environmental conditions
T=20.9;      %[C] temperature
SF_depth=20e-3;      %[m] seafloor depth
S=28.3;   %[psu] salinity
P=1;    %[bar] pressure
rho_sw = gsw_rho(S,T,P);    %[kg/m^3] in situ seawater density computed from GSW toolbox
%[OS: need to update this section to have pressure be computed automatically as a function of seafloor depth]

%% bottom-water values of dissolved species
dtalkw=(681)*1e-6*rho_sw; %[mol/m3] 
dtCO2w=(672)*1e-6*rho_sw; %[mol/m3]  
dtPO4w=(0.75)*1e-6*rho_sw; %[mol/m3] 
dCaw=0.02128./40.087.*(S./1.80655)*rho_sw;  %[mol/m3] Ca, computed from salinity using Riley CG(1967)
dtSiw=(5)*1e-6*rho_sw; %[mol/m3] typical surface waters (Cermeno et al., PNAS 2015)

%% depth-dependent porosity 
phiBeta = 33;   %porosity attenuation coefficient
phiInf = 0.4;   %porosity at infinite depth
phi0 = 0.4;    %porosity at interface
phi = (phi0 - phiInf)*exp(-phiBeta*depths) + phiInf;   %porosity profile (porewater bulk fraction) fitted from station7 mooring3 of cruise NBP98-2 by Sayles et al DSR 2001
phiS=1-phi;   %solid volume fraction
tort=(1-2*log(phi)).^0.5;   %tortuosity from Boudreau (1996, GCA)
tort2=tort.^2;   %tortuosity squared
%[OS: porosity is currently low, as per Steve's suggestion, let's check that again] 

%% pteropod layer 
%we need to create a vector (a mask) with zeros where there is no olivine and ones where there is olivine
mask=zeros(1,ndepths);
%mask(1,1:4)=1; %[mol/m3] olivine is in the three first layers, i.e., over the top centimeter
%olivine concentration computed as the total sediment density divided by olivine molar mass

%% diffusive boundary layer
dbl=SF_depth./20;            %[m] thickness at location, random value for now

saved_time_res=16000;

rerun = 0;
