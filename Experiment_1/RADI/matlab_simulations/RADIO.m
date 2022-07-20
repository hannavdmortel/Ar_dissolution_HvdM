%% 1-D Reaction-Advection-Diffusion-Irrigation with Olivine (RADIO) Diagenetic Sediment Module
%% Source code by O. Sulpis, M.P. Humphreys, M. Wilhelmus, D. Carroll & D. Cole
%% Uses: CO2SYS

%Beach

disp("RADIO is running the following experiment:")
disp(Station)
if rerun==1
    disp("it is a rerun: you can stop the run at any time to visualize the evolution of the system")
    disp("years after start of simulation:")
elseif rerun==0
    disp("it is not a rerun: you can stop the run at any time to visualize the evolution of the system")
    disp("when you stop, set rerun to '1' to continue the analysis")
    disp("years after start of simulation")
else
    disp('initial conditions loaded')
    disp("years after start of simulation")
end

tStart = tic;

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% carbonate system initialization this is used only at the first time step to initialize the calc_pco2 program
CO2SYS_data = CO2SYS(dtalkw*1e6/rho_sw,dtCO2w*1e6/rho_sw,1,2,S,T,T,P,P,dtSiw*1e6/rho_sw,dtPO4w*1e6/rho_sw,1,10,1);
k1(1,1:ndepths) = CO2SYS_data(1,67);    %carbonic acid first dissociation constant
k2(1,1:ndepths) = CO2SYS_data(1,68);     %carbonic acid second dissociation constant
k1p(1,1:ndepths) = CO2SYS_data(1,75);     %phosphate constant 1
k2p(1,1:ndepths) = CO2SYS_data(1,76);      %phosphate constant 2
k3p(1,1:ndepths) = CO2SYS_data(1,77);     %phosphate constant 3
kb(1,1:ndepths) = CO2SYS_data(1,72);      %boron constant 
kw(1,1:ndepths) = CO2SYS_data(1,71);       %water dissociation constants
ksi(1,1:ndepths) = CO2SYS_data(1,78);       %silica constants
bt(1,1:ndepths) = CO2SYS_data(1,79);      %[umol/kg] total boron 
omegaC = CO2SYS_data(1,30);          %calcite saturation state
omegaA = CO2SYS_data(1,31);          %aragonite saturation state
co3 = CO2SYS_data(1,22) .* 10^-6;          %[mol/kg] CO3 
hco3 = CO2SYS_data(1,21) .* 10^-6;          %[mol/kg] HCO3
ph = CO2SYS_data(1,37) ;         %pH on the total scale
%pco2 = CO2SYS_data(1,19) ;         %pCO2
Ca_ini = dCaw ./ rho_sw;         %[mol/kg] Ca concentration 
fg(1,1:ndepths)=dtalkw./rho_sw-hco3-2*co3;    %sum of all alkalinity species that are not carbon
kspc = (co3 .* Ca_ini) ./ omegaC;        %[mol2/kg2] calcite in situ solubility
kspa = (co3 .* Ca_ini) ./ omegaA;      %[mol2/kg2] aragonite in situ solubility
ff(1,1:ndepths) = 1;      %random parameter needed for calc_pco2
H(1,1:ndepths) = 10^-ph;       %[mol/kg] H concentration first guess
clear co3 hco3 omegaC omegaA Ca_ini CO2SYS_data pco2 
sit = dtSiw./rho_sw;        %[mol/kg] convert silica concentration
bt = bt .* 10^-6;        %[mol/kg] convert boron concentration

%% temperature dependent "free solution" diffusion coefficients
D_dtalk=1*(0.015169+0.000793*T);       %[m2/a] approximted to bicarbonate diffusion coefficient from Hulse et al (2018)
D_dtCO2=1*(0.015169+0.000793*T);       %[m2/a] approximted to bicarbonate diffusion coefficient from Hulse et al (2018)
D_dtPO4=1*(0.011291+0.000559*T);      %[m2/a] phosphate diffusion coefficient from Li and Gregory (1974)
D_dCa=1*(0.0107+0.001677*T);     %[m2/a] calcium diffusion coefficient from Li and Gregory (1974)

%% depth-dependent porosity and diffusion coefficient loss
delta_phi = -phiBeta.*(phi0 - phiInf).*exp(-phiBeta*depths); % depth-dependent porosity loss
delta_phiS = -delta_phi;   % depth-dependent solid fraction gain
delta_tort2 = -2*delta_phi./phi;   % depth-dependent tortuosity gain [not used in Julia]

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% R.A.D.I. main loop %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if rerun == 0 %concentrations set to random values for solids, bottom water values for not-too-sensitive solutes
    dtalk(1,1:ndepths)=dtalkw; %[mol/m3] random value
    dtCO2(1,1:ndepths)=dtCO2w; %[mol/m3] random value
    dtPO4(1,1:ndepths)=dtPO4w;      
    dCa(1,1:ndepths)=dCaw;  
    dtSi(1,1:ndepths)=dtSiw;  
    %pcalcite(1,1:4)=0.259*27000; %[mol/m3]
    pcalcite(1,1:ndepths)=0.973*27000; %[mol/m3]
    paragonite(1,1:ndepths)=0.4725*29300.*mask; %[mol/m3]

    % saving variables every XY time steps 
    i=1;
    idx=1;
    plot_number=0:t_length/stoptime/saved_time_res:t_length;  %good to look at seasonal signal
%    plot_number=0:t_length/stoptime/32000:t_length;  %good to look at tidal signal
    plot_number(1)=1;
    
elseif rerun==1 %if it is a rerun, initial conditions are concentrations from last time step
    dtalk=dtalkf(:,idx-1)';            %[mol/m3]
    dtCO2=dtCO2f(:,idx-1)';      %[mol/m3]
    dtPO4=dtPO4f(:,idx-1)';            %[mol/m3]
    dCa=dCaf(:,idx-1)';            %[mol/m3]
    dtSi=dtSif(:,idx-1)';            %[mol/m3]
    pcalcite=pcalcitef(:,idx-1)';            %[mol/m3]
    paragonite=paragonitef(:,idx-1)';            %[mol/m3]
    
    plot_number=0:t_length/stoptime/saved_time_res:t_length;  %good to look at seasonal signal
%    plot_number=0:t_length/stoptime/32000:t_length;  %good to look at tidal signal
    i=plot_number(idx-1);
    
else
    % initial condition for solutes: bottom-water value
    dtalk=dtalkic';                %[mol/m3]
    dtCO2=dtCO2ic';                %[mol/m3]
    dtPO4=dtPO4ic';            %[mol/m3]
    dCa=dCaic';            %[mol/m3]
    dtSi=dtSiic';             %[mol/m3]
    pcalcite=pcalciteic';%.*mask;
    paragonite=paragoniteic';%.*mask;
    
    % variable saving
    i=1;
    idx=1;
    plot_number=0:t_length/stoptime/saved_time_res:t_length;  %good to look at seasonal signal
%    plot_number=0:t_length/stoptime/32000:t_length;  %good to look at tidal signal
    plot_number(1)=1;        
    
end

%% short-cut transport variables
DFF=(tort2.* delta_phi./ phi - delta_tort2)./ (tort2.^2);
TR=(2*z_res.* (tort.^2)./ dbl);

%% Prepare for timestep calculations

% Set indices for depth-varying reactions
j = 2:ndepths-1;
jp1 = j + 1;
jm1 = j - 1;
z_res2 = z_res.^2;

% Subsample variables that do not change from timestep to timestep
z_res_j = z_res(j);
z_res2_j = z_res2(j);
DFF_j = DFF(j);
tort2_j = tort2(j);
phi_j = phi(j);
phiS_j = phiS(j);

%% Begin timestep loop
% Start with some O2 and OC
if rerun==0

% Preallocate saving arrays
dtalkf = NaN(ndepths, stoptime+1);
dtCO2f = NaN(ndepths, stoptime+1);
dtPO4f = NaN(ndepths, stoptime+1);
dCaf = NaN(ndepths, stoptime+1);
pcalcitef = NaN(ndepths, stoptime+1);
paragonitef = NaN(ndepths, stoptime+1);
end

for i=i:t_length-1

%     disp(i)
    
    %% carbonate system solver    
    DIC_molPerKg = dtCO2 ./ rho_sw;     %convert DIC to [mol/kg]
    TA_molPerKg = dtalk ./ rho_sw;         %convert TA to [mol/kg]
    PO4_molPerKg = dtPO4 ./ rho_sw;   %convert PO4 to [mol/kg]
    Si_molPerKg = dtSi ./ rho_sw;   %convert dSi to [mol/kg]
    hg = H;

    CO2SYS_data = CO2SYS(TA_molPerKg.*1e6,DIC_molPerKg.*1e6,1,2,S,T,T,P,P,Si_molPerKg.*1e6,PO4_molPerKg.*1e6,1,10,1);
    co3 = CO2SYS_data(:,22)' .* 10^-6.* rho_sw;          %[mol/kg] CO3 
    H = 10.^-CO2SYS_data(:,37)' ;         %pH on the total scale
    
    %% CaCO3 reactions
    OmegaC = dCa.*co3./ (kspc.* rho_sw.^2); %[no unit] calcite saturation state
    OmegaA = dCa.*co3./ (kspa.* rho_sw.^2); %[no unit] aragonite saturation state
    
     if isreal(OmegaC(1,1))==0
         i %this is a little trick to see at which time step the carbonate system crashes
        break
     end
     
     % Aragonite dissolution rate from Dong et al. (2019) EPSL
       for jj=1:ndepths
          if OmegaA(1,jj)<=1 && OmegaA(1,jj)>0.835 %this is the omega value for which both laws are equals
          Rd_aragonite(1,jj)=2*2*paragonite(1,jj).*200*7.6e-5*((1-OmegaA(1,jj)).^0.13);   %  0.0083 %/day * [Aragonite] * (1-OmegaA)^0.13
          elseif OmegaA(1,jj)<=0.835 %this is the omega value for which both laws are equals
          Rd_aragonite(1,jj)=2*2*paragonite(1,jj).*200*8.4e-4*((1-OmegaA(1,jj)).^1.46);   %  0.09 %/day * [Aragonite] * (1-OmegaA)^1.46
          else            
          Rd_aragonite(1,jj)=0;    
          end
       end
        
    %Calcite dissolution from Naviaux et al. (2019) Marine Chemistry
    for jj=1:ndepths
          if OmegaC(1,jj)<=1 && OmegaC(1,jj)>0.8275 %this is the omega value for which both laws are equals
          Rd_calcite(1,jj)=(1/40)*pcalcite(1,jj).*400*6.32e-5*((1-OmegaC(1,jj)).^0.11); %     1.7E-4 %/day * [Calcite] * (1-OmegaC)^0.11
          elseif OmegaC(1,jj)<=0.8275 %this is the omega value for which both laws are equals
          Rd_calcite(1,jj)=(1/40)*pcalcite(1,jj).*400*0.2*((1-OmegaC(1,jj)).^4.7);    %     0.54%/day * [Calcite] * (1-OmegaC)^4.7
          else            
          Rd_calcite(1,jj)=0;    
          end
    end
    
    %Calcite precipitation rate from Zuddas and Mucci, GCA (1998)
    %normalized to the same surface area than for dissolution (4m2/g)
    for jj=1:ndepths
         if OmegaC(1,jj)>1
         Rp_calcite(1,jj)=100000*((OmegaC(1,jj)-1).^1.76);  % 1E5 mol/m3/a * (OmegaC - 1)^1.76
         else 
         Rp_calcite(1,jj)=0;
         end
    end
        
    %% Time-dependent bottom water chemistry
    
    %set up initial bottom water values
    if i==1
        if rerun==0 
            dtalkw_ini=dtalkw;
            dtCO2w_ini=dtCO2w; 
            dtPO4w_ini=dtPO4w; 
            dCaw_ini=dCaw;
        end
    end
      
    %compute diffusive fluxes across sediment-water interface 
    F_dtalki=D_dtalk(1)*phi(1)*(dtalk(:,1)-dtalkw)./dbl; %[mol/m2/a]
    F_dtCO2i=D_dtCO2(1)*phi(1)*(dtCO2(:,1)-dtCO2w)./dbl; %[mol/m2/a]
    F_dtPO4i=D_dtPO4(1)*phi(1)*(dtPO4(:,1)-dtPO4w)./dbl; %[mol/m2/a]
    F_dCai=D_dCa(1)*phi(1)*(dCa(:,1)-dCaw)./dbl; %[mol/m2/a]

    %% compute new bottom-water concentrations as a function of diffusive sediment fluxes, residence time and seafloor depth
    dtalkw=dtalkw+interval.*F_dtalki./SF_depth;
    dtCO2w=dtCO2w+interval.*F_dtCO2i./SF_depth;
    dtPO4w=dtPO4w+interval.*F_dtPO4i./SF_depth;
    dCaw=dCaw+interval.*F_dCai./SF_depth;
    
    %% air-sea CO2 exchange
    
    CO2SYS_data = CO2SYS(dtalkw/rho_sw*1e6,dtCO2w/rho_sw*1e6,1,2,S,T,T,P,P,Si_molPerKg(1,1).*1e6,PO4_molPerKg(1,1).*1e6,1,10,1);
    pCO2w=CO2SYS_data(:,19);
    
    K0=2.9e-5; %[mol/m3/uatm] Weiss, Mar Chem 1974
    kCO2=200; %[m/a] Typical value, Dinauer and Mucci BG 2017
    F_pCO2=kCO2*K0*(pCO2w-480); %[mol/m2/a]  Flux of CO2 from the water to the air 
    pCO2w=pCO2w-interval*((F_pCO2/SF_depth)/K0); %[uatm]

    CO2SYS_data = CO2SYS(dtalkw/rho_sw*1e6,pCO2w,1,4,S,T,T,P,P,Si_molPerKg(1,1).*1e6,PO4_molPerKg(1,1).*1e6,1,10,1);
    dtCO2w=(CO2SYS_data(:,6)+CO2SYS_data(:,7)+CO2SYS_data(:,8))./1e6.*rho_sw;    
    
    %% Calculate all reactions (24 species, units: [mol/m3/a])
    TotR_dtalk = + phiS./ phi.*2.* (Rd_calcite + Rd_aragonite - Rp_calcite);
    TotR_dtCO2 = phiS./phi.*(Rd_calcite + Rd_aragonite - Rp_calcite);
    TotR_dCa = phiS./ phi.* (Rd_calcite + Rd_aragonite - Rp_calcite); 
    TotR_pcalcite = -Rd_calcite + Rp_calcite; 
    TotR_paragonite = -Rd_aragonite;
    
    %% top boundary condition: prescribed solid fluxes, diffusive boundary layer control on solutes, constant flux for olivine
    
     dtalk_1 = dtalk(1) + interval * ( D_dtalk(1) / tort2(1) * (2*dtalk(2) - 2*dtalk(1) + TR(1) * (dtalkw - dtalk(1))) / (z_res(1)^2) ... %diffusion
        + TotR_dtalk(1)); %reaction
    
    dtCO2_1 = dtCO2(1) + interval * ( D_dtCO2(1) / tort2(1) * (2*dtCO2(2) - 2*dtCO2(1) + TR(1) * (dtCO2w - dtCO2(1))) / (z_res(1)^2) ... %diffusion
        + TotR_dtCO2(1)); %reaction
    
    dtPO4_1 = dtPO4(1) + interval * ( D_dtPO4(1) / tort2(1) * (2*dtPO4(2) - 2*dtPO4(1) + TR(1) * (dtPO4w - dtPO4(1))) / (z_res(1)^2)); %diffusion

    dCa_1 = dCa(1) + interval * ( D_dCa(1) / tort2(1) * (2*dCa(2) - 2*dCa(1) + TR(1) * (dCaw - dCa(1))) / (z_res(1)^2) ... %diffusion
        + TotR_dCa(1)); %reaction
 
    pcalcite_1 = pcalcite(1) + interval * (TotR_pcalcite(1)); %reaction

    paragonite_1 = paragonite(1) + interval * (TotR_paragonite(1)); %reaction
        
    %% bottom boundary condition: gradients disappear
    %[OS: should try to implement a 'constant flux' boundary condition instead]

    dtalk_z = dtalk(ndepths) + interval * (D_dtalk / tort2(ndepths) * 2 * ((dtalk(ndepths-1) - dtalk(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        +TotR_dtalk(ndepths));

%    dtalk_z = dtalk(ndepths) + interval * (TotR_dtalk(ndepths));
    
    dtCO2_z = dtCO2(ndepths) + interval * (D_dtCO2 / tort2(ndepths) * 2 * ((dtCO2(ndepths-1) - dtCO2(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        +TotR_dtCO2(ndepths));

%    dtCO2_z = dtCO2(ndepths) + interval * (TotR_dtCO2(ndepths));
    
    dtPO4_z = dtPO4(ndepths) + interval * (D_dtPO4 / tort2(ndepths) * 2 * ((dtPO4(ndepths-1) - dtPO4(ndepths)) / z_res(ndepths).^2));  %diffusion
    
%    dtPO4_z = dtPO4(ndepths);
    
    dCa_z = dCa(ndepths) + interval * (D_dCa / tort2(ndepths) * 2 * ((dCa(ndepths-1) - dCa(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        +TotR_dCa(ndepths));

%    dCa_z = dCa(ndepths) + interval * (TotR_dCa(ndepths));

    pcalcite_z = pcalcite(ndepths) + interval * (TotR_pcalcite(ndepths));

    paragonite_z = paragonite(ndepths) + interval * (TotR_paragonite(ndepths));
    
    %% all other depths
      
    % Total alkalinity
    dtalk_j = dtalk(j);
    dtalk_jp1 = dtalk(jp1);
    dtalk_jm1 = dtalk(jm1);
    dtalk(j) = dtalk_j + interval*(TotR_dtalk(j) + ...
        (D_dtalk./tort2_j).*((dtalk_jp1 - 2*dtalk_j + dtalk_jm1)./z_res2_j));

    % Dissolved inorganic carbon
    dtCO2_j = dtCO2(j);
    dtCO2_jp1 = dtCO2(jp1);
    dtCO2_jm1 = dtCO2(jm1);
    dtCO2(j) = dtCO2_j + interval*(TotR_dtCO2(j) + ...
        (D_dtCO2./tort2_j).*((dtCO2_jp1 - 2*dtCO2_j + dtCO2_jm1)./z_res2_j));

    % Phosphate
    dtPO4_j = dtPO4(j);
    dtPO4_jp1 = dtPO4(jp1);
    dtPO4_jm1 = dtPO4(jm1);
    dtPO4(j) = dtPO4_j + interval*( ...
        (D_dtPO4./tort2_j).*((dtPO4_jp1 - 2*dtPO4_j + dtPO4_jm1)./z_res2_j));

    % Dissolved calcium
    dCa_j = dCa(j);
    dCa_jp1 = dCa(jp1);
    dCa_jm1 = dCa(jm1);
    dCa(j) = dCa_j + interval*(TotR_dCa(j) + ...
        (D_dCa./tort2_j).*((dCa_jp1 - 2*dCa_j + dCa_jm1)./z_res2_j));
    
    % Calcite
    pcalcite_j = pcalcite(j);
    pcalcite_jp1 = pcalcite(jp1);
    pcalcite_jm1 = pcalcite(jm1);
    pcalcite(j) = pcalcite_j + interval*(TotR_pcalcite(j));
    
    % Aragonite
    paragonite_j = paragonite(j);
    paragonite_jp1 = paragonite(jp1);
    paragonite_jm1 = paragonite(jm1);
    paragonite(j) = paragonite_j + interval*(TotR_paragonite(j)); 
   
    %% Set top and bottom conditions in arrays
    dtalk(1) = dtalk_1; 
    dtCO2(1) = dtCO2_1; 
    dtPO4(1) = dtPO4_1;
    dCa(1) = dCa_1;
    pcalcite(1) = pcalcite_1;
    paragonite(1) = paragonite_1;
    dtalk(ndepths) = dtalk_z; 
    dtCO2(ndepths) = dtCO2_z; 
    dtPO4(ndepths) = dtPO4_z;
    dCa(ndepths) = dCa_z;
    pcalcite(ndepths) = pcalcite_z;
    paragonite(ndepths) = paragonite_z;
    
    %% set very small or negative concentration to zero
    dtalk(dtalk<0)=0;
    dtCO2(dtCO2<0)=0;
    dtPO4(dtPO4<0)=0;
    dCa(dCa<0)=0;
    pcalcite(pcalcite<0)=0;    
    paragonite(paragonite<0)=0;    

    %% save data 

     if i == plot_number(idx)
       disp(plot_number(idx)*interval)
       dtalkf(:, idx) = dtalk;
       dtCO2f(:, idx) = dtCO2; 
       dtPO4f(:, idx) = dtPO4;
       dCaf(:, idx) = dCa;
       OmegaCf(:,idx)=OmegaC;
       pcalcitef(:, idx) = pcalcite; 
       paragonitef(:, idx) = paragonite; 
       dtalkwf(1,idx)=dtalkw;
       dtCO2wf(1,idx)=dtCO2w;
       dtPO4wf(1,idx)=dtPO4w;
       dCawf(1,idx)=dCaw;
       pHf(:,idx)=-log10(H)';
       pCO2wf(:,idx)=pCO2w;
       TotR_paragonitef(:,idx)=TotR_paragonite;
       TotR_pcalcitef(:,idx)=TotR_pcalcite;
       Rd_calcitef(:,idx)=Rd_calcite;
       Rp_calcitef(:,idx)=Rp_calcite;
       idx=idx+1;
     end  
end

tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));
