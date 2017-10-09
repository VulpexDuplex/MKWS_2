% Shock and Detonation Toolboox
% http://www.galcit.caltech.edu/EDL/public/cantera/html/SD_Toolbox/


clear; clc;
display('Initial Temperature Series')

mech = 'h2air_highT.cti';   
gas = importPhase(mech); 
gas1 = importPhase(mech);
nsp = nSpecies(gas);
T1 = [];

% find Hydrogen nitrogen, and oxygen indices
ih2 = speciesIndex(gas,'H2');
io2  = speciesIndex(gas,'O2');
in2  = speciesIndex(gas,'N2');
x = zeros(nsp, 1);
x(ih2,1) = 2.0;
x(io2,1) = 1.0;
x(in2,1) = 3.76;
To = 300;
P1 = 100000;
% fig_num & fname are for 'znd' - use '0' for no output
fig_num = 0;
fname = 0;
% plots: 1 = make plots, otherwise no plots
plots = 1;
display('Initial Conditions')
display('   Equivalence Ratio = 1, Temperature = 300 K')

npoints=15;
   disp(['For ', num2str(npoints), ' values of P1'])
for i = 1:npoints
   T1(i) = To+50*(i-1);
   disp([' ', num2str(i),  ': T1 = ', num2str(T1(i))])
   set(gas,'Temperature',T1(i),'Pressure',P1,'MoleFractions',x);
   
   %%%Constant Volume Explosion Data%%%
   % FIND POST SHOCK STATE FOR GIVEN SPEED
   [cj_speed(i)] = CJspeed(P1,T1(i),x,mech,0); 
   set(gas,'Temperature',T1(i),'Pressure',P1,'MoleFractions',x);
   [gas] = PostShock_fr(cj_speed(i), P1, T1(i), x, mech);
   % SOLVE CONSTANT VOLUME EXPLOSION ODES
   [CVout] = explosion(gas,fig_num);
   ind_time_CV(i) = CVout.ind_time;
   exo_time_CV(i) = CVout.exo_time;   

   %%%%%ZND Detonation Data%%%%%
   % FIND POST SHOCK STATE FOR GIVEN SPEED
   set(gas1, 'T', T1(i), 'P', P1, 'X', x);
   [gas] = PostShock_fr(cj_speed(i), P1, T1(i), x, mech);
   Ts(i) = temperature(gas); %frozen shock temperature   
   Ps(i) = pressure(gas); %frozen shock pressure
   
   %%Calculate CJstate Properties%%%
   [gas] = PostShock_eq(cj_speed(i),P1,T1(i),x, mech);
   T2(i) = temperature(gas);
   P2(i) = pressure(gas);
   rho2(i) = density(gas);

   %Approximate the effective activation energy using finite differences
    Ta = Ts(i)*(1.02);
    set(gas, 'T', Ta, 'P', Ps(i), 'X', x);
    [CVouT1(i)] = explosion(gas,0);
    Tb = Ts(i)*(0.98);
    set(gas, 'T', Tb, 'P', Ps(i), 'X', x);
    [CVout2] = explosion(gas,0); 
    taua = CVouT1(i).ind_time;
    taub = CVout2.ind_time;
    % Approximate effective activation energy for CV explosion
    if(taua==0 || taub==0)
        theta_effective_CV(i) = 0;
    else
        theta_effective_CV(i) = 1/Ts(i)*((log(taua)-log(taub))/((1/Ta)-(1/Tb))); 
    end
      
    [gas] = PostShock_fr(cj_speed(i)*1.02,P1,T1(i),x, mech);
    Ta= temperature(gas);
    [ZNDouT1(i)] = znd(gas,gas1,fig_num,cj_speed(i)*1.02,fname);
    [gas] = PostShock_fr(cj_speed(i)*0.98,P1,T1(i),x, mech);
    Tb = temperature(gas);
    [ZNDout2] = znd(gas,gas1,fig_num,cj_speed(i)*0.98,fname);
    ind_lena = ZNDouT1(i).ind_len_ZND;
    ind_lenb = ZNDout2.ind_len_ZND;
end

if(plots==1)
    % make plots
    close all;
    fontsize=12;
    figure;
    subplot(2,2,1);
    plot(T1(:),T2(:),'k');
    xlabel('Initial Temperature (K)','FontSize',fontsize);
    ylabel('Temperature (K)','FontSize',fontsize);
    title('Post CJ State Temperature','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    subplot(2,2,2);
    plot(T1(:),P2(:),'k');
    xlabel('Initial Temperature (K)','FontSize',fontsize);
    ylabel('Pressure (Pa)','FontSize',fontsize);
    title('Post CJ State Pressure','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    subplot(2,2,3);
    plot(T1(:),rho2(:),'k');
    xlabel('Initial Temperature (K)','FontSize',fontsize);
    ylabel('Density (kg/m^3)','FontSize',fontsize);
    title('Post CJ State Density','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    subplot(2,2,4);
    plot(T1(:),cj_speed(:),'k');
    xlabel('Initial Temperature (K)','FontSize',fontsize);
    ylabel('Ucj (m/s)','FontSize',fontsize);
    title('CJ Speed','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    %%% Plots for the Induction Zone (times and lengths)
    figure;
    subplot(2,2,1);
    plot(T1(:),ind_time_CV(:),'k');
    xlabel('Initial Temperature (K)','FontSize',fontsize);
    ylabel('tau_{CV_i} (s)','FontSize',fontsize);
    title('Induction time for CV explosion','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

      %%% Plots for the Exothermic Zone (times and lengths)
    figure;
    subplot(2,2,1);
    plot(T1(:),exo_time_CV(:),'k');
    xlabel('Initial Temperature (K)','FontSize',fontsize);
    ylabel('tau_{CV_e} (s)','FontSize',fontsize);
    title('Exothermic time for CV explosion','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

   
    %Ts and Tf for ZND detonation
    figure;
    subplot(1,2,1);
    plot(T1(:),Ts(:),'k');
    xlabel('Initial Temperature (K)','FontSize',fontsize);
    ylabel('T_s (K)','FontSize',fontsize);
    title('Frozen shock Temperature for ZND detonation','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    %Approximation of the effective activation energy for CV explosion
    figure;
    plot(T1, theta_effective_CV);
    xlabel('Initial Temperature (K)','FontSize',fontsize);
    ylabel('Theta_{CV} (J)','FontSize',fontsize);
    title('Effective activation energy for CV explosion','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE OUTPUT TEXT FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn = ['Temperature_Series.plt'];
d = date;
P = P1/oneatm;
	
fid = fopen(fn, 'w');
fprintf(fid, '# ZND DETONATION STRUCTURE CALCULATION AND CONSTANT VOLUME EXPLOSION\n');
fprintf(fid, '# CALCULATION RUN ON %s\n\n', d);
	
fprintf(fid, '# INITIAL CONDITIONS\n');
fprintf(fid, '# TEMPERATURE (K) %4.1f\n', T1(i));
%fprintf(fid, '# PRESSURE (ATM) %2.1f\n', P);
fprintf(fid, '# DENSITY (KG/M^3) %1.4e\n', density(gas1));
q = 'H2:2 O2:1 N2:3.76';    
fprintf(fid, ['# SPECIES MOLE FRACTIONS: ' q '\n']);

%fprintf(fid, '# REACTION ZONE STRUCTURE\n\n');

fprintf(fid, '# THE OUTPUT DATA COLUMNS ARE:\n');
fprintf(fid, 'Variables = "Temperature state 1", "Temperature state 2", "Pressure state 2", "density state 2", "CJ Speed", "Temperature Post Shock", "Pressure Post Shock", "Induction time CV", "Exothermic time CV",  "Effective Activation Energy CV ",  "Effective Activation Energy ZND ", "Induction  time ZND", "Induction length ZND", "Exothermic time ZND", "Exothermic length ZND", "Final Temperature ZND"\n');

z = [T1; T2; P2; rho2; cj_speed; Ts; Ps; ind_time_CV; exo_time_CV; theta_effective_CV];
fprintf(fid, '%14.5e \t %14.5e \t %14.5e \t %14.5e \t %14.5e \t %14.5e \t %14.5e \t %14.5e  \t %14.5e  \t %14.5e  \n',z);

fclose(fid);
