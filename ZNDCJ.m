% Shock and Detonation Toolboox
% http://www.galcit.caltech.edu/EDL/public/sdt/SD_Toolbox/ 
%
% Generate plots and output files for a ZND detonation with the shock front
% traveling at the CJ speed.

clear;clc;
display('ZNDCJ')

P1 = 101300; T1 = 293; 
mech = 'h2air_highT.cti';  

for i = [0.4875, 0.5224, 0.5583, 0.595, 0.6327, 0.7109, 0.7933, 0.8803, 0.9721, 1, 1.02, 1.12, 1.3388, 1.5867, 1.87, 2.1969, 2.38, 2.5783, 2.7939, 2.9089, 3.0291, 3.2867, 3.57]
   j=num2str(2*i); %because at stoichiometry h2/o2 = 2 and "i" is phi here.
   k=num2str(i); %phi
   q = strcat('H2:',j,' O2:1 N2:3.76');
   display(strcat('phi= ',k))
   display(q)
   [cj_speed,gas2] = znd_CJ(1, P1, T1, q, mech, 'h2air', 2);
end