%-------------------------------------------------------------------------
% Author: Alvaro Gomez Inesta (UPC, MIT) & Thomas Christensen (MIT).
% January 2017. (Last updated: January 20, 2018)
% Massachusetts Institute of Technology (MIT) & UPC (CFIS).
%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 4: select the mode to compute Lambdas in a twin spheres system
% by plotting all of them in descending product of dipole moments.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear all

%% Coupled spheres
diam = 1; % diameter [nm]
ab = 5; % diameter/gap
Nelem = 500;

disp('Creating coupled spheres...')
[p, op] = AgCoupledSpheres(diam,ab,Nelem);

%% Compute Lambdas
disp('   -Computing lambdas...')

mode = -1;
dirs = [0 0 1;
        0 0 -1];
[Lambda_0, Lambda_T, Lambda_II, Sig] = Lambdas(p, op, mode, dirs);

plot(p,Sig)
colormap(bluewhitered_mod())
disp(['Lambda_0 = ' num2str(Lambda_0)])
disp(['Lambda_T = ' num2str(Lambda_T)])
disp(['Lambda_II = ' num2str(Lambda_II)])