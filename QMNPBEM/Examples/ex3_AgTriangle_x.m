%-------------------------------------------------------------------------
% Author: Alvaro Gomez Inesta (UPC, MIT) & Thomas Christensen (MIT).
% January 2017. (Last updated: January 20, 2018)
% Massachusetts Institute of Technology (MIT) & UPC (CFIS).
%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 3: compute the Lambdas and the charge distribution for a 
% triangle (mode with larger dipole moment in x direction).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear all

%% Create triangle
disp('Computing triangle...')
edge = 1;
ab = 5;
Nelem = 500;
[p, op] = AgTriangle(edge,ab,Nelem);

%% Compute Lambdas
disp('   -Computing lambdas...')

mode = 'dipolar_max';
dir = [1 0 0];
[Lambda_0, Lambda_T, Lambda_II, Sig] = Lambdas(p, op, mode, dir);

disp(['Lambda_0 = ' num2str(Lambda_0)])
disp(['Lambda_T = ' num2str(Lambda_T)])
disp(['Lambda_II = ' num2str(Lambda_II)])

%% Plot
figure()
plot(p,'EdgeColor','b')
hold on
plot(p,Sig)
colormap(bluewhitered_mod())




