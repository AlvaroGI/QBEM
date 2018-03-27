%-------------------------------------------------------------------------
% Author: Alvaro Gomez Inesta (UPC, MIT) & Thomas Christensen (MIT).
% January 2017. (Last updated: January 20, 2018)
% Massachusetts Institute of Technology (MIT) & UPC (CFIS).
%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 1: display the modes for an Ag nanorod and manually select
% one of them to compute the Lambdas. The modes will be ordered in
% descending dipole moment in a direction specified by the user.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear all

%% Create rod
disp('Computing rod...')

diam = 1; % Diameter [nm]
ab = 3;
Nelem = 500; % Approximated number of elements

[p, op] = AgRod(diam, ab, Nelem);

%% Compute Lambdas
disp('   -Computing Lambdas...')

mode = -1;
dir = [0 0 1];
[Lambda_0, Lambda_T, Lambda_II, Sig] = Lambdas(p, op, mode, dir);

disp(['Lambda_0 = ' num2str(Lambda_0)])
disp(['Lambda_T = ' num2str(Lambda_T)])
disp(['Lambda_II = ' num2str(Lambda_II)])