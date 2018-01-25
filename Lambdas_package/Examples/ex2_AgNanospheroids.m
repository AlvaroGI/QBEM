%-------------------------------------------------------------------------
% Author: Alvaro Gomez Inesta (UPC, MIT) & Thomas Christensen (MIT).
% January 2017. (Last updated: January 20, 2018)
% Massachusetts Institute of Technology (MIT) & UPC (CFIS).
%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 2: display the charge distribution for the mode with larger
% dipole moment in z direction of spheroids of different aspect ratio.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear all

%% Data common to all spheroids:
diam = 1; % Diameter [nm]
Nelem = 500; % Approximated number of elements
mode = 'dipolar_max';
dir = [1 0 0];

for ab = [4 2 1 0.5 0.1]
    %% Create spheroid
    disp('Computing spheroid...')

    [p, op] = AgSpheroid(diam, ab, Nelem);

    %% Compute Lambdas
    disp('   -Computing Lambdas...')

    [Lambda_0, Lambda_T, Lambda_II, Sig] = Lambdas(p, op, mode, dir);
    
    figure()
    plot(p,Sig)
    colormap(bluewhitered_mod())

end