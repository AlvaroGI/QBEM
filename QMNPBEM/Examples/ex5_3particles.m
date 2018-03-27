%-------------------------------------------------------------------------
% Author: Alvaro Gomez Inesta (UPC, MIT) & Thomas Christensen (MIT).
% January 2017. (Last updated: January 23, 2018)
% Massachusetts Institute of Technology (MIT) & UPC (CFIS).
%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 5: select the mode to compute Lambdas in a three particles system
% by plotting all of them.
% Estimated time to run: <30s.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear all

%% Particles
disp('Creating particles...')

op = bemoptions( 'sim', 'stat', 'waitbar', 0, 'interp', 'curv');
    
% Dielectric functions
eps_d = 1;
epstab = { epsconst( eps_d ), epsdrude( 'Ag' ), epsdrude( 'Ag' ) };

% Spheres
diam = 1; % diameter [nm]
ab = 1; % diameter/gap
Nelem = 300;

sphere = trisphere( Nelem, diam );

gap = diam/ab;
sphere1 = shift( sphere, [-(diam+gap)/2, 0, 0]);
sphere2 = shift( sphere, [(diam+gap)/2, 0, 0]);
sphere2 = scale( sphere2, [1,2,3] );

% Rod
ab_rod = 3;
height = ab_rod*diam; % [nm]
    % Empirical formulas for the rod mesh:
    n1 = (Nelem/2)^.5;
    n2 = (1/(2+ab_rod))*(Nelem/n1);
    n3 = n2*ab_rod/1.5;

rod = trirod( diam, height, [ n1, n2, n3 ], 'triangles' );

rod = shift( rod, [0, (gap^2-(gap/2)^2)^.5, 0]);  
    
% System of particles
p = comparticle( epstab, { sphere1, sphere2, rod }, [ 2, 1; 3, 1 ], 1, 2, op );
p = closed(p,[1, 2, 3]);


%% Compute Lambdas
disp('   -Computing lambdas...')

mode = 0;
dirs = [0 0 1;
        0 0 1];
[Lambda_0, Lambda_T, Lambda_II, Sig] = Lambdas(p, op, mode, dirs);

plot(p,Sig)
colormap(bluewhitered_mod())

disp(['Lambda_0 = ' num2str(Lambda_0)])
disp(['Lambda_T = ' num2str(Lambda_T)])
disp(['Lambda_II = ' num2str(Lambda_II)])