%-------------------------------------------------------------------------
% Author: Alvaro Gomez Inesta (UPC, MIT) & Thomas Christensen (MIT).
% December 2017. (Last updated: March 2, 2018)
% Massachusetts Institute of Technology (MIT) & UPC (CFIS).
%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shortcut to create Ag twin rods (same diameter) embedded in air 
% using MNPBEM for non-retarded static simulations and 'curv' interpolation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% INPUTS:
%   diam    : diameter of the spheres [nm].
%   ab      : height/diameter.
%   gap     : gap between rods [nm].
%   Nelem   : approximated number of elements PER SPHERE.

% OUTPUTS:
%   p       : composite particle (see comparticle in MNPBEM).
%   op      : options (see MNPBEM).


function [p, op] = AgCoupledRods(diam,ab,gap,Nelem)

    op = bemoptions( 'sim', 'stat', 'waitbar', 0, 'interp', 'curv');

    % Dielectric functions
    eps_d = 1;
    epstab = { epsconst( eps_d ), epsdrude( 'Ag' ), epsdrude( 'Ag' ) };

    % Rod
    height = ab*diam; % [nm]
        % Empirical formulas for the rod mesh:
        n1 = (Nelem/2)^.5;
        n2 = (1/(2+ab))*(Nelem/n1);
        n3 = n2*ab/1.5;

    rod = trirod( diam, height, [ n1, n2, n3 ], 'triangles' );

    rod1 = shift( rod, [(gap+diam)/2, 0, 0]);  
    rod2 = shift( rod, [-(gap+diam)/2, 0, 0]);  

    % System of particles
    p = comparticle( epstab, { rod1, rod2 }, [ 2, 1; 3, 1 ], 1, 2, op );
    p = closed(p,[1, 2]);
    
    disp(['Number of elements per particle: ' num2str(length(p.p{1}.faces))]);

end