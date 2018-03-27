%-------------------------------------------------------------------------
% Author: Alvaro Gomez Inesta (UPC, MIT) & Thomas Christensen (MIT).
% December 2017. (Last updated: January 20, 2018)
% Massachusetts Institute of Technology (MIT) & UPC (CFIS).
%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shortcut to create Ag twin spheres (same diameter) embedded in air 
% using MNPBEM for non-retarded static simulations and 'curv' interpolation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% INPUTS:
%   diameter: diameter of the spheres [nm].
%   ab      : diameter/gap.
%   Nelem   : approximated number of elements PER SPHERE.

% OUTPUTS:
%   p       : composite particle (see comparticle in MNPBEM).
%   op      : options (see MNPBEM).


function [p, op] = AgCoupledSpheres(diameter,ab,Nelem)

    op = bemoptions( 'sim', 'stat', 'waitbar', 0, 'interp', 'curv');
    
    % Dielectric functions
    eps_d = 1;
    epstab = { epsconst( eps_d ), epsdrude( 'Ag' ), epsdrude( 'Ag' ) };
    
    %% Single sphere
    % Nelem to sizefactor:
    Nelem_vec = [508 644 796 964 1054 1246 1348 1564 1678 1796 2044 2446 2884];
    old = Nelem;
    Nelem = interp1(Nelem_vec,Nelem_vec,Nelem,'nearest');
    if isnan(Nelem) && old<Nelem_vec(1)
        Nelem = Nelem_vec(1);
    elseif isnan(Nelem)
        Nelem = Nelem_vec(end);
    end
    indexN = find(Nelem_vec==Nelem);
    sf_vec = [0.5 0.65 0.80 0.950 1.10 1.250 1.40 1.55 1.70 1.85 2 2.30  2.750];
    sizefactor = sf_vec(indexN);    
        
    n = 500*sizefactor;

    % Sphere:
    particle = trisphere( n, diameter );
    
    %% Coupled spheres
    gap = diameter / (ab);
    p1 = shift( particle, [-(diameter+gap)/2, 0, 0]);
    p2 = shift( particle, [(diameter+gap)/2, 0, 0]);
    
    %  Initialize nanoellipsoid:
    p = comparticle( epstab, { p1, p2 }, [ 2, 1; 3, 1 ], 1, 2, op );
    p = closed(p,[1, 2]);
    
    disp(['Number of elements per particle: ' num2str(length(p.p{1}.faces))]);

end