%-------------------------------------------------------------------------
% Author: Alvaro Gomez Inesta (UPC, MIT) & Thomas Christensen (MIT).
% December 2017. (Last updated: January 20, 2018)
% Massachusetts Institute of Technology (MIT) & UPC (CFIS).
%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shortcut to create an Ag spheroid embedded in air using MNPBEM for
% non-retarded static simulations and 'curv' interpolation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% INPUTS:
%   diameter: diameter of the spheroid [nm].
%   ab      : (total length of the spheroid)/diameter.
%   Nelem   : approximated number of elements.

% OUTPUTS:
%   p       : composite particle (see comparticle in MNPBEM).
%   op      : options (see MNPBEM).


function [p, op] = AgSpheroid(diameter,ab,Nelem)

    op = bemoptions( 'sim', 'stat', 'waitbar', 0, 'interp', 'curv');

    % Dielectric functions
    eps_d = 1;
    epstab = { epsconst( eps_d ), epstable( 'silver.dat' ) };    
        
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
    
    % Geometry:
    height = ab*diameter; % [nm]
    n = 500*sizefactor; % n is the number of boundary elements
    axiss = [diameter, diameter, height];

    particle = scale( trisphere( n, 1 ), axiss );
    
    %  Initialize spheroid:
    p = comparticle( epstab, { particle }, [ 2, 1 ], 1, op );
    p = closed(p,1);
    
    disp(['Number of elements: ' num2str(length(p.p{1}.faces))]);

end
