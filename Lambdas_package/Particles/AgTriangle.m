%-------------------------------------------------------------------------
% Author: Alvaro Gomez Inesta (UPC, MIT) & Thomas Christensen (MIT).
% December 2017. (Last updated: January 20, 2018)
% Massachusetts Institute of Technology (MIT) & UPC (CFIS).
%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shortcut to create an Ag triangle embedded in air using MNPBEM for
% non-retarded static simulations and 'curv' interpolation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% INPUTS:
%   edge    : edge of the equilateral triangle [nm].
%   ab      : edge/height.
%   Nelem   : approximated number of elements.

% OUTPUTS:
%   p       : composite particle (see comparticle in MNPBEM).
%   op      : options (see MNPBEM).


function [p, op] = AgTriangle(edge,ab,Nelem)
    
    op = bemoptions( 'sim', 'stat', 'waitbar', 0, 'interp', 'curv');
  
    % Dielectric function:
    eps_d = 1;
    epstab = { epsconst( eps_d ), epstable( 'silver.dat' ) };

    % Nelem to sizefactor:
    Nelem_vec = [450 548 924 1004 1088 1504 3112 3220 3328 3476 3720 3892 4228 5188 6296 12736 12952];
    old = Nelem;
    Nelem = interp1(Nelem_vec,Nelem_vec,Nelem,'nearest');
    if isnan(Nelem) && old<Nelem_vec(1)
        Nelem = Nelem_vec(1);
    elseif isnan(Nelem)
        Nelem = Nelem_vec(end);
    end
    indexN = find(Nelem_vec==Nelem);
    sf_vec = [0.5:0.15:3];
    sizefactor = sf_vec(indexN);  
    
    % Geometry:
    height = edge/ab;
    nz = sizefactor*8;
    %  Create polygon
    h = (edge^2-(edge/2)^2)^.5;
    poly = round( polygon( 3, 'size', [ h,  2 / sqrt( 3 ) * h ], 'dir', 1 ) );  % 'dir' = 1 -> normal vectors outwards    
    %  Extrude polygon to particle
    zd = edgeprofile( height, nz ); % round edges   
    particle = tripolygon( poly, zd, 'hdata',   struct( 'hmax', 0.05/sizefactor ));

    %  Initialize comparticle:
    p = comparticle( epstab, { particle }, [ 2, 1 ], 1, op );
    p = closed(p,1);
    
    disp(['Number of elements: ' num2str(length(p.p{1}.faces))]);

  
end