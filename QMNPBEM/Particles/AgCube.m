%-------------------------------------------------------------------------
% Author: Alvaro Gomez Inesta (UPC, MIT) & Thomas Christensen (MIT).
% December 2017. (Last updated: January 20, 2018)
% Massachusetts Institute of Technology (MIT) & UPC (CFIS).
%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shortcut to create an Ag cube embedded in air using MNPBEM for
% non-retarded static simulations and 'curv' interpolation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% INPUTS:
%   len     : edge size [nm].
%   ab      : len/(fraction of the edge that is rounded). ab>2; ab=2 for a
%             sphere.
%   Nelem   : approximated number of elements.
%   shape   : 'tria' for triangular elements, 'cuad' for cuadrilateral
%             elements.

% OUTPUTS:
%   p       : composite particle (see comparticle in MNPBEM).
%   op      : options (see MNPBEM).


function [p, op] = AgCube(len,ab,Nelem,shape)
    % Dielectric functions
    eps_d = 1;
    epstab = { epsconst( eps_d ), epstable( 'silver.dat' ) };
    
    % Nelem to sizefactor:
    Nelem_vec = [384 600 864 1176 1536 1944 2400 2904 3456 4056 4704 5400 6144 6936 7776 8664 9600 10584 11616 12696 13824];
    Nelem = interp1(Nelem_vec,Nelem_vec,Nelem,'nearest');
    old = Nelem;
    Nelem = interp1(Nelem_vec,Nelem_vec,Nelem,'nearest');
    if isnan(Nelem) && old<Nelem_vec(1)
        Nelem = Nelem_vec(1);
    elseif isnan(Nelem)
        Nelem = Nelem_vec(end);
    end
    indexN = find(Nelem_vec==Nelem);
    sf_vec = [0.9:0.2:5];
    sizefactor = sf_vec(indexN);


    % n for grid size (number of elements in each dimension),
    % 'e' for the round-off parameter for the edges (e=1 for a sphere)
    n = sizefactor*10;
    e = 2/ab;
    
    particle = tricube( n, len, 'e', e );
    
    if shape == 'tria'
        % Only triangular elements:
        ff = totriangles(particle);
        particle.faces = [ff ones(length(ff),1)*NaN];
        op = bemoptions( 'sim', 'stat', 'waitbar', 0, 'interp', 'flat');
    elseif shape == 'cuad'
        op = bemoptions( 'sim', 'stat', 'waitbar', 0, 'interp', 'flat');
    else
        error('"shape" variable should be "tria" or "cuad"')
    end
    
    %  Initialize cube:
    p = comparticle( epstab, { particle }, [ 2, 1 ], 1, op );
    p = closed(p,1);
    
    disp(['Number of elements: ' num2str(length(p.p{1}.faces))]);

end