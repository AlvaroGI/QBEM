%-------------------------------------------------------------------------
% Author: Alvaro Gomez Inesta (UPC, MIT) & Thomas Christensen (MIT).
% December 2017. (Last updated: January 20, 2018)
% Massachusetts Institute of Technology (MIT) & UPC (CFIS).
%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shortcut to create an Ag rod embedded in air using MNPBEM for
% non-retarded static simulations and 'curv' interpolation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% INPUTS:
%   diameter: diameter of the rod [nm].
%   ab      : (total length of the rod)/diameter.
%   Nelem   : approximated number of elements.

% OUTPUTS:
%   p       : composite particle (see comparticle in MNPBEM).
%   op      : options (see MNPBEM).


function [p, op] = AgRod(diameter,ab,Nelem)

    op = bemoptions( 'sim', 'stat', 'waitbar', 0, 'interp', 'curv');
    
    % Dielectric functions:
    eps_d = 1;
    epstab = { epsconst( eps_d ), epstable( 'silver.dat' ) };
    
    % Geometry:
    height = ab*diameter; % [nm]
    % n1 for the circumference of the rod.   
    % n2 for the polar angles of the rod caps,
    % n3 for the cylinder-shaped middle part of the rod, and
    n1 = (Nelem/2)^.5;
    n2 = (1/(2+ab))*(Nelem/n1);
    n3 = n2*ab/1.5;

    % 'triangles' is a keyword that indicates that the surface is 
    %discretized through triangles rather than quadrilaterals.
    particle = trirod( diameter, height, [ n1, n2, n3 ], 'triangles' );
    
    %  Initialize rod:
    p = comparticle( epstab, { particle }, [ 2, 1 ], 1, op );
    p = closed(p,1);
    
    disp(['Number of elements: ' num2str(length(p.p{1}.faces))]);
    
end