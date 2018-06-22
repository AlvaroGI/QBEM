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
%   edge    : edge of the triangles [nm].
%   height  : height of the triangles [nm].
%   gap     : gap between triangles [nm].
%   round_R : round radius [1/edge].
%   Nelem   : approximated number of elements PER SPHERE.

% OUTPUTS:
%   p       : composite particle (see comparticle in MNPBEM).
%   op      : options (see MNPBEM).


function [p, op] = AgCoupledBowtie(edge,height,gap,round_R,max_element_size)

    op = bemoptions( 'sim', 'stat', 'waitbar', 0, 'interp', 'curv');
    
    % Dielectric functions
    eps_d = 1;
    epstab = { epsconst( eps_d ), epsdrude( 'Ag' ), epsdrude( 'Ag' ) };
    
    %% Single triangle
    % Geometry:
    nz = ceil(edge/height)*2;
    %  Create polygon
    h = (edge^2-(edge/2)^2)^.5;
    poly = round( polygon( 3, 'size', [ h,  2 / sqrt( 3 ) * h ], 'dir', 1 ), 'rad', round_R );  % 'dir' = 1 -> normal vectors outwards    
    %poly = polygon( 3, 'size', [ h,  2 / sqrt( 3 ) * h ], 'dir', 1 )  % 'dir' = 1 -> normal vectors outwards    
    %  Extrude polygon to particle
    zd = edgeprofile( height, nz ); % round edges   
    particle = tripolygon( poly, zd, 'hdata',   struct( 'hmax', max_element_size ));
    
    %% Coupled triangles
    xdisp = edge/3^.5;
    p1 = shift( particle, [xdisp+gap/2+0.0056, 0, 0]); % MANUALLY ADJUSTED 0.0056
    p2 = shift(flip(particle,1), [-xdisp-gap/2-0.0056, 0, 0]); % MANUALLY ADJUSTED 0.0056
    
    %  Initialize bowtie:
    p = comparticle( epstab, { p1, p2 }, [ 2, 1; 3, 1 ], 1, 2, op );
    p = closed(p,[1, 2]);
    
    disp(['Number of elements per particle: ' num2str(length(p.p{1}.faces))]);

end