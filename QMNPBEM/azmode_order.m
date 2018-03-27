%-------------------------------------------------------------------------
% Author: Alvaro Gomez Inesta (UPC, MIT) & Thomas Christensen (MIT).
% March 2017. (Last updated: March 3, 2018)
% Massachusetts Institute of Technology (MIT) & UPC (CFIS).
%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliar function to Lambdas(). 
% Select the mode with larger azimuthal symmetry in particle 'azpart'.
% This is a BETA FUNCTION, may not work as desired in some cases.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% INPUTS:
%   p       : composite particle (see comparticle in MNPBEM).
%   Sig     : matrix containing the charge distribution of the modes
%             in columns.
%   dirs    : direction of longitudinal axis in particle 'azpart'.
%   azpart  : number of the particle to analyze.

% OUTPUTS:
%   mode    : mode with larger value (see definition below).


function index = azmode_order(p,Sig,dirs,azpart)
    
    % TRIANGULAR ELEMENTS
    
    % Geometric data:
    geom = p.p{azpart};
    faces = geom.faces(:,1:3);
    verts = geom.verts;
    
    % Matrix with the x-y-z coordinates of each vertex:
    % Xv(i,:) = [x1i y1i z1i x2i y2i z2i x3i y3i z3i], i-th element
    Xv = [];
    for i = 1:length(faces)
        Xv(i,:) = [verts(faces(i,1),:) verts(faces(i,2),:) verts(faces(i,3),:)];
    end
    
    % Coordinates matrix (centroids of the elements):
    X = [(Xv(:,1)+Xv(:,4)+Xv(:,7))/3 (Xv(:,2)+Xv(:,5)+Xv(:,8))/3 (Xv(:,3)+Xv(:,6)+Xv(:,9))/3];
    
    % 'dirs' CAN ONLY BE Z FOR NOW
    
    % Elements ordered in ascending coordinate of dirs:
    [~,ind] = sort(X(:,3));
    X = X(ind,:);
    
    % 'Sig' with same order as X:
    Sig = Sig(ind,:);
    
    % For each mode (columns in Sig), compute the sum of charge density for
    % all elements with the same coordinate in 'dirs' direction:
    B = [];
    g = 0; % group number
    zold = X(1,3)+1; % old z (+1 to make it different to the first 'z')
    for jj = 1:length(X)
        z = X(jj,3);
        if z == zold
            B(g,:) = B(g,:) + Sig(jj,:);
        else
            B(g+1,:) = Sig(jj,:);
            g = g+1;
        end
        zold = z;
    end
    
    % Val(i) = value of mode i. This value represents how likely is the
    % mode to have an azimuthal symmetry.
    Val = sum(B.^2,1); % Sum elements of B squared to avoid cancellations.
    
    % Order modes in decreasing Val.
    [~,index] = sort(real(Val),'ascend');
    
end
