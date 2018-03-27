%-------------------------------------------------------------------------
% Author: Alvaro Gomez Inesta (UPC, MIT) & Thomas Christensen (MIT).
% December 2017. (Last updated: January 23, 2018)
% Massachusetts Institute of Technology (MIT) & UPC (CFIS).
%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliar function to Lambdas_coupled() and Lambdas_coupled_vec(). 
% Select the mode with larger dipolar moments in dir1 and dir2 direction.
% This function maximizes the product of dipolar moments.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% INPUTS:
%   p       : composite particle (see comparticle in MNPBEM).
%   Sig     : matrix containing the charge distribution of the modes
%             in columns.
%   dir1    : direction of dipole in particle 1.
%   dir2    : direction of dipole in particle 2.

% OUTPUTS:
%   mode    : mode with larger dipolar moment product.


function mode = dipmode(p,Sig,dirs)

    % Number of particles:
    P = length(p.p);
    
    % Store areas of the elements in a cell array
    Areas = cell(P,1);
    for ii = 1:P
        Areas{ii} = spdiag(p.p{1,ii}.area);
    end
    
    % Geometry
    r = cell(P,1); % Position of elements
    rp = cell(P,1); % Dipole direction
    rrp = cell(P,1); % r dot rp
    for ii = 1:P
        r{ii} = p.p{1,ii}.pos;
        rp{ii} = repmat(dirs(ii,:),[length(r{ii}),1]);
        rrp{ii} = dot(r{ii},rp{ii},2);
    end
    
    % Dipolar moment along each direction
    dip = zeros(P,length(Sig)); % Dipole moments of the particles
    elem = 0; % Number of elements of the previous particle
    for ii = 1:P
        dip(ii,:) = rrp{ii}'*Areas{ii}*Sig(elem+1:length(r{ii})+elem,:);
        elem = length(r{ii});
    end
    
    % Ignore non-physical modes (charged particles)
    elem = 0;
    for ii = 1:P % Loop over particles
        C = sum(Sig(elem+1:length(r{ii}),:)); % Total charge in the particle for each mode
        for m = 1:length(C) % Loop over modes
            if abs(C(m))>mean(abs(C)) % If the charge of the mode is larger than the mean value
                                      % (this is a way to check if the
                                      % particle is charged)
                dip(ii,m) = 0; % The dipole moment of the mode wont be taken into account
            end    
        end
        elem = length(r{ii});
    end
    
    Prod = prod(dip,1);
        
    
    [~, mode] = max(sign(real(Prod)) .* abs(Prod));
    
end