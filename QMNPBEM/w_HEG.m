%-------------------------------------------------------------------------
% Author: Alvaro Gomez Inesta (UPC, MIT) & Thomas Christensen (MIT).
% March 2018. (Last updated: March 6, 2018)
% Massachusetts Institute of Technology (MIT) & UPC (CFIS).
%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical computation of the classical renonance frequency and the quantum
% first order correction for a system of particles, assuming an homogeneous 
% electron gas behavior (HEG hypothesis).
%
% For more information about these quantum corrections see:
%  [1] Christensen, T., Yan, W., Jauho, A. P., Soljacic, M., & Mortensen, 
%      N. A. (2017). Quantum corrections in nanoplasmonics: shape, scale, 
%      and material. Physical Review Letters, 118(15), 157402.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% INPUTS:
%   Lambda_0  : classical eigenvalue (see [1]). Units: adimensional.
%   Lambda_T  : perpendicular Lambda (see [1]). Units (recommended):
%               [1/nm].
%   Lambda_II : parallel Lambda (see [1]). Units (recommended): [1/nm].
%   wp        : plasma frequency of the material. Units (recommended):
%               [rad/s].
%   d_T       : 2 columns matrix containing the perpendicular d-parameter
%               (first column) at each frequency (second column). Use 'HDM' 
%               to employ the built-in hydrodynamic (HDM) approximation.
%               Units (recommended): [m].
%   d_II      : 2 columns matrix containing the parallel d-parameter
%               (first column) at each frequency (second column). Use 'HDM' 
%               to employ the built-in hydrodynamic (HDM) approximation.
%               Units (recommended): [m].
%   beta      : hydrodynamic parameter of the material. Only required if
%               d_T or d_II = 'HDM'. Units (recommended): [m/s].

% OUTPUTS:
%   w0    : renonance frequency for the classical approach. Units: [rad/s]
%           if recommended units for inputs were used.
%   w1    : first order quantum correction to the classical approach.
%           Units: [rad/s] if recommended units for inputs were used.
    
function [w0,w1] = w_HEG(Lambda_0, Lambda_T, Lambda_II, wp, d_T, d_II, beta)
    %% Pre-process input values
    if nargin > 7
        error('Error: Too many inputs in Lambdas().');
    end
    
    % Fill in unset optional values:
    switch nargin
        case 1
            error('Error: input(s) missing.')
        case 2
            error('Error: input(s) missing.')
        case 3
            error('Error: input(s) missing.')
        case 4
            error('Error: input(s) missing (d-parameters).')
        case 5
            error('Error: input(s) missing (d-parameters).')
        case 6
            if strcmp(d_T,'HDM') || strcmp(d_II,'HDM')
                error('Error: input(s) missing (beta).')
            end
    end
    
    %% Classical resonance (w0)
    w0 = wp*((1+Lambda_0)/2).^.5;
    
    %% d-parameters at frequency w0
    % Perpendicular d-parameter
    if strcmp(d_T,'HDM') % HDM value
        d_T0 = -beta ./ (wp^2-w0.^2).^.5;
    else
        [~,index] = min(abs( d_T(:,2)-w0 )); % Find in d_T(:,2) the
                                             % closest frequency to w0
        d_T0 = d_T(index,1);
    end
    % Parallel d-parameter
    if strcmp(d_II,'HDM') % HDM value
        d_II0 = 0;
    else
        [~,index] = min(abs( d_II(:,2)-w0 )); % Find in d_T(:,2) the
                                              % closest frequency to w0
        d_II0 = d_II(index,1);
    end
    
    %% Quantum first order correction
    Lambda_T = Lambda_T*1e9; % [1/nm] to [1/m]
    Lambda_II = Lambda_II*1e9; % [1/nm] to [1/m]
    w1 = 0.25*(Lambda_T.*d_T0 + Lambda_II.*d_II0) .* (wp^2 ./ w0);
    
end
    
    
    
    
    
    
    
    
    
    
    
    