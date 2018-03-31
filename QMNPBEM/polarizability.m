%-------------------------------------------------------------------------
% Author: Alvaro Gomez Inesta (UPC, MIT) & Thomas Christensen (MIT).
% March 2018. (Last updated: March 26, 2018)
% Massachusetts Institute of Technology (MIT) & UPC (CFIS).
%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical computation of the polarizability of a 
% system of particles (Footnote [71] from ref. [1]).
%
% For more information about these quantum corrections see:
%  [1] Christensen, T., Yan, W., Jauho, A. P., Soljacic, M., & Mortensen, 
%      N. A. (2017). Quantum corrections in nanoplasmonics: shape, scale, 
%      and material. Physical Review Letters, 118(15), 157402.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% INPUTS:
%   w         : compute the polarizability at this frequency/ies. It should
%               be a row vector. Units (recommended): [rad/s].
%   p         : composite particle (see comparticle in MNPBEM).
%   op        : options (see MNPBEM).
%   ind_vec   : induced dipole direction [x y z].
%   per_vec   : pertubation direction [x y z].
%   eps       : row vector containing the dielectric function at the
%               frequencies given by the input w. Units: adimensional.
%   d_T       : row vector containing the perpendicular d-parameter at the 
%               frequencies given by the input w. Use 'HDM' to employ the
%               built-in hydrodynamic (HDM) approximation.
%               Units (recommended): [m].
%   d_II      : row vector containing the parallel d-parameter at the 
%               frequencies given by the input w. Use 'HDM' to employ the
%               built-in hydrodynamic (HDM) approximation.
%               Units (recommended): [m].
%   Lambda_0  : classical eigenvalue (see [1]). Units: adimensional.
%   Lambda_T  : perpendicular Lambda (see [1]). Units (recommended):
%               [1/nm].
%   Lambda_II : parallel Lambda (see [1]). Units (recommended): [1/nm].
%   Sig       : charge distribution for the selected mode. Irrelevant
%               units.
%   beta      : hydrodynamic parameter of the material. Only required if
%               d_T or d_II = 'HDM'. Units (recommended): [m/s].

% OUTPUTS:
%   alpha     : polarizability at the desired values of frequency.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% FURTHER IMPROVEMENTS %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%warning('GET eps_d AND eps OF EACH PARTICLE FROM p')
%warning('w IS ONLY NEEDED IF WE COMPUTE d_T OR eps')
%warning('THE MAGNITUDE OF ALPHA HAS NOT BEEN VERIFIED BUT THE FREQ-DEPENDENCE IS OK')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
function [alpha] = polarizability(w, p, op, eps_d, eps, ind_vec, per_vec, d_T, d_II, Lambda_0, Lambda_T, Lambda_II, Sig, g, beta, wp)
    warning('This is not the last version of polarizability() function (some features not implemented here).')
    %% Check arguments
    if nargin > 16
        error('Error: Too many inputs in polarizability().');
    end
    
    % Fill in unset optional values:
    switch nargin
        case 1
            error('Error in polarizability(): input(s) missing.')
        case 2
            error('Error in polarizability(): input(s) missing.')
        case 3
            error('Error in polarizability(): input(s) missing.')
        case 4
            error('Error in polarizability(): input(s) missing.')
        case 5
            error('Error in polarizability(): input(s) missing.')
        case 6
            error('Error in polarizability(): input(s) missing.')
        case 7
            error('Error in polarizability(): input(s) missing.')
        case 8
            error('Error in polarizability(): input(s) missing.')
        case 9
            error('Error in polarizability(): input(s) missing.')
        case 10
            error('Error in polarizability(): input(s) missing.')
        case 11
            error('Error in polarizability(): input(s) missing.')
        case 12
            error('Error in polarizability(): input(s) missing.')
        case 13
            error('Error in polarizability(): input(s) missing.')
        case 14
            if strcmp(d_T,'HDM') || strcmp(d_II,'HDM')
                error('Error: input(s) missing (beta, wp).')
            end
        case 15
            if strcmp(d_T,'HDM') || strcmp(d_II,'HDM')
                error('Error: input(s) missing (wp).')
            else
                warning('beta provided but not needed.')
            end
        case 16
            if strcmp(d_T,'HDM') || strcmp(d_II,'HDM')
            else
                warning('beta and wp provided but not needed.')
            end
    end
    
    %% Pre-process input values
    % Number of particles:
    P = length(p.p);
    
    %Convert Lambdas from [nm] to [m]
    Lambda_T = Lambda_T*1e9; % [1/nm] to [1/m]
    Lambda_II = Lambda_II*1e9; % [1/nm] to [1/m]
   
    % Set w as row vector
    sizew = size(w);
    if sizew(1)~=1 && sizew(2)==1
        w = w.';
    elseif sizew(1)~=1 && sizew(2)~=1
        error('Error in polarizability(): "w" is a matrix and should be a vector')
    end
    
    % Normalize induced dipole and perturbation directions
    ind_vec = ind_vec/norm(ind_vec);
    per_vec = per_vec/norm(per_vec);
    
    % Perpendicular d-parameter
    if strcmp(d_T,'HDM') % HDM value
        %d_T = (-beta ./ (wp^2-w.^2).^.5);
        d_T = (-beta ./ (wp^2-w.^2).^.5);
    else
        if length(w)~=length(d_T)
            error('Error in polarizability(): "w" and "d_T" must have same length')
        end
    end
    % Parallel d-parameter
    if strcmp(d_II,'HDM') % HDM value
        d_II = zeros(1,length(w));
    else
        if length(w)~=length(d_T)
            error('Error in polarizability(): "w" and "d_T" must have same length')
        end
    end
    
    %% Previous definitions
    % Number of modes
    sizeSig = size(Sig);
    M = sizeSig(2);
    
    % Number of frequency points
    Nw = length(w);
    
    % Areas of the boundary elements
    Areas_vec = [];
    for part = 1:P
        Areas_vec = [Areas_vec; p.p{1,part}.area];
    end
    Areas = spdiag(Areas_vec);
    
    % Dielectric function in vacuum (Electrostatic CGS)
    eps0 = 1/(4*pi); 
    
    % Position of each boundary element
    rvec = [];
    for part = 1:P
        rvec = [rvec; p.p{1,part}.pos];
    end
    
    % Normal vectors to the boundary elements
    aux1 = zeros(length(Areas_vec),1); % Auxiliar field
    [~,~,t1,t2] = deriv(p,aux1);
    nvec = cross( t1, t2 );
    % To plot these vectors, use the following commands:
        %>>plot(p,'cone',nvec,'scale',0.1)
        %>>plot(p,'EdgeColor','b')
        
    % Lambda
    Lambda = (eps_d + eps)./(eps_d-eps);
    
    %% alpha^0_n: column vector with alpha^0 for each mode in each row
    %alpha0 is INDEPENDENT of the frequency
    ind_mat = repmat(ind_vec,[M,1]);
    per_mat = repmat(per_vec,[M,1]);
    ri = dot(rvec,ind_mat,2);
    nj = dot(nvec,per_mat,2);

    phi_sig = (1/(4*pi*eps0))*Sig'*g'*Areas*Sig; % Bra-ket potential-charge 
                                                 %(matrix: sig modes in
                                                 %each column and phi modes
                                                 %in each row)
    phi_nj = ((1/(4*pi*eps0))*Sig'*g'*Areas*nj);
    ri_sig = (ri'*Areas*Sig).';
    phi_sig_nn = diag(phi_sig); % Bra-ket potential_n-charge_n
    sig_sig = Sig'*Areas*Sig; % Bra-ket charge-charge
    
    alpha0 = phi_nj.*ri_sig./phi_sig_nn;

    %% alpha^1_n: column vector with alpha^1 for each mode in each row
    %alpha1 DEPENDS on the frequency
    
    sum_n_T = zeros(M,Nw); % each element is the sum over m=~n for a
                           %particular mode n
    sum_n_II = zeros(M,Nw); % each element is the sum over m=~n for a
                            %particular mode n
                                  
    for m = 1:M % loop over the total number of modes
        Lambda_0_n = repmat([Lambda_0(1:m-1); Lambda_0(m+1:end)],[1,Nw]);
        phi_sig_mm_n = repmat([phi_sig_nn(1:m-1); phi_sig_nn(m+1:end)],[1,Nw]);
        sig_sig_nm_n = repmat([sig_sig(1:m-1,m); sig_sig(m+1:end,m)],[1,Nw]);
        phi_nj_n = repmat([phi_nj(1:m-1); phi_nj(m+1:end)],[1,Nw]);
        Lambda_n = repmat(Lambda,[M-1,1]);
        
        sum_mat_T = sig_sig_nm_n.*phi_nj_n./((Lambda_0_n-Lambda_n).*phi_sig_mm_n);
        sum_n_T(m,:) = sum(sum_mat_T)*(Lambda_0(m)-1)*(Lambda_0(m)+1)/(2*eps0);
        sum_n_II(m,:) = sum_n_II(m,:)*0;
    end

% This is a try to reduce the computational cost of the previous loop. It
%didnt work. Left here since it may be useful for somebody:
%     Lambda_0_n = zeros(M-1,Nw,M);
%     phi_sig_mm_n = zeros(M-1,Nw,M);
%     sig_sig_nm_n = zeros(M-1,Nw,M);
%     phi_nj_n = zeros(M-1,Nw,M);
%     Lambda_n = zeros(M-1,Nw,M);
%     for m = 1:M % loop over the total number of modes
%         Lambda_0_n(:,:,m) = repmat([Lambda_0(1:m-1); Lambda_0(m+1:end)],[1,Nw]);
%         phi_sig_mm_n(:,:,m) = repmat([phi_sig_nn(1:m-1); phi_sig_nn(m+1:end)],[1,Nw]);
%         sig_sig_nm_n(:,:,m) = repmat([sig_sig(1:m-1,m); sig_sig(m+1:end,m)],[1,Nw]);
%         phi_nj_n(:,:,m) = repmat([phi_nj(1:m-1); phi_nj(m+1:end)],[1,Nw]);
%         Lambda_n(:,:,m) = repmat(Lambda,[M-1,1]);
%     end
%         sum_mat_T = sig_sig_nm_n.*phi_nj_n./((Lambda_0_n-Lambda_n).*phi_sig_mm_n);
%         sum_n_T = reshape(sum(sum_mat_T,1),[M,Nw])*(Lambda_0(m)-1)*(Lambda_0(m)+1)/(2*eps0);
%         sum_n_II = sum_n_II*0;
    
    alpha1_T = (repmat(ri_sig,[1,Nw]).*sum_n_T)./repmat(phi_sig_nn,[1,Nw]);
    alpha1_II = (repmat(ri_sig,[1,Nw]).*sum_n_II)./repmat(phi_sig_nn,[1,Nw]);

    %% Polarizability (alpha)
    num = repmat(alpha0,[1,Nw])+repmat(d_T,[M,1]).*alpha1_T+repmat(d_II,[M,1]).*alpha1_II;
%     den = repmat(Lambda_0,[1,Nw])+repmat(d_T,[M,1]).*repmat(Lambda_T,[1,Nw])+...
%           repmat(d_II,[M,1]).*repmat(Lambda_II,[1,Nw])-repmat(Lambda,[M,1]); 
    den = repmat([0; Lambda_0(2:end)],[1,Nw])+repmat(d_T,[M,1]).*repmat(Lambda_T,[1,Nw])+...
          repmat(d_II,[M,1]).*repmat(Lambda_II,[1,Nw])-repmat(Lambda,[M,1]); 
    alpha = 2*sum(num./den);
    
    %% Study discrepancies with Apell
    F = repmat([0; Lambda_0(2:end)],[1,Nw])+repmat(d_T,[M,1]).*repmat(Lambda_T,[1,Nw])-...
        -repmat(Lambda,[M,1]);
    figure()
    semilogy(w*6.5821e-16,abs(real(F))); hold on
    title('Re(F)')
    figure()
    semilogy(w*6.5821e-16,abs(imag(F))); hold on
    title('Im(F)')

end