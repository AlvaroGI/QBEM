%-------------------------------------------------------------------------
% Author: Alvaro Gomez Inesta (UPC, MIT) & Thomas Christensen (MIT).
% November 2017. (Last updated: March 19, 2018)
% Massachusetts Institute of Technology (MIT) & UPC (CFIS).
%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical computation of the Lambdas and charge distribution for a 
% system of particles.
%
% For more information about these quantum corrections see:
%  [1] Christensen, T., Yan, W., Jauho, A. P., Soljacic, M., & Mortensen, 
%      N. A. (2017). Quantum corrections in nanoplasmonics: shape, scale, 
%      and material. Physical Review Letters, 118(15), 157402.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% INPUTS:
%   p       : composite particle (see comparticle in MNPBEM).
%   op      : options (see MNPBEM).
%   mode    : selected mode (values from 1 to the number of nodes).
%             Set mode=0 to display the charge distribution for every
%             mode, the code will then ask you to select one of those modes
%             in the Command Window. 
%             Set mode=-1 for the same but showing the modes in decreasing
%             dipole moment (product of dipole moments of each particle)
%             along particular directions that the code will ask for.
%             Set mode='dipolar_max' to use the mode with larger product
%             of dipole moments in the directions specified by the user.
%             Set mode='azimuthal' to show the modes in decreasing order of
%             'stronger' azimuthal symmetry in particle 'azpart'. See 
%             azmode_order.m for further information. This is a BETA
%             FUNCTION, may not work as desired in some cases.
%             Set mode='all' to get the outputs for every mode (each column
%             in the output arrays corresponds to a different mode).
%   dirs    : direction of the total dipole moment of particle i in row
%             i (Px3 matrix, P = number of particles) for mode = -1 or
%             'dipolar_max. Direction of longitudinal axis for mode = 
%             'azimuthal'. Optional input, only required if mode = -1, 
%             'dipolar_max' or 'azimuthal'.
%   azpart  : particle whose modes' azimuthal symmetry will be taken into 
%             account. Input only required if mode = 'azimuthal'.
%¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡
% IMPORTANT: the maximization of the dipole moment is performed by
% maximizing the product of the dipole moments of the particles.
% We recommend to use mode=-1 to check if the desired mode is
% indeed the first mode so one can use mode='dipolar_max' instead.
%¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡

% OUTPUTS:
%   Lambda_0    : classical eigenvalue (see [1]). Units: adimensional.
%   Lambda_T    : perpendicular Lambda (see [1]). Units: [1/nm].
%   Lambda_II   : parallel Lambda (see [1]). Units: [1/nm].
%   Sig         : charge distribution for the selected mode.
%   g           : Green function over the finite elements of the system.

    
function [Lambda_0, Lambda_T, Lambda_II, Sig, g] = Lambdas(p, op, mode, dirs, azpart)
    %% Pre-process input values
    if nargin > 5
        error('Error: Too many inputs in Lambdas().');
    end

    % Number of particles:
    P = length(p.p);
    
    % Fill in unset optional values:
    switch nargin
        case 2
            error('Error in Lambdas(): mode not selected.')
        case 3
            if sum(mode)==-1 || strcmp(mode,'dipolar_max')
                warning('Directions of dipoles not specified. Set to [0 0 1] as default.');
                dirs = repmat([0 0 1],[P,1]);
            end
            if strcmp(mode,'azimuthal')
                warning('Direction of longitudinal axis not specified. Set to [0 0 1] as default.');
                dirs = [0 0 1];
                warning('Particle to analyze not specified. Set to 1 as default.');
                azpart = 1;
            end
        case 4
            if strcmp(mode,'azimuthal')
                if size(dirs,1)>P
                    error('Error in Lambdas(): too many directions specified.');
                end
                warning('Particle to analyze not specified. Set to 1 as default.');
                azpart = 1;
            elseif strcmp(mode,'all') || sum(mode)==0
                warning('Direction of dipoles specified but not required for this mode choice');
            else
                if size(dirs,1)<P
                    dirs = [dirs; repmat([0 0 1],[(P-size(dirs,1)),1])];
                    warning('Directions of all dipoles not specified. Set the remaining to [0 0 1] as default.');
                elseif size(dirs,1)>P
                    error('Error in Lambdas(): too many directions specified.');
                end
            end
        case 5
            if strcmp(mode,'azimuthal')
                if size(dirs,1)>P
                    error('Error in Lambdas(): too many directions specified.');
                end
                if length(azpart)~=1 || floor(azpart)~=azpart
                    error('Error in Lambdas(): "azpart" must be an integer.');
                end
            else
                if size(dirs,1)<P
                    dirs = [dirs; repmat([0 0 1],[(P-size(dirs,1)),1])];
                    warning('Directions of all dipoles not specified. Set the remaining to [0 0 1] as default.');
                elseif size(dirs,1)>P
                    error('Error in Lambdas(): too many directions specified.');
                end
                warning('Variable "azpart" specified but not needed.');
            end
            
    end
    
    if sum(mode)==-1 || strcmp(mode,'dipolar_max') || strcmp(mode,'azimuthal')
        % Normalize directions (required to compare dipole moments later):
        dirsn = (dirs(:,1).^2+dirs(:,2).^2+dirs(:,3).^2).^.5;
        dirs = dirs./repmat(dirsn,[1,3]);
    end
    
    %% Previous definitions
    Areas_vec = [];
    for part = 1:P
        Areas_vec = [Areas_vec; p.p{1,part}.area];
    end
    Areas = spdiag(Areas_vec);
    
    eps0 = 1/(4*pi); % Dielectric function in vacuum (Electrostatic CGS)
    
    %% Compute only up to the requested mode
    % This condition speeds up the solver dramatically when a particular
    % mode is used. However, it may cause trouble to find that particular
    % mode if mode=0 is not used. If the user wants to ensure the same mode
    % order as in mode=0, replace this condition by "nev = length(Areas);"
    % for any mode option.
    if sum(mode) == 0 || sum(mode) == -1 || strcmp(mode,'dipolar_max')|| strcmp(mode,'azimuthal') || strcmp(mode,'all')
        nev = length(Areas); % The max number of modes is the number of nodes
    else
        nev = mode;
    end

    %% Solver
    warning off
    bem = bemstateig( p, op, 'nev', nev );
    Lambdas_0 =  bem.ene;

    Sig = bem.ur; % The charge distribution are the right eigenvectors
    warning on

    % Biorthonormality: bem.ul*Areas*Sig = eye(nev)
    
    %% Order modes by decreasing classical eigenvalue (only if mode==0, -1,
    % 'dipolar_max', 'azimuthal' or 'all')
    if sum(mode) == 0 || sum(mode) == -1 || strcmp(mode,'dipolar_max')|| strcmp(mode,'azimuthal') || strcmp(mode,'all')
        Lambdas_0 = diag(flipud(diag(Lambdas_0)));
        Sig = fliplr(Sig);
    end
    
    %% Select mode
    while mode == 0
        fghd = figure();
        fignum = fghd.Number;
        plot(p,Sig)
        colormap(bluewhitered_mod())
        colorbar
        mode = input(['Select mode (1-' num2str(nev) '): ']);
        if floor(mode)~=mode || abs(mode)>nev || mode<=0
            disp(['Incorrect answer. Mode value should be integer between 1 and ' num2str(nev) '.'])
            mode = 0;
        end
        figure(fignum)
        close(fignum)
    end
    
    while mode == -1
        %Sort with decreasing product of dipole moments
        index = dipmode_order(p,Sig,dirs);
        Sig = Sig(:,index);
        
        fghd = figure();
        fignum = fghd.Number;
        plot(p,Sig)
        colormap(bluewhitered_mod())
        colorbar
        mode = input(['Select mode (1-' num2str(nev) '): ']);
        if floor(mode)~=mode || abs(mode)>nev || mode<=0
            disp(['Incorrect answer. Mode value should be integer between 1 and ' num2str(nev) '.'])
            mode = 0;
        end
        figure(fignum)
        close(fignum)
    end

    if strcmp(mode,'dipolar_max')
        mode = dipmode(p, Sig, dirs);
    end
 
    while strcmp(mode,'azimuthal')
        %Sort with decreasing product of dipole moments
        index = azmode_order(p,Sig,dirs,azpart);
        Sig = Sig(:,index);
        
        fghd = figure();
        fignum = fghd.Number;
        plot(p,Sig)
        colormap(bluewhitered_mod())
        colorbar
        mode = input(['Select mode (1-' num2str(nev) '): ']);
        if floor(mode)~=mode || abs(mode)>nev || mode<=0
            disp(['Incorrect answer. Mode value should be integer between 1 and ' num2str(nev) '.'])
            mode = 0;
        end
        figure(fignum)
        close(fignum)
    end
    
    %% Get all solutions or only a particular mode
    if strcmp(mode,'all')
        %% Classical Lambda0
        disp('   -Computing Lambda0...')
        Lambda_0 = diag(Lambdas_0)./(2*pi);

        %% Quantum Lambdas
        disp('   -Computing Lambda_T...')
        g = bem.g.G; % Green function    
        sig_sig = diag(Sig'*Areas*Sig); % Bra-ket charge-charge
        phi_sig = diag((1/(4*pi*eps0))*Sig'*g'*Areas*Sig); % Bra-ket potential-charge

        Lambda_T = (Lambda_0.^2-1).*sig_sig./(phi_sig*2*eps0);

        %% Lambda_II
        disp('   -Computing Lambda_II...')
        Lambda_II = ones(nev,1);
        for modenum = 1:nev % THIS LOOP COULD BE VECTORIZED TO SPEED UP CALCULATIONS
            Sigg = Sig(:,modenum);
            pot = (1/(4*pi*eps0))*g*Sigg; % Potential at the elements
            potv = interp(p,pot,'area'); % Potential at the vertices

            %  Computation of tangential gradient 
                 [dp1x,dp2x,t1,t2] = deriv(p,potv);
                 %  normal vector
                 nvec = cross( t1, t2 );
                 %  decompose into norm and unit vector
                 h = sqrt( dot( nvec, nvec, 2 ) );  nvec = bsxfun( @rdivide, nvec, h );
                 %  tangential gradient of V
                 grad = outer( bsxfun( @rdivide, cross( t2, nvec, 2 ), h ), dp1x ) -  ...
                        outer( bsxfun( @rdivide, cross( t1, nvec, 2 ), h ), dp2x );

            int_dpot2 = grad(:,1)'*Areas*grad(:,1) + grad(:,2)'*Areas*grad(:,2) + ...
                        grad(:,3)'*Areas*grad(:,3); % Bra-ket (tangential grad of
                                                    % potential)-(tangential
                                                    % grad of potential)

            Lambda_II(modenum) = 8*pi*eps0^2 * (int_dpot2) / (Sigg'*g'*Areas*Sigg);    
        end
    else
        Sig = Sig(:,mode); % Charge distribution of the mode
        %% Classical Lambda0
        Lambda_0 = diag(Lambdas_0)/(2*pi);
        Lambda_0 = Lambda_0(mode);   

        %% Quantum Lambdas
        g = bem.g.G; % Green function    
        sig_sig = Sig'*Areas*Sig; % Bra-ket charge-charge
        phi_sig = (1/(4*pi*eps0))*Sig'*g'*Areas*Sig; % Bra-ket potential-charge

        Lambda_T = (Lambda_0^2-1)*sig_sig/(phi_sig*2*eps0);

        %% Lambda_II
        pot = (1/(4*pi*eps0))*g*Sig; % Potential at the elements
        potv = interp(p,pot,'area'); % Potential at the vertices

        %  Computation of tangential gradient 
             [dp1x,dp2x,t1,t2] = deriv(p,potv);
             %  normal vector
             nvec = cross( t1, t2 );
             %  decompose into norm and unit vector
             h = sqrt( dot( nvec, nvec, 2 ) );  nvec = bsxfun( @rdivide, nvec, h );
             %  tangential gradient of V
             grad = outer( bsxfun( @rdivide, cross( t2, nvec, 2 ), h ), dp1x ) -  ...
                    outer( bsxfun( @rdivide, cross( t1, nvec, 2 ), h ), dp2x );

        int_dpot2 = grad(:,1)'*Areas*grad(:,1) + grad(:,2)'*Areas*grad(:,2) + ...
                    grad(:,3)'*Areas*grad(:,3); % Bra-ket (tangential grad of
                                                % potential)-(tangential
                                                % grad of potential)

        Lambda_II = 8*pi*eps0^2 * (int_dpot2) / (Sig'*g'*Areas*Sig);    
    end
end