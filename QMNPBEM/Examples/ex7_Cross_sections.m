%-------------------------------------------------------------------------
% Author: Alvaro Gomez Inesta (UPC, MIT) & Thomas Christensen (MIT).
% March 2017. (Last updated: March 19, 2018)
% Massachusetts Institute of Technology (MIT) & UPC (CFIS).
%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 7: compute cross section maps (frequency vs size) for a single 
% sphere. Each step in the loop takes about 15s.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear all

%% Constants
h = 6.5821e-16; % [eV*s] (h-bar)
c = 3e8; %[m/s]

%% Geometrical data:
diamvec = [1:1:25]; % Diameters [nm]
ab = 1; % Aspect ratio
ind_vec = [-0.144 -0.334 0.343]; % Direction of induced dipole
per_vec = [-0.144 -0.334 0.343]; % Direction of perturbation

%% Material:
wp = 2.2e15; %[rad/s]
w = [0.2:0.001:1.1]*wp; % Frequencies [rad/s]
k = w./c; % Wavenumber in vacuum [1/m]
beta = 1e6; % HDM parameter [m/s]
gamma = 1.8e13; % Loss rate [1/s]
eps_B = 1;
eps = eps_B - wp^2./(w.*(w+1i*gamma)); % Dielectric function: Drude metal
eps_d = 1; % Dielectric media

%% d-parameters:
d_T = 'HDM';
d_II = 'HDM';

%% Numerical data:
Nelem = 500; % Estimated number of elements
mode = 'all';
cs_scamat = [];
cs_absmat = [];


for ii = 1:length(diamvec)
    tic
    diam = diamvec(ii);
    disp(['Diameter: ' num2str(diam) ' nm'])

    %% Create spheroid
    disp('   -Creating sphere...')
    [p, op] = AgSpheroid(diam, ab, Nelem);

    %% Compute Lambdas
    disp('   -Computing Lambdas...')
    [Lambda_0, Lambda_T, Lambda_II, Sig, g] = Lambdas(p, op, mode);

    %% Compute polarizability
    disp('   -Computing polarizability...')
    alpha = polarizability(w, p, op, eps_d, eps, ind_vec,per_vec, d_T, d_II, Lambda_0, Lambda_T, Lambda_II, Sig, g, beta, wp);

    %% Absorption cross section
    cs_absmat = [cs_absmat; k.*eps_d^.5.*imag(alpha)];
    cs_scamat = [cs_scamat; (eps_d^.5*k).^4.*abs(alpha).^2./(6*pi)];
    toc
end
    
%% Plots
figure()
surf(w*h,diamvec,log(abs(cs_absmat)))
    view(2)
    colormap(morgenstemning)
    shading interp
    colorbar
    xlim([w(1)*h w(end)*h])
    ylim([diamvec(1) diamvec(end)])
    xlabel('Energy (eV)','FontSize',14)
    ylabel('2R (nm)','FontSize',14)
    title('Absorption cross section','FontSize',14)

figure()
surf(w*h,diamvec,log(abs(cs_scamat)))
    view(2)
    colormap(morgenstemning)
    shading interp
    colorbar
    xlim([w(1)*h w(end)*h])
    ylim([diamvec(1) diamvec(end)])    
    xlabel('Energy (eV)','FontSize',14)
    ylabel('2R (nm)','FontSize',14)   
    title('Scattering cross section','FontSize',14)

        
        
        
        
        