%-------------------------------------------------------------------------
% Author: Alvaro Gomez Inesta (UPC, MIT) & Thomas Christensen (MIT).
% March 2017. (Last updated: March 26, 2018)
% Massachusetts Institute of Technology (MIT) & UPC (CFIS).
%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 6: compute polarizability and cross section of a single sphere. 
% The discrepancy between the classical resonance and the first order 
% corrected is larger for smaller radius.
% Running time with default parameters: ~15s.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear all

%% Constants
h = 6.5821e-16; % [eV*s] (h-bar)
c = 3e8; %[m/s]

%% Geometrical data:
diam = 10; % Diameter [nm]
ab = 1; % Sphere
Nelem = 500; % Approximated number of elements

%% Material data:
wp = 2.2e15;
w = [0.2:0.001:1.1]*wp;
beta = 1e6;
gamma = 1.8e13;
eps_B = 1;
eps = eps_B - wp^2./(w.*(w+1i*gamma));
eps_d = 1;

%% Create spheroid
disp('   -Creating sphere...')
[p, op] = AgSpheroid(diam, ab, Nelem);

%% Compute Lambdas
mode = 'all';
disp('   -Computing Lambdas...')

[Lambda_0, Lambda_T, Lambda_II, Sig, g] = Lambdas(p, op, mode);

figure()
plot(p,Sig)
colormap(bluewhitered_mod())

%% Compute polarizability (QMNPBEM)
w = [0.2:0.001:1.1]*wp;

d_T = 'HDM';
d_II = 'HDM';

ind_vec = [-0.144 -0.334 0.343];
per_vec = [-0.144 -0.334 0.343];

disp('   -Computing polarizability...')
alpha = polarizability(w, p, op, eps_d, eps, ind_vec, per_vec, d_T, d_II, Lambda_0, Lambda_T, Lambda_II, Sig, g, beta, wp);

%% Compute polarizability (d-param)
%Analytical expression using a modified dielectric function, derived from
%Feibelman's d-parameters theory. See [1] for further details.
%[1] Apell, P., & Ljungbert, Å. (1982). A general non-local theory for the 
%electromagnetic response of a small metal particle. Physica Scripta, 
%26(2), 113.

d_T = (-beta ./ (wp^2-w.^2).^.5 * 1e9); % [nm]
d_II = 0;
l = 1;

% WARNING: the d-parameters used in [1] have a negative sign w.r.t. the
%definition we are using, so we need to include sign changes for the 
%d-parameters in the expression of the modified dielectric function.
%Moreover, d_II~=d_{theta}, so we cannot use the same definition (we omit
%the terms with d_{theta} since we are assuming d_II=0):
eps_mod = eps./(1+(l*2/diam)*(eps-1).*(-d_T));
alpha_dp = (diam/2)^(2*l+1)*(eps_mod - 1)./(eps_mod + (l+1)/l);
%alpha_dp = 4*pi*(diam/2)^(2*l+1)*(l*eps_mod-l*eps_d)./(l*eps_mod+(l+1)*eps_d);


%% Compute polarizability (classical)
%See [2] for further info.
%[2] Novotny, L., & Hecht, B. (2012). Principles of nano-optics. Cambridge 
%university press.

eps0 = 1/(4*pi); 
eps_d = 1;
l = 1;
alpha_clas = 4*pi*(diam/2)^(2*l+1)*(l*eps-l*eps_d)./(l*eps+(l+1)*eps_d);

%% Compute polarizability (HDM all orders)
%See [3] for further info.
%[3] Christensen, T., Yan, W., Raza, S., Jauho, A. P., Mortensen, N. A., & 
%Wubs, M. (2014). Nonlocal response of metallic nanospheres probed by 
%light, electrons, and atoms. Acs Nano, 8(2), 1745-1758.
l = 1;
k0 = w./c;
kM = k0.*eps.^.5;
kNL = (wp./beta) .* (eps ./ (eps_B*(eps_B-eps))).^.5;

xM = kM*diam/2*1e-9;
xNL = kNL*diam/2*1e-9;

jlM = (pi./(2*xM)).^.5.*besselj(l+0.5,xM);
jlNL = (pi./(2*xNL)).^.5.*besselj(l+0.5,xNL);
jl1NL = (pi./(2*xNL)).^.5.*besselj(l+0.5+1,xNL);
djlNL = (l./xNL).*jlNL-jl1NL;

Dl = l*(l+1)*jlM.*(eps-eps_B).*jlNL./(eps_B.*xNL.*djlNL);
dl = Dl ./ (jlM*(l+1));

num = l*(eps-(1+dl)*eps_d);
den = l*eps + (l+1)*(1+dl)*eps_d;
alpha_HDM = 4*pi*(diam/2)^(2*l+1)*(num)./(den);

%% Plot polarizabilities real and imag
figure()
plot(w*h,real(alpha_clas)/max(real(alpha_clas)),'color',[0 0.8 0],'LineWidth',2); hold on
plot(w*h,real(alpha_dp)/max(real(alpha_dp)),'*','color',[0.2 0.5 0.9],'LineWidth',1); hold on
plot(w*h,real(alpha)/max(real(alpha)),'k','LineWidth',2); hold on
plot(w*h,real(alpha_HDM)/max(real(alpha_HDM)),'-.','color',[1 0 0],'LineWidth',2); hold on

xlabel('Energy (eV)','FontSize',14)
ylabel('Normalized Re(\alpha(\omega))','FontSize',14)   
title(['Polarizability REAL - single sphere (R = ' num2str(diam/2) ' nm)'],'FontSize',14)
legend({'Classical','Analytical (d-param)','QMNPBEM','HDM'},'FontSize',14)

figure()
plot(w*h,imag(alpha_clas)/max(imag(alpha_clas)),'color',[0 0.8 0],'LineWidth',2); hold on
plot(w*h,imag(alpha_dp)/max(imag(alpha_dp)),'*','color',[0.2 0.5 0.9],'LineWidth',1); hold on
plot(w*h,imag(alpha)/max(imag(alpha)),'k','LineWidth',2); hold on
plot(w*h,imag(alpha_HDM)/max(imag(alpha_HDM)),'-.','color',[1 0 0],'LineWidth',2); hold on

xlabel('Energy (eV)','FontSize',14)
ylabel('Normalized Im(\alpha(\omega))','FontSize',14)   
title(['Polarizability IMAG - single sphere (R = ' num2str(diam/2) ' nm)'],'FontSize',14)
legend({'Classical','Analytical (d-param)','QMNPBEM','HDM'},'FontSize',14)

%% Plot polarizabilities
figure()
plot(w*h,abs(alpha_clas)/max(abs(alpha_clas)),'color',[0 0.8 0],'LineWidth',2); hold on
plot(w*h,abs(alpha_dp)/max(abs(alpha_dp)),'*','color',[0.2 0.5 0.9],'LineWidth',1); hold on
plot(w*h,abs(alpha)/max(abs(alpha)),'k','LineWidth',2); hold on
plot(w*h,abs(alpha_HDM)/max(abs(alpha_HDM)),'-.','color',[1 0 0],'LineWidth',2); hold on

xlabel('Energy (eV)','FontSize',14)
ylabel('Normalized |\alpha(\omega)|','FontSize',14)   
title(['Polarizability - single sphere (R = ' num2str(diam/2) ' nm)'],'FontSize',14)
legend({'Classical','Analytical (d-param)','QMNPBEM','HDM'},'FontSize',14)

%% Absorption cross section
figure()
k = w./c;
cs_abs = k.*eps_d^.5.*imag(alpha);
cs_abs_dp = k.*eps_d^.5.*imag(alpha_dp);
cs_abs_clas = k.*eps_d^.5.*imag(alpha_clas);
cs_abs_HDM = k.*eps_d^.5.*imag(alpha_HDM);

plot(w*h,abs(cs_abs_clas)/max(abs(cs_abs_clas)),'color',[0 0.8 0],'LineWidth',2); hold on
plot(w*h,abs(cs_abs_dp)/max(abs(cs_abs_dp)),'*','color',[0.2 0.5 0.9],'LineWidth',1); hold on
plot(w*h,abs(cs_abs)/max(abs(cs_abs)),'k','LineWidth',2); hold on
plot(w*h,abs(cs_abs_HDM)/max(abs(cs_abs_HDM)),'-.','color',[1 0 0],'LineWidth',2); hold on

xlabel('Energy (eV)','FontSize',14)
ylabel('Normalized \sigma_{abs}(\omega)','FontSize',14)   
title(['Absorption cross section - single sphere (R = ' num2str(diam/2) ' nm)'],'FontSize',14)
legend({'Classical','Analytical (d-param)','QMNPBEM','HDM'},'FontSize',14)
set(gca,'YScale','log')
 
%% Scattering cross section
figure()
k = w./c;
cs_sca = (eps_d^.5*k).^4.*abs(alpha).^2./(6*pi);
cs_sca_dp = (eps_d^.5*k).^4.*abs(alpha_dp).^2./(6*pi);
cs_sca_clas = (eps_d^.5*k).^4.*abs(alpha_clas).^2./(6*pi);
cs_sca_HDM = (eps_d^.5*k).^4.*abs(alpha_HDM).^2./(6*pi);

plot(w*h,abs(cs_sca_clas)/max(abs(cs_sca_clas)),'color',[0 0.8 0],'LineWidth',2); hold on
plot(w*h,abs(cs_sca_dp)/max(abs(cs_sca_dp)),'*','color',[0.2 0.5 0.9],'LineWidth',1); hold on
plot(w*h,abs(cs_sca)/max(abs(cs_sca)),'k','LineWidth',2); hold on
plot(w*h,abs(cs_sca_HDM)/max(abs(cs_sca_HDM)),'-.','color',[1 0 0],'LineWidth',2); hold on

xlabel('Energy (eV)','FontSize',14)
ylabel('Normalized \sigma_{sca}(\omega)','FontSize',14)   
title(['Scattering cross section - single sphere (R = ' num2str(diam/2) ' nm)'],'FontSize',14)
legend({'Classical','Analytical (d-param)','QMNPBEM','HDM'},'FontSize',14)
set(gca,'YScale','log')




%% Test
% %  BEM simulation
% bem = bemsolver( p, op );
% %  Plane wave excitation
% exc = planewave( [ 1, 0, 0 ], [ 0, 0, 1 ], op );
% %  Light wavelength in vacuum [nm]:
% enei = linspace(800, 1500, 100 );
% 
% sca_vec = zeros( length( enei ), 1 );
% ext_vec = zeros( length( enei ), 1 );
% abs_vec = zeros( length( enei ), 1 );
% 
% %  Loop over wavelengths
% polariz = [];
% for ien = 1 : length( enei )
%     %  surface charge
%     sig = bem \ exc( p, enei( ien ) );
%     sigma = sig.sig;
%     r = p.p{1}.pos;
%     dipp = (r'*spdiag(p.p{1}.area)*sigma)
%     polariz = [polariz dipp(1)/1];
%     
%     plot(p,sigma)
%     colorbar
%     pause
%     drawnow
% end
% 
% w = 2*pi*3e8 ./ enei *1e9;
% hw = w*h; % [eV]
% 
% plot(hw,polariz)

 
    