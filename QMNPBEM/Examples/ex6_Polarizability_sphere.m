%-------------------------------------------------------------------------
% Author: Alvaro Gomez Inesta (UPC, MIT) & Thomas Christensen (MIT).
% March 2018. (Last updated: March 26, 2018)
% Massachusetts Institute of Technology (MIT) & UPC (CFIS).
%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 6: compute polarizability of a single sphere. 
% The discrepancy between the classical spectrum and the quantum-corrected 
% spectra is larger for smaller radius.
% Running time with default parameters: ~15s.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear all

%% Constants
h = 6.5821e-16; % [eV*s] (h-bar)
c = 3e8; %[m/s]

%% Geometrical data:
diam = 40; % Diameter [nm]
R = diam/2; % Radius [nm]
ab = 1; % Aspect ratio
Nelem = 1000; % Approximated number of elements

%% Material data:
wp = 2.2e15;
w = (0.2:0.01:1.1)*wp;
beta = 1e6;
gamma = 1.8e13;
eps_B = 1;
eps = eps_B - wp^2./(w.*(w+1i*gamma));
Lambda = (1+eps)./(1-eps);
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
d_T = 'HDM';
d_II = 'HDM';

ind_vec = [0 0 1]; % Induced dipole direction
per_vec = [0 0 1]; % Perturbation direction

disp('   -Computing polarizability...')
alpha = polarizability(w, p, op, eps_d, eps, ind_vec, per_vec, d_T, d_II, Lambda_0, Lambda_T, Lambda_II, Sig, g, beta, wp);

%% Compute polarizability (Apell)
%Analytical expression using a modified dielectric function, derived from
%Feibelman's d-parameters theory. See [1] for further details.
%[1] Apell, P., & Ljungbert, Å. (1982). A general non-local theory for the 
%electromagnetic response of a small metal particle. Physica Scripta, 
%26(2), 113.

d_T = (-beta ./ (wp^2-w.^2).^.5 * 1e9); % [nm]
d_II = 0;

% Each semiaxis of the ellipsoid:
a1 = diam/2;
a2 = diam/2;
a3 = ab*diam/2;

% WARNING: the d-parameters used in [1] have a negative sign w.r.t. the
%definition we are using, so we need to include sign changes for the 
%d-parameters in the expression of the modified dielectric function.
%Moreover, d_II~=d_{theta}, so we cannot use the same definition (we omit
%the terms with d_{theta} since we are assuming d_II=0).

% Introducing Lambda=(1+eps)/(1-eps) in the expression for the 
%polarizability from [1], it yields:
alpha_apell = (4*pi/3)*a1*a2*a3*(2+2*d_T/R)./(-1/3-4*d_T/(3*R)-Lambda);

%% Compute polarizability (classical)
%See [2] for further info.
%[2] Novotny, L., & Hecht, B. (2012). Principles of nano-optics. Cambridge 
%university press.

% Again, rewriting the expression for the dipolar polarizability from [2]
%in terms of Lambda, one gets:
alpha_clas = (4*pi/3)*a1*a2*a3*2./(-1/3-Lambda);

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

eps_d_mod = (1+dl)*eps_d; % Modified dielectric function of the external medium

alpha_HDM = (4*pi/3)*a1*a2*a3*3*(eps-eps_d_mod)./(eps+2*eps_d_mod);


%% Plot polarizabilities and error QMNPBEM vs Apell (real)    
figure()
ax1 = axes('YAxisLocation','right','Color','none');
errvec = abs((real(alpha)-real(alpha_apell))./real(alpha_apell)*100);
line(w*h,errvec,'color','r','Parent',ax1)
ax1.YColor = 'r';
ax1.YScale = 'log';
ylabel('Error (%)','FontSize',14)   

ax2 = axes('Position',ax1.Position,'YAxisLocation','left','Color','none');
ylabel('Error (%)','FontSize',14)     
line(w*h,real(alpha_clas),'color',[0 0.8 0],'LineWidth',2,'Parent',ax2); hold on
line(w*h,real(alpha),'color','k','LineWidth',2,'Parent',ax2); hold on
line(w*h,real(alpha_apell),'color',[0.2 0.5 0.9],'LineStyle',':','LineWidth',3,'Parent',ax2); hold on
line(w*h,real(alpha_HDM),'color',[0.7 0 0],'LineStyle','-.','LineWidth',2,'Parent',ax2); hold on
xlabel('Energy (eV)','FontSize',14)
ylabel('Re(\alpha(\omega))','FontSize',14)   
title(['Polarizability REAL - ab = ' num2str(ab) ', R = ' num2str(diam/2) ' nm, Nelem = ' num2str(Nelem)],'FontSize',14)
legend({'Classical','QMNPBEM','Apell','HDM'},'FontSize',14,'Location','NorthWest')

linkprop([ax1, ax2],{'XLim'}); % This allows zooming x-axis for both axis

%% Plot polarizabilities and error QMNPBEM vs Apell (imag)
figure()
ax1 = axes('YAxisLocation','right','Color','none');
errvec = abs((imag(alpha)-imag(alpha_apell))./imag(alpha_apell)*100);
line(w*h,errvec,'color','r','Parent',ax1)
ax1.YColor = 'r';
ax1.YScale = 'log';
ylabel('Error (%)','FontSize',14)   

ax2 = axes('Position',ax1.Position,'YAxisLocation','left','Color','none');
ylabel('Error (%)','FontSize',14)     
line(w*h,imag(alpha_clas),'color',[0 0.8 0],'LineWidth',2,'Parent',ax2); hold on
line(w*h,imag(alpha),'color','k','LineWidth',2,'Parent',ax2); hold on
line(w*h,imag(alpha_apell),'color',[0.2 0.5 0.9],'LineStyle',':','LineWidth',3,'Parent',ax2); hold on
line(w*h,imag(alpha_HDM),'color',[0.7 0 0],'LineStyle','-.','LineWidth',2,'Parent',ax2); hold on
xlabel('Energy (eV)','FontSize',14)
ylabel('Im(\alpha(\omega))','FontSize',14)   
title(['Polarizability IMAG - ab = ' num2str(ab) ', R = ' num2str(diam/2) ' nm, Nelem = ' num2str(Nelem)],'FontSize',14)
legend({'Classical','QMNPBEM','Apell','HDM'},'FontSize',14,'Location','NorthWest')
ax2.YScale = 'log';

linkprop([ax1, ax2],{'XLim'}); % This allows zooming x-axis for both axis