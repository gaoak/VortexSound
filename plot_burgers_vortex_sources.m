%% Burgers vortex: source-term distributions
% Output: ./burgers_vortex_sources.pdf

clear; clc; close all;

%% --------------------------- plot style -------------------------------
set(0,'defaultlinelinewidth',2);
set(0,'defaultaxeslinewidth',2);
set(0,'defaultaxesfontsize',16);
set(0,'defaultTextFontName', 'Times New Roman');
set(0,'defaultAxesFontName', 'Times New Roman');
set(0,'defaulttextfontsize',16);
set(0,'DefaultLineMarkerSize',9);

%% --------------------------- parameters -------------------------------
sigma = 1.0;
k     = 0.5;
Gamma = 4.0;
c0    = 10.0;

%% ------------------------- vortex fields ------------------------------
% avoid r = 0 singular point in intermediate expressions
r = linspace(1e-4*sigma, 6*sigma, 2200).';
xi = r/sigma;

V   = Gamma./(2*pi*r) .* (1 - exp(-(r.^2)/(sigma^2)));     % v_theta
Vp  = Gamma/(2*pi) .* ( ...
      (-1./r.^2).*(1 - exp(-(r.^2)/(sigma^2))) ...
    + (2/sigma^2).*exp(-(r.^2)/(sigma^2)) );
Vpp = Gamma/(2*pi) .* ( ...
      (2./r.^3).*(1 - exp(-(r.^2)/(sigma^2))) ...
    - (2./(sigma^2*r)).*exp(-(r.^2)/(sigma^2)) ...
    - (4*r/sigma^4).*exp(-(r.^2)/(sigma^2)) );

omega  = Vp + V./r;                      % omega_z
omegap = Vpp + Vp./r - V./(r.^2);        % d(omega_z)/dr

% material acceleration of base flow
ar = (k^2*r/4) - V.^2./r;
at = -(k/2).*(V + r.*Vp);

%% -------------------------- source terms ------------------------------
% Phillips equation RHS
S_phi = (3/2)*k^2 - 2*(V./r).*Vp;

% Howe equation RHS
S_howe_div = -(1./r).*gradient(r.*omega.*V, r);
S_howe_acc = -(1/c0^2) .* ( ar.*(-omega.*V) + at.*(omega.*(-k*r/2)) );
S_howe = S_howe_div + S_howe_acc;

% Dilatation equation RHS
trA2 = S_phi;
Q = -0.5*trA2;
dQdr = gradient(Q, r);
S_dil_4DQdt = 4*(-k*r/2).*dQdr;

Srth = 0.5*(Vp - V./r);
detS = k*(k^2/4 - Srth.^2);
S_dil_detS = -6*detS;
S_dil_stretch = -(3/2)*k*(omega.^2);
S_dil_acc = at.*(-omegap);
S_dil = S_dil_4DQdt + S_dil_detS + S_dil_stretch + S_dil_acc;

%% ----------------------------- plotting -------------------------------
position = [ ...
    0.10 0.70 0.85 0.24
    0.10 0.39 0.85 0.24
    0.10 0.08 0.85 0.24];

handle = figure;
set(handle,'Position',[0 0 980 900]);
set(gcf,'color','w');
hold on; axis off;

ax1 = axes('Position', position(1,:));
box on; hold on; grid off;
plot(xi, S_phi, 'k-');
ylabel('$S_{\mathrm{Phillips}}$','Interpreter','latex');
set(gca,'TickDir','out');
xlim([0 6]);
ax1.XMinorTick = 'on';
ax1.YMinorTick = 'on';
text(0.01,0.90,'$(a)$','Units','normalized','Interpreter','latex','FontSize',18);

ax2 = axes('Position', position(2,:));
box on; hold on; grid off;
plot(xi, S_howe_div, 'b-');
plot(xi, S_howe_acc, 'r--');
plot(xi, S_howe, 'k-.');
ylabel('$S_{\mathrm{Howe}}$','Interpreter','latex');
set(gca,'TickDir','out');
xlim([0 6]);
ax2.XMinorTick = 'on';
ax2.YMinorTick = 'on';
legend({'$\nabla\cdot(\omega\times u)$','$-c^{-2}a\cdot(\omega\times u)$','total'}, ...
    'Interpreter','latex','NumColumns',3,'Location','northoutside');
text(0.01,0.90,'$(b)$','Units','normalized','Interpreter','latex','FontSize',18);

ax3 = axes('Position', position(3,:));
box on; hold on; grid off;
plot(xi, S_dil_4DQdt, 'Color',[0.10 0.60 0.10]);
plot(xi, S_dil_detS, 'm--');
plot(xi, S_dil_stretch, 'c-.');
plot(xi, S_dil_acc, 'r-');
plot(xi, S_dil, 'k-');
xlabel('$r/\sigma$','Interpreter','latex');
ylabel('$S_{\mathrm{dilatation}}$','Interpreter','latex');
set(gca,'TickDir','out');
xlim([0 6]);
ax3.XMinorTick = 'on';
ax3.YMinorTick = 'on';
legend({'$4\mathrm{D}Q/\mathrm{D}t$','$-6\det(\mathbf{S})$', ...
        '$-(3/2)\,\omega\cdot\mathbf{S}\cdot\omega$', ...
        '$a\cdot(\nabla\times\omega)$','total'}, ...
       'Interpreter','latex','NumColumns',2,'Location','northoutside');
text(0.01,0.90,'$(c)$','Units','normalized','Interpreter','latex','FontSize',18);

exportgraphics(gcf, 'burgers_vortex_sources.pdf','ContentType','vector');
disp('Done: ./burgers_vortex_sources.pdf');
