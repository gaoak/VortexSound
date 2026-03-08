%% Burgers stretching shear layer: source-term distributions
% Output: ./shear_layer_sources.pdf
% Translation-velocity and translation-acceleration terms are shown separately.

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
Uinf  = 1.0;
sigma = 1.0;
k     = 0.5;
c0    = 10.0;

% background translation and acceleration
Ux = 0.30; Uy = 0.20; Uz = 0.00; %#ok<NASGU>
dUx = 0.06; dUy = 0.04; dUz = 0.00; %#ok<NASGU>

%% ------------------------ shear-layer fields --------------------------
y = linspace(-4*sigma, 4*sigma, 2001).';
eta = y/sigma;

v   = Uinf * erf(eta);
vp  = 2*Uinf/(sqrt(pi)*sigma) .* exp(-(y.^2)/(sigma^2));
vpp = -(2*y/sigma^2).*vp;

%% -------------------------- source terms ------------------------------
% Phillips equation RHS
S_phi = 2*k^2 * ones(size(y));

% Howe equation RHS decomposition
S_howe_div_base = -(vp.^2) - v.*vpp;
S_howe_div_U    = -Ux.*vpp;

S_howe_acc_base = (1/c0^2).*vp.*(k^2*y.*v - k^2*y.^2.*vp);
S_howe_acc_U    = (1/c0^2).*vp.*(Ux*dUy + k*y.*k*Ux + k*y.*Uy.*vp);
S_howe_acc_dU   = (1/c0^2).*vp.*(-dUx*Uy + k*y.*dUx + dUy.*v);

S_howe_base = S_howe_div_base + S_howe_acc_base;
S_howe_U    = S_howe_div_U + S_howe_acc_U;
S_howe_dU   = S_howe_acc_dU;
S_howe      = S_howe_base + S_howe_U + S_howe_dU;

% Dilatation equation RHS decomposition
S_dil_4DQdt   = zeros(size(y));
S_dil_detS    =  (3/2)*k*(vp.^2);
S_dil_stretch = -(3/2)*k*(vp.^2);
S_dil_acc_base= (k*y.*vp).*vpp;
S_dil_acc_dU  = -(dUx).*vpp;

S_dil_base = S_dil_4DQdt + S_dil_detS + S_dil_stretch + S_dil_acc_base;
S_dil_dU   = S_dil_acc_dU;
S_dil      = S_dil_base + S_dil_dU;

%% ----------------------------- plotting -------------------------------
position = [ ...
    0.10 0.78 0.85 0.18
    0.10 0.54 0.85 0.18
    0.10 0.30 0.85 0.18
    0.10 0.06 0.85 0.18];

handle = figure;
set(handle,'Position',[0 0 980 1050]);
set(gcf,'color','w');
hold on; axis off;

ax1 = axes('Position', position(1,:));
box on; hold on; grid off;
plot(eta, S_phi, 'k-');
ylabel('$S_{\mathrm{Phillips}}$','Interpreter','latex');
set(gca,'TickDir','out'); xlim([-4 4]);
ax1.XMinorTick = 'on'; ax1.YMinorTick = 'on';
text(0.01,0.86,'$(a)$','Units','normalized','Interpreter','latex','FontSize',18);

ax2 = axes('Position', position(2,:));
box on; hold on; grid off;
plot(eta, S_howe_base, 'k-');
plot(eta, S_howe_U, 'b--');
plot(eta, S_howe_dU, 'r-.');
plot(eta, S_howe, 'm-');
ylabel('$S_{\mathrm{Howe}}$','Interpreter','latex');
set(gca,'TickDir','out'); xlim([-4 4]);
ax2.XMinorTick = 'on'; ax2.YMinorTick = 'on';
legend({'base-flow terms','$U$-related terms','$\dot{U}$-related terms','total'}, ...
    'Interpreter','latex','NumColumns',2,'Location','northoutside');
text(0.01,0.86,'$(b)$','Units','normalized','Interpreter','latex','FontSize',18);

ax3 = axes('Position', position(3,:));
box on; hold on; grid off;
plot(eta, S_dil_base, 'k-');
plot(eta, S_dil_dU, 'r--');
plot(eta, S_dil, 'm-.');
ylabel('$S_{\mathrm{dilatation}}$','Interpreter','latex');
set(gca,'TickDir','out'); xlim([-4 4]);
ax3.XMinorTick = 'on'; ax3.YMinorTick = 'on';
legend({'base-flow terms','$\dot{U}$-related terms','total'}, ...
    'Interpreter','latex','NumColumns',3,'Location','northoutside');
text(0.01,0.86,'$(c)$','Units','normalized','Interpreter','latex','FontSize',18);

ax4 = axes('Position', position(4,:));
box on; hold on; grid off;
plot(eta, S_howe_div_U, 'b-');
plot(eta, S_howe_acc_U, 'b--');
plot(eta, S_howe_acc_dU, 'r-.');
plot(eta, S_dil_acc_dU, 'r-');
xlabel('$y/\sigma$','Interpreter','latex');
ylabel('translation terms','Interpreter','latex');
set(gca,'TickDir','out'); xlim([-4 4]);
ax4.XMinorTick = 'on'; ax4.YMinorTick = 'on';
legend({'Howe: $U$ in $\nabla\cdot(\omega\times u)$', ...
        'Howe: $U$ in acceleration term', ...
        'Howe: $\dot{U}$ in acceleration term', ...
        'Dilatation: $\dot{U}_x$ term'}, ...
       'Interpreter','latex','NumColumns',2,'Location','northoutside');
text(0.01,0.86,'$(d)$','Units','normalized','Interpreter','latex','FontSize',18);

exportgraphics(gcf, 'shear_layer_sources.pdf','ContentType','vector');
disp('Done: ./shear_layer_sources.pdf');
