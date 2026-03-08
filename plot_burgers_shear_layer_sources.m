%% Burgers stretching shear layer: source-term distributions
% Output: figures/shear_layer_sources.pdf

clear; clc; close all;

%% --------------------------- plot style -------------------------------
set(0,'defaultlinelinewidth',2);
set(0,'defaultaxeslinewidth',2);
set(0,'defaultaxesfontsize',16);
set(0,'defaultTextFontName', 'Times New Roman');
set(0,'defaultAxesFontName', 'Times New Roman');
set(0,'defaulttextfontsize',16);
set(0,'DefaultLineMarkerSize',9);

outDir = 'figures';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% --------------------------- parameters -------------------------------
Uinf  = 1.0;
sigma = 1.0;
k     = 0.5;
c0    = 10.0;

% background translation and acceleration
Ux = 0.0; Uy = 0.0; Uz = 0.0; %#ok<NASGU>
dUx = 0.0; dUy = 0.0; dUz = 0.0; %#ok<NASGU>

%% ------------------------ shear-layer fields --------------------------
y = linspace(-4*sigma, 4*sigma, 2001).';
eta = y/sigma;

v   = Uinf * erf(eta);
vp  = 2*Uinf/(sqrt(pi)*sigma) .* exp(-(y.^2)/(sigma^2));
vpp = -(2*y/sigma^2).*vp;

%% -------------------------- source terms ------------------------------
% Phillips equation RHS
S_phi = 2*k^2 * ones(size(y));

% Howe equation RHS
S_howe_div = -(vp.^2) - (Ux + v).*vpp;
S_howe_acc = (1/c0^2) .* vp .* ( ...
      Ux*dUy - dUx*Uy ...
    + k*y.*(dUx + k*Ux) ...
    + (dUy + k^2*y).*v ...
    + k*y.*(Uy - k*y).*vp );
S_howe = S_howe_div + S_howe_acc;

% Dilatation equation RHS
S_dil_4DQdt   = zeros(size(y));
S_dil_detS    =  (3/2)*k*(vp.^2);
S_dil_stretch = -(3/2)*k*(vp.^2);
S_dil_acc     = -(dUx - k*y.*vp).*vpp;
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
plot(eta, S_phi, 'k-');
ylabel('$S_{\mathrm{Phillips}}$','Interpreter','latex');
set(gca,'TickDir','out');
xlim([-4 4]);
ax1.XMinorTick = 'on';
ax1.YMinorTick = 'on';
text(0.01,0.90,'$(a)$','Units','normalized','Interpreter','latex','FontSize',18);

ax2 = axes('Position', position(2,:));
box on; hold on; grid off;
plot(eta, S_howe_div, 'b-');
plot(eta, S_howe_acc, 'r--');
plot(eta, S_howe, 'k-.');
ylabel('$S_{\mathrm{Howe}}$','Interpreter','latex');
set(gca,'TickDir','out');
xlim([-4 4]);
ax2.XMinorTick = 'on';
ax2.YMinorTick = 'on';
legend({'$\nabla\cdot(\omega\times u)$','$-c^{-2}a\cdot(\omega\times u)$','total'}, ...
    'Interpreter','latex','NumColumns',3,'Location','northoutside');
text(0.01,0.90,'$(b)$','Units','normalized','Interpreter','latex','FontSize',18);

ax3 = axes('Position', position(3,:));
box on; hold on; grid off;
plot(eta, S_dil_4DQdt, 'Color',[0.10 0.60 0.10]);
plot(eta, S_dil_detS, 'm--');
plot(eta, S_dil_stretch, 'c-.');
plot(eta, S_dil_acc, 'r-');
plot(eta, S_dil, 'k-');
xlabel('$y/\sigma$','Interpreter','latex');
ylabel('$S_{\mathrm{dilatation}}$','Interpreter','latex');
set(gca,'TickDir','out');
xlim([-4 4]);
ax3.XMinorTick = 'on';
ax3.YMinorTick = 'on';
legend({'$4\mathrm{D}Q/\mathrm{D}t$','$-6\det(\mathbf{S})$', ...
        '$-(3/2)\,\omega\cdot\mathbf{S}\cdot\omega$', ...
        '$a\cdot(\nabla\times\omega)$','total'}, ...
       'Interpreter','latex','NumColumns',2,'Location','northoutside');
text(0.01,0.90,'$(c)$','Units','normalized','Interpreter','latex','FontSize',18);

exportgraphics(gcf, fullfile(outDir,'shear_layer_sources.pdf'),'ContentType','vector');
disp('Done: figures/shear_layer_sources.pdf');
