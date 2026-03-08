%% Burgers vortex: source-term distributions (numerical differentiation)
% Output: ./burgers_vortex_sources.pdf
% Notes:
% 1) Spatial derivatives are evaluated by 6th-order central finite difference
%    on a dense uniform grid.
% 2) Terms associated with translation velocity U and translation acceleration
%    dU are separated and plotted explicitly.

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

% Translation velocity and acceleration (local cylindrical components)
Ur  = 0.20;  Ut  = 0.15;  Uz  = 0.00; %#ok<NASGU>
dUr = 0.05;  dUt = 0.04;  dUz = 0.00; %#ok<NASGU>

%% ------------------------- vortex fields ------------------------------
% dense uniform grid for high-order finite difference
% r = 0 is singular in cylindrical coordinates; start from a small positive
% radius r0 = eps_r*sigma instead of 0.
eps_r = 1e-3;
r0 = eps_r*sigma;
rMax = 6.0*sigma;
N = 6001;
r = linspace(r0, rMax, N).';
xi = r/sigma;

% Optional sanity check
if r(1) <= 0
    error('Radial grid must start from a positive value to avoid r=0 singularity.');
end

% Base velocity components (v_r, v_theta, v_z) in cylindrical coordinates
vr = -0.5*k*r;
vt = Gamma./(2*pi*r) .* (1 - exp(-(r.^2)/(sigma^2)));
vz = zeros(size(r));

% Numerical derivatives (6th-order finite difference)
dvtdr = d1_6th_uniform(vt, r);

% Vorticity: omega_z = (1/r)d(r v_theta)/dr = dv_theta/dr + v_theta/r
omega = dvtdr + vt./r;
domegadr = d1_6th_uniform(omega, r);

% Curl of vorticity for axisymmetric omega_z(r):
% (curl omega)_theta = - d(omega_z)/dr
curlw_t = -domegadr;

% Material acceleration of base flow Dv/Dt
ar = (k^2*r/4) - vt.^2./r;
at = -(k/2).*(vt + r.*dvtdr);
az = zeros(size(r));

%% -------------------------- source terms ------------------------------
% Phillips equation RHS
S_phi = (3/2)*k^2 - 2*(vt./r).*dvtdr;

% ===== Howe equation RHS decomposition =====
% 1) divergence part: div(omega x u) = div(omega x v) + U.(curl omega)
S_howe_div_base = -(1./r).*d1_6th_uniform(r.*omega.*vt, r);      % intrinsic
S_howe_div_U    = Ut .* curlw_t;                                 % U-related

% 2) acceleration part: -c^{-2} a.(omega x u)
% Following manuscript decomposition:
%   + c^{-2} v.(omega x Dv/Dt)
%   + c^{-2} U.(omega x Dv/Dt)
%   - c^{-2} dU.(omega x v)
%   + (omega/c^2) (dU x U).e_z
omega_cross_Dv_r = -omega .* at;
omega_cross_Dv_t =  omega .* ar;

omega_cross_v_r = -omega .* vt;
omega_cross_v_t =  omega .* vr;

S_howe_acc_base = (1/c0^2) .* (vr.*omega_cross_Dv_r + vt.*omega_cross_Dv_t);
S_howe_acc_U    = (1/c0^2) .* (Ur.*omega_cross_Dv_r + Ut.*omega_cross_Dv_t);
S_howe_acc_dU   = -(1/c0^2) .* (dUr.*omega_cross_v_r + dUt.*omega_cross_v_t) ...
                  + (omega/c0^2) .* (dUr*Ut - dUt*Ur);

S_howe_base = S_howe_div_base + S_howe_acc_base;
S_howe_U    = S_howe_div_U + S_howe_acc_U;
S_howe_dU   = S_howe_acc_dU;
S_howe      = S_howe_base + S_howe_U + S_howe_dU;

% ===== Dilatation equation RHS decomposition =====
trA2 = S_phi;
Q = -0.5*trA2;
dQdr = d1_6th_uniform(Q, r);
S_dil_4DQdt = 4*vr.*dQdr;

Srth = 0.5*(dvtdr - vt./r);
detS = k*(k^2/4 - Srth.^2);
S_dil_detS = -6*detS;
S_dil_stretch = -(3/2)*k*(omega.^2);

S_dil_acc_base = at .* curlw_t;
S_dil_acc_dU   = dUt .* curlw_t;

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
plot(xi, S_phi, 'k-');
ylabel('$S_{\mathrm{Phillips}}$','Interpreter','latex');
set(gca,'TickDir','out');
xlim([0 6]);
ax1.XMinorTick = 'on'; ax1.YMinorTick = 'on';
text(0.01,0.86,'$(a)$','Units','normalized','Interpreter','latex','FontSize',18);

ax2 = axes('Position', position(2,:));
box on; hold on; grid off;
plot(xi, S_howe_base, 'k-');
plot(xi, S_howe_U, 'b--');
plot(xi, S_howe_dU, 'r-.');
plot(xi, S_howe, 'm-');
ylabel('$S_{\mathrm{Howe}}$','Interpreter','latex');
set(gca,'TickDir','out');
xlim([0 6]);
ax2.XMinorTick = 'on'; ax2.YMinorTick = 'on';
legend({'base-flow terms','$U$-related terms','$\dot{U}$-related terms','total'}, ...
    'Interpreter','latex','NumColumns',2,'Location','northoutside');
text(0.01,0.86,'$(b)$','Units','normalized','Interpreter','latex','FontSize',18);

ax3 = axes('Position', position(3,:));
box on; hold on; grid off;
plot(xi, S_dil_base, 'k-');
plot(xi, S_dil_dU, 'r--');
plot(xi, S_dil, 'm-.');
ylabel('$S_{\mathrm{dilatation}}$','Interpreter','latex');
set(gca,'TickDir','out');
xlim([0 6]);
ax3.XMinorTick = 'on'; ax3.YMinorTick = 'on';
legend({'base-flow terms','$\dot{U}$-related terms','total'}, ...
    'Interpreter','latex','NumColumns',3,'Location','northoutside');
text(0.01,0.86,'$(c)$','Units','normalized','Interpreter','latex','FontSize',18);

ax4 = axes('Position', position(4,:));
box on; hold on; grid off;
plot(xi, S_howe_div_U, 'b-');
plot(xi, S_howe_acc_U, 'b--');
plot(xi, S_howe_acc_dU, 'r-.');
plot(xi, S_dil_acc_dU, 'r-');
xlabel('$r/\sigma$','Interpreter','latex');
ylabel('translation terms','Interpreter','latex');
set(gca,'TickDir','out');
xlim([0 6]);
ax4.XMinorTick = 'on'; ax4.YMinorTick = 'on';
legend({'Howe: $U\cdot(\nabla\times\omega)$', ...
        'Howe: $c^{-2}U\cdot(\omega\times Dv/Dt)$', ...
        'Howe: $\dot{U}$ terms', ...
        'Dilatation: $\dot{U}_\theta(\nabla\times\omega)_\theta$'}, ...
       'Interpreter','latex','NumColumns',2,'Location','northoutside');
text(0.01,0.86,'$(d)$','Units','normalized','Interpreter','latex','FontSize',18);

exportgraphics(gcf, 'burgers_vortex_sources.pdf','ContentType','vector');
disp('Done: ./burgers_vortex_sources.pdf');

%% --------------------------- local function ---------------------------
function dfdx = d1_6th_uniform(f, x)
% First derivative on uniform grid:
% interior: 6th-order central; near boundaries: 2nd-order one-sided/centered.

n = numel(f);
h = x(2)-x(1);
dfdx = zeros(size(f));

% interior points (6th-order central)
for i = 4:n-3
    dfdx(i) = (f(i-3) - 9*f(i-2) + 45*f(i-1) - 45*f(i+1) + 9*f(i+2) - f(i+3)) / (60*h);
end

% boundaries (2nd-order formulas)
dfdx(1)   = (-3*f(1) + 4*f(2) - f(3)) / (2*h);
dfdx(2)   = (f(3) - f(1)) / (2*h);
dfdx(3)   = (f(4) - f(2)) / (2*h);
dfdx(n)   = (3*f(n) - 4*f(n-1) + f(n-2)) / (2*h);
dfdx(n-1) = (f(n) - f(n-2)) / (2*h);
dfdx(n-2) = (f(n-1) - f(n-3)) / (2*h);
end
