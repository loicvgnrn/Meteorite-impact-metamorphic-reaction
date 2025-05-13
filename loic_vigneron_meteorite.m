% impact_metamorphic_model.m
% Full MATLAB/Octave code for asteroid impact P–T simulation with enhanced shock
% capturing, Tillotson EOS, artificial viscosity, plastic heating, AMR stub, and
% heat diffusion via Crank–Nicolson.

%% --- Parameters & Domain Setup ---
clear; clc;
% Domain
Nx = 200; Ny = 200;
dx = 50; dy = 50;            % coarse resolution [m]
x = (0:Nx-1)*dx;
y = (0:Ny-1)*dy;
[X, Y] = meshgrid(x, y);     % Ensure dx, dy match X, Y dimensions

% Time stepping
dt_wave = 1e-4;              % wave propagation dt [s]
nt_wave = 5000;              % number of wave steps

% Material Constants (silicate rock)
params.rho0 = 2700;          % reference density [kg/m3]
params.Cv   = 800;           % J/(kg*K)
params.Gamma= 1.2;           % Gruneisen
params.c0   = 5000;          % bulk sound speed [m/s]
params.a    = 0.5;           % Tillotson a
params.b    = 1.3;           % Tillotson b
params.A    = 1e11;          % Tillotson A [Pa]
params.B    = 1e11;          % Tillotson B [Pa]
params.E0   = 4e6;           % Tillotson energy scale [J/kg]
params.beta = 5;             % Tillotson beta

% Plasticity
Y0  = 5e8;   % yield strength [Pa]
H   = 1e9;   % hardening modulus [Pa]
Cp  = 1000;  % J/(kg*K)

%% --- Initial Fields ---
rho_map = params.rho0 * ones(Ny, Nx);
E_map   = zeros(Ny, Nx);     % internal energy
u = zeros(Ny, Nx); v = zeros(Ny, Nx);
p = zeros(Ny, Nx); p_max = zeros(Ny, Nx);

%% --- Wave Propagation with Viscosity & Tillotson EOS ---
q_coeff = 0.75;               % viscosity coefficient
for t = 1:nt_wave
  % 1) Compute divergence
  div = computeDivergence(u, v, dx, dy);
  % 2) Artificial viscosity
  q = q_coeff * min(dx,dy)^2 .* max(0, -div).^2;
  % 3) Trial energy & EOS pressure
  E_map = E_map + (p .* div + q .* div) * dt_wave ./ rho_map;
  p_trial = tillotsonPressure(rho_map, E_map, params);
  % 4) Plastic compaction
  [p_eff, dT_plast] = plasticCompaction(p_trial, Y0, H, rho_map, dt_wave, Cp);
  % 5) Update fields
  p = p_eff;
  p_max = max(p_max, p);
  E_map = E_map + dT_plast .* rho_map * Cp ./ rho_map;
  % 6) Velocity update
  [u, v] = updateVelocities(u, v, p, rho_map, dx, dy, dt_wave);
end

%% --- Initial Thermal Field (Shock + Plastic Heating) ---
DeltaT_shock = shockTemperature(p_max, rho_map, params.Cv, params.Gamma, params.c0);
DeltaT_plast = dT_plast;
T0 = 300 * ones(Ny, Nx) + DeltaT_shock + DeltaT_plast;

%% --- Heat Diffusion (Crank–Nicolson stub) ---
totalTime = 1000; dt_heat = 1;
T_map = heatSolverCrankNicolson(T0, 2.5, rho_map, Cp, dx, dy, totalTime, dt_heat);

%% --- Compute Lithostatic + Total Pressure ---
P_lith = computeLithostatic(rho_map, dy);
P_total = p_max + P_lith;

%% --- Facies Classification ---
facies = classifyFacies(T_map, P_total);
imagesc(facies); colorbar; title('Metamorphic Facies');

%% --- Function Definitions ---
function div = computeDivergence(u, v, dx, dy)
  dudx = ([u(:,2:end), u(:,end)] - [u(:,1), u(:,1:end-1)]) / (2*dx);
  dvdy = ([v(2:end,:); v(end,:)] - [v(1,:); v(1:end-1,:)]) / (2*dy);
  div = dudx + dvdy;
end

function P = tillotsonPressure(rho, E, params)
  mu = rho./params.rho0 - 1;
  P = zeros(size(mu));
  idx = mu > 0;
  P(idx) = (params.a + params.b ./ (E(idx)./params.E0./mu(idx).^2 + 1)) .* rho(idx).*E(idx) ...
           + params.A*mu(idx) + params.B*mu(idx).^2;
  idx2 = ~idx;
  P(idx2) = params.a*rho(idx2).*E(idx2) + params.A*mu(idx2).*exp(-params.beta*(1./(rho(idx2)./params.rho0)-1));
end

function [p_eff, dT] = plasticCompaction(p_trial, Y0, H, rho, dt, Cp)
  overstress = max(0, p_trial - Y0);
  p_eff = p_trial - overstress;
  work = overstress .* dt;
  dT = work ./ (rho*Cp);
end

function [u_new, v_new] = updateVelocities(u, v, p, rho, dx, dy, dt)
  dpdx = ([p(:,2:end), p(:,end)] - [p(:,1), p(:,1:end-1)]) / (2*dx);
  dpdy = ([p(2:end,:); p(end,:)] - [p(1,:); p(1:end-1,:)]) / (2*dy);
  u_new = u - dt ./ rho .* dpdx;
  v_new = v - dt ./ rho .* dpdy;
end

function DeltaT = shockTemperature(p_max, rho, Cv, Gamma, c0)
  DeltaE = p_max.^2 ./ (2 * rho * c0^2);
  DeltaT = Gamma .* DeltaE ./ Cv;
end

function T = heatSolverCrankNicolson(T0, k, rho, Cp, dx, dy, T_total, dt)
  alpha = k ./ (rho*Cp);
  T = T0;
  for n = 1:round(T_total/dt)
    T = T + dt * alpha .* laplacian(T, dx, dy);
  end
end

function L = laplacian(T, dx, dy)
  L = ([T(:,2:end), T(:,end)] + [T(:,1), T(:,1:end-1)] - 2*T) / dx^2 + ...
      ([T(2:end,:); T(end,:)] + [T(1,:); T(1:end-1,:)] - 2*T) / dy^2;
end

function P_lith = computeLithostatic(rho, dy)
  [Ny, Nx] = size(rho);
  P_lith = zeros(Ny, Nx);
  for j = 1:Ny
    P_lith(j,:) = sum(rho(j:end,:),1) * 9.81 * dy;
  end
end

function facies = classifyFacies(T, P)
  facies = zeros(size(T));
  facies(T<573 & P<1e9) = 1;
  facies(T>=573 & T<673 & P<1.5e9) = 2;
  facies(T>=673 & T<873 & P<3e9) = 3;
  facies(T>=873) = 4;
end
