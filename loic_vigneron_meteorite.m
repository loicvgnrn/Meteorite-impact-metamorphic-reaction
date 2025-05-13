% Enhanced Meteorite Impact Simulation with Higher-Order Accuracy, Heat Conduction, and Metamorphic Mapping
% Author: ChatGPT (Refactor)

function enhanced_meteorite_simulation()
    %% === User Input and Validation ===
    params = getParameters();
    validateParameters(params);

    %% === Grid Setup ===
    [x, y, dx, dy, Nx, Ny, ix_impact, jy_impact] = setupGrid(params);

    %% === Geological Property Mapping (Smooth) ===
    [rho_map, beta_map, k_map] = smoothPropertyMap(params, y, Ny, Nx);

    %% === Wave Propagation Setup (4th-Order FD) ===
    [dt, Nt, cmax] = computeTimeStep(dx, dy, rho_map, beta_map);
    [u_prev, u_now] = deal(zeros(Ny, Nx));
    source = defineRickerSource(params, Nt, dt);

    %% === Pressure and Visualization Setup ===
    p_max = zeros(Ny, Nx);
    [fig, ax1, ax2] = initLivePlot(x, y);
    hWait = waitbar(0, 'Running Wave Propagation...');

    %% === Wave Propagation Loop ===
    for n = 1:Nt
        t = n * dt;
        u_now = injectGaussianSource(u_now, source(n), params, ix_impact, jy_impact, dx, dy);
        u_next = waveStep4thOrder(u_now, u_prev, dt, dx, dy, rho_map, beta_map);
        u_prev = u_now; u_now = u_next;

        p = -1 ./ beta_map .* divergence(u_now, dy, dx);
        p_max = max(p_max, p);

        if mod(n, max(1, round(Nt/50))) == 0
            updateLivePlot(ax1, ax2, x, y, u_now, p, t);
        end

        waitbar(n/Nt, hWait, sprintf('Wave Propagation: %.1f%% (%.2f s)', 100*n/Nt, params.Tmax - t));
    end

    %% === Heat Conduction (Implicit Crank-Nicolson) ===
    waitbar(0, hWait, 'Running Heat Conduction...');
    T0 = initialTemperature(x, y, params, ix_impact, jy_impact, dx, dy);
   T_map = heatSolverExplicit(T0, k_map, dx, dy, 1000, 1);


    %% === Total Pressure ===
    P_lith = computeLithostatic(rho_map, dy);
    P_total = p_max + P_lith;

    %% === Facies Classification ===
    waitbar(1, hWait, 'Classifying Metamorphic Facies...');
    classifyFacies(x, y, T_map, P_total);
    close(hWait);
end

%% === Supporting Functions ===

function params = getParameters()
    params.meteorite_radius = 100;
    params.impact_velocity = 2000;
    params.simulation_width = 2000;
    params.simulation_depth = 6000;
    params.Nx = 300; params.Ny = 300;
    params.Tmax = 0.5;
    params.impact_x = 1000; params.impact_y = 500;
    params.f0 = 100;
    % Layers: [thicknesses; rho; beta; k]
    params.layers = [1000,2000,3000; 2500,2700,3000; 5e-10,4e-10,3e-10; 2.5,3.0,3.5];
end

function validateParameters(p)
    assert(all(p.layers(:) > 0), 'Layer parameters must be positive.');
end

function [x, y, dx, dy, Nx, Ny, ix, jy] = setupGrid(p)
    Nx = p.Nx; Ny = p.Ny;
    dx = p.simulation_width / (Nx-1);
    dy = p.simulation_depth / (Ny-1);
    x = linspace(0, p.simulation_width, Nx);
    y = linspace(0, p.simulation_depth, Ny);
    ix = round(p.impact_x / dx) + 1;
    jy = round(p.impact_y / dy) + 1;
end

function [rho_map, beta_map, k_map] = smoothPropertyMap(p, y, Ny, Nx)
    rho_map = zeros(Ny,Nx); beta_map = zeros(Ny,Nx); k_map = zeros(Ny,Nx);
    bounds = [0 cumsum(p.layers(1,:))];
    for j = 1:Ny
        depth = y(j);
        idx = find(depth < bounds,1) - 1; if isempty(idx), idx = size(p.layers,2); end
        rho_map(j,:) = p.layers(2,idx);
        beta_map(j,:) = p.layers(3,idx);
        k_map(j,:) = p.layers(4,idx);
    end
end

function [dt, Nt, cmax] = computeTimeStep(dx, dy, rho, beta)
    c = sqrt(1 ./ (rho .* beta));
    cmax = max(c(:));
    dt = 0.4 * min(dx, dy) / (sqrt(2) * cmax);
    Nt = floor(0.5 / dt);
end

function source = defineRickerSource(p, Nt, dt)
    f0 = p.f0; t0 = 3 / f0;
    t = (1:Nt)*dt;
    source = (1 - 2*(pi*f0)^2*(t-t0).^2).*exp(-(pi*f0)^2*(t-t0).^2);
    KE = 0.5 * 3000 * (4/3*pi*p.meteorite_radius^3) * p.impact_velocity^2;
    source = 1e-12 * KE * source;
end

function u = injectGaussianSource(u, val, p, ix, jy, dx, dy)
    [Ny, Nx] = size(u);
    [X, Y] = meshgrid((0:Nx-1)*dx, (0:Ny-1)*dy);
    cx = (ix-1)*dx; cy = (jy-1)*dy;
    r2 = (X-cx).^2 + (Y-cy).^2;
    sigma2 = (p.meteorite_radius/2)^2;
    gauss = exp(-r2/(2*sigma2));
    gauss = gauss / sum(gauss(:));
    u = u + val * gauss;
end

function u_next = waveStep4thOrder(u, u_prev, dt, dx, dy, rho, beta)
    c2 = 1./(rho .* beta);
    kernel_x = [-1/12, 4/3, -5/2, 4/3, -1/12]/dx^2;
    kernel_y = kernel_x';
    u_xx = conv2(u, kernel_x, 'same');
    u_yy = conv2(u, kernel_y, 'same');
    u_next = 2*u - u_prev + dt^2 * c2 .* (u_xx + u_yy);
end

function [fig, ax1, ax2] = initLivePlot(x, y)
    fig = figure('Name','Live Simulation','Position',[100 100 1200 500]);
    ax1 = subplot(1,2,1); imagesc(x,y,zeros(length(y),length(x))); axis xy; colorbar; title('Displacement');
    ax2 = subplot(1,2,2); imagesc(x,y,zeros(length(y),length(x))); axis xy; colorbar; title('Pressure');
end

function updateLivePlot(ax1, ax2, x, y, u, p, t)
    imagesc(ax1, x, y, u); title(ax1,sprintf('u, t=%.3f',t)); axis(ax1,'xy'); colorbar(ax1);
    imagesc(ax2, x, y, p); title(ax2,'Pressure'); axis(ax2,'xy'); colorbar(ax2); drawnow;
end

function div_u = divergence(u, dy, dx)
    [dudy, dudx] = gradient(u, dy, dx);
    div_u = dudx + dudy;
end

function T0 = initialTemperature(x, y, p, ix, jy, dx, dy)
    [X,Y] = meshgrid(x,y);
    T0 = 15 + 30*(Y/1000);
    cx = (ix-1)*dx; cy = (jy-1)*dy;
    r2 = (X-cx).^2 + (Y-cy).^2;
    T0 = T0 + 1000*exp(-r2/(2*p.meteorite_radius^2));
end

function T = heatSolverExplicit(T0, k_map, dx, dy, T_total, dt_heat)
    [Ny,Nx] = size(T0);
    T = T0;
    alpha = mean(k_map(:)) / (2700 * 1000);
    Nt = round(T_total / dt_heat);
    for n = 1:Nt
        Txx = ([T(:,2:end), T(:,end)] - 2*T + [T(:,1), T(:,1:end-1)]) / dx^2;
        Tyy = ([T(2:end,:); T(end,:)] - 2*T + [T(1,:); T(1:end-1,:)]) / dy^2;
        T = T + dt_heat * alpha * (Txx + Tyy);
    end
end


function P = computeLithostatic(rho, dy)
    g = 9.81; [Ny,Nx] = size(rho);
    P = zeros(Ny,Nx);
    for j = 2:Ny
        P(j,:) = P(j-1,:) + 0.5*(rho(j,:) + rho(j-1,:))*g*dy;
    end
end

function classifyFacies(x, y, T, P)
    pkbar = P / 1e8;
    facies = {'Zeolite','Greenschist','Amphibolite','Granulite','Eclogite'};
    Pr = [0 3; 2 6; 4 9; 7 12; 10 20];
    Tr = [100 300; 200 500; 500 750; 700 1000; 600 1000];
    F = zeros(size(T));
    for i = 1:length(facies)
        mask = pkbar >= Pr(i,1) & pkbar <= Pr(i,2) & T >= Tr(i,1) & T <= Tr(i,2);
        F(mask) = i;
    end
    figure; imagesc(x, y, F); axis xy;
    colormap(lines(length(facies))); colorbar('Ticks',1:length(facies),'TickLabels',facies);
    title('Metamorphic Facies'); xlabel('x (m)'); ylabel('y (m)');
end
