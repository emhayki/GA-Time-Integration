
clc; clear; close all;

% -------------------------------------------------
% 2D Taylor-Green vortex with P2/P1 finite elements
% -------------------------------------------------

cfg.L = 2 * pi;
cfg.meshSize = 0.225;

cfg.nu = 0.1;

cfg.finalTime = 3;
cfg.dt = 0.01;
cfg.rhoInfinity = 0;

cfg.gradDiv = 0.03;

cfg.newtonTolerance = 1e-9;
cfg.maxNewtonIterations = 22;

cfg.velocityDecay = -2 * cfg.nu;
cfg.omega = 2 * pi;

cfg.pressureNode = 1;
cfg.pressureValue = 0;

cfg.outputFolder = 'OscillatingTaylorGreen';

cfg.saveCombinedVideo = true;
cfg.videoFrameRate = 24;
cfg.tracerVideoSeconds = 12;

cfg.particlesPerSide = 20;
cfg.tracerGridSize = 220;
cfg.tracerTrailFrames = 220;
cfg.movieMethod = 'GA-234';

if (isfolder(cfg.outputFolder) == false)
    mkdir(cfg.outputFolder);
end

% ------------------------------------------
% Finite Element Mesh and Initial Conditions
% ------------------------------------------

fprintf('\n-------------------------------------------------------');
fprintf('\nTaylor-Green vortex | P2/P1 finite element simulation\n');
fprintf('-------------------------------------------------------\n');

fem = buildMesh(cfg);
fem = setupBoundaryDOFs(fem, cfg);

xu = fem.velocityNodes(:,1);
yu = fem.velocityNodes(:,2);

initialU =  sin(xu) .* cos(yu);
initialV = -cos(xu) .* sin(yu);

xp = fem.nodes(:,1);
yp = fem.nodes(:,2);

initialP = 0.25 * (cos(2 * xp) + cos(2 * yp));
initialP = initialP - initialP(fem.pinnedPressure) + cfg.pressureValue;

fem.initialState.u = initialU;
fem.initialState.v = initialV;
fem.initialState.p = initialP;

d = real((cfg.velocityDecay + 1i * cfg.omega).^(0:3));

fem.initialState.du = d(2) * initialU;
fem.initialState.dv = d(2) * initialV;

fem.initialState.ddu = d(3) * initialU;
fem.initialState.ddv = d(3) * initialV;

fem.initialState.dddu = d(4) * initialU;
fem.initialState.dddv = d(4) * initialV;

fem.boundaryU0 = initialU(fem.fixedVelocity);
fem.boundaryV0 = initialV(fem.fixedVelocity);

fem = assembleFEMMatrices(fem);

fem.forceU = fem.M * initialU;
fem.forceV = fem.M * initialV;

% ------------------------------------------
% Method Comparison
% ------------------------------------------
methods   = {'GA-2', 'GA-23', 'GA-234'};
solutions = cell(numel(methods), 1);

errors     = zeros(numel(methods), 1);
movieIndex = find(strcmp(methods, cfg.movieMethod), 1);

for k = 1:numel(methods)

    fprintf('\nMethod: %s (rho Infinity = %.2f)| dt = %.4f\n', methods{k}, cfg.rhoInfinity, cfg.dt);

    if(k == movieIndex)
        snapshotTimes = (0:cfg.dt:cfg.finalTime)';
    else
        snapshotTimes = [];
    end

    solutions{k} = solveFlow(fem, cfg, methods{k}, cfg.dt, snapshotTimes);
    errors(k)    = sqrt(mean((solutions{k}.energy - exactEnergyAt(solutions{k}.time, cfg)).^2));

    fprintf('  Energy RMS vs exact: %.6e\n', errors(k));
end
 
% ------------------------------------------
% Plotting
% ------------------------------------------

for k = 1:numel(methods)
    plotEnergyComparison(cfg, solutions{k}, methods{k});
    plotVelocityComparison(fem, cfg, solutions{k}, methods{k}, k == 1);
end

if (cfg.saveCombinedVideo)
   createFlowAnimation(fem, cfg, solutions{movieIndex});
end


fprintf('\nSimulation completed successfully. Results saved to: %s\n', cfg.outputFolder);

% ============================================================
% FINITE ELEMENT SETUP
% ============================================================

function fem = buildMesh(cfg)

    n = max(4, round(cfg.L/cfg.meshSize));
    x = linspace(0, cfg.L, n+1);

    [X, Y] = meshgrid(x, x);
    nodes = [X(:), Y(:)];

    triangles = zeros(2 * n * n, 3);
    row = 0;

    for i = 1:n
        for j = 1:n

            a = (i-1) * (n+1) + j;
            b = i * (n+1) + j;

            row = row + 1;
            triangles(row,:) = [a, b, b + 1];

            row = row + 1;
            triangles(row,:) = [a, b + 1, a + 1];
        end
    end


    nPressure = size(nodes, 1);
    nElements = size(triangles, 1);

    edgeList = sort([triangles(:,[1 2]); triangles(:,[2 3]); triangles(:,[3 1])], 2);
    [edges, ~, edgeIndex] = unique(edgeList, 'rows');

    midpoints = (nodes(edges(:,1),:) + nodes(edges(:,2),:))/2;
    velocityNodes = [nodes; midpoints];

    velocityTriangles = [triangles, ...
                         nPressure + edgeIndex(1:nElements), ...
                         nPressure + edgeIndex(nElements + 1:2 * nElements), ...
                         nPressure + edgeIndex(2 * nElements + 1:end)];

    fem.nodes = nodes;
    fem.edges = edges; 
    fem.triangles = triangles;

    fem.velocityNodes = velocityNodes;
    fem.velocityTriangles = velocityTriangles;

    fem.nPressure = nPressure;
    fem.nVelocity = size(velocityNodes, 1);

    fem.area = cfg.L * cfg.L;
    fem.nElements = nElements;
end

% -------
function fem = setupBoundaryDOFs(fem, cfg)

    boundaryEdges = freeBoundary(triangulation(fem.triangles, fem.nodes));
    [~, boundaryEdgeIndex] = ismember(sort(boundaryEdges,2), fem.edges, 'rows');

    fem.fixedVelocity = unique([boundaryEdges(:); fem.nPressure+boundaryEdgeIndex]);
    fem.pinnedPressure = cfg.pressureNode;

    fixed = [fem.fixedVelocity; ...
             fem.nVelocity + fem.fixedVelocity; ...
             2 * fem.nVelocity + fem.pinnedPressure];

    fem.freeDofs = setdiff((1:2 * fem.nVelocity + fem.nPressure)', fixed);
end

% -------
function fem = assembleFEMMatrices(fem)

    triangles = fem.triangles;
    velocityTriangles = fem.velocityTriangles;

    nodes = fem.nodes;
    nElements = fem.nElements;

    fem.elements = repmat(struct('nodes', [], 'area', [], 'gradL', []), nElements, 1);

    nu = fem.nVelocity;
    np = fem.nPressure;

    fem.M = sparse(nu,nu);
    fem.K = sparse(nu,nu);

    fem.Bx = sparse(np,nu);
    fem.By = sparse(np,nu);

    fem.Gxx = sparse(nu,nu);
    fem.Gxy = sparse(nu,nu);
    fem.Gyy = sparse(nu,nu);

    [quadPoints, quadWeights] = triangleQuadrature();

    for e = 1:nElements
        pNodes = triangles(e,:);
        uNodes = velocityTriangles(e,:);

        [area, gradL] = triangleGeometry(nodes(pNodes,:));
        fem.elements(e) = struct('nodes', uNodes, 'area', area, 'gradL', gradL);

        Me = zeros(6);
        Ke = zeros(6);

        Bxe = zeros(3,6);
        Bye = zeros(3,6);

        Gxxe = zeros(6);
        Gxye = zeros(6);
        Gyye = zeros(6);

        for q = 1:numel(quadWeights)
            
            L = quadPoints(q,:);
            [N, Nx, Ny] = quadraticShape(L, gradL);

            weight = area * quadWeights(q);

            Me = Me + weight * (N * N');
            Ke = Ke + weight * (Nx * Nx' + Ny * Ny');

            Bxe = Bxe + weight * (L' * Nx');
            Bye = Bye + weight * (L' * Ny');

            Gxxe = Gxxe + weight * (Nx * Nx');
            Gxye = Gxye + weight * (Nx * Ny');
            Gyye = Gyye + weight * (Ny * Ny');
        end

        fem.M(uNodes,uNodes) = fem.M(uNodes,uNodes) + Me;
        fem.K(uNodes,uNodes) = fem.K(uNodes,uNodes) + Ke;

        fem.Bx(pNodes,uNodes) = fem.Bx(pNodes,uNodes) + Bxe;
        fem.By(pNodes,uNodes) = fem.By(pNodes,uNodes) + Bye;

        fem.Gxx(uNodes,uNodes) = fem.Gxx(uNodes,uNodes) + Gxxe;
        fem.Gxy(uNodes,uNodes) = fem.Gxy(uNodes,uNodes) + Gxye;
        fem.Gyy(uNodes,uNodes) = fem.Gyy(uNodes,uNodes) + Gyye;
    end

    fprintf('\nFinite element assembly completed\n');
    fprintf('  Total degrees of freedom = %d\n', 2*nu + np);
    fprintf('  Total mesh elements      = %d\n', nElements);
    fprintf('  Total mesh nodes         = %d\n', nu);
    fprintf('  Velocity nodes (P2)      = %d\n', nu);
    fprintf('  Pressure nodes (P1)      = %d\n', np);
end

% -------
function [area, gradL] = triangleGeometry(P)

    J = [(P(2,:)-P(1,:))', (P(3,:)-P(1,:))'];
    area  = abs(det(J))/2;
    gradL = [P(2,2) - P(3,2), P(3,2) - P(1,2), P(1,2) - P(2,2); ...
             P(3,1) - P(2,1), P(1,1) - P(3,1), P(2,1) - P(1,1)]/det(J);
end

% -------
function [N, Nx, Ny] = quadraticShape(L, gradL)

    a = L(1);
    b = L(2);
    c = L(3);

    N = [a*(2*a-1); b*(2*b-1); c*(2*c-1); 4*a*b; 4*b*c; 4*c*a];

    dNdL = [4*a-1, 0, 0; 0, 4*b-1, 0; 0, 0, 4*c-1; ...
            4*b, 4*a, 0; 0, 4*c, 4*b; 4*c, 0, 4*a];

    derivatives = dNdL * gradL';

    Nx = derivatives(:,1);
    Ny = derivatives(:,2);
end

% -------
function [points, weights] = triangleQuadrature()
    points = [            1/3,             1/3,             1/3; ...
              .05971587178977, .47014206410511, .47014206410511; ...
              .47014206410511, .05971587178977, .47014206410511; ...
              .47014206410511, .47014206410511, .05971587178977; ...
              .79742698535309, .10128650732346, .10128650732346; ...
              .10128650732346, .79742698535309, .10128650732346; ...
              .10128650732346, .10128650732346, .79742698535309];

    weights = [.225; repmat(.13239415278851,3,1); repmat(.12593918054483,3,1)];
end


% ============================================================
% TIME INTEGRATION
% ============================================================

function result = solveFlow(fem, cfg, method, dt, snapshotTimes)

    steps = round(cfg.finalTime/dt);
    coeff = timeCoefficients(method, cfg.rhoInfinity);
    state = fem.initialState;

    result.time = (0:steps)' * dt;

    result.energy = zeros(steps+1, 1);
    result.energy(1) = kineticEnergy(fem, state);

    result.snapTimes = snapshotTimes(:);
    result.uSnap = zeros(fem.nVelocity, numel(snapshotTimes));
    result.vSnap = zeros(fem.nVelocity, numel(snapshotTimes));
    result.newtonIterations = 0;

    if (isempty(snapshotTimes) == false)
        result.uSnap(:,1) = state.u;
        result.vSnap(:,1) = state.v;
    end

    nextSnapshot = 2;

    for step = 1:steps

        previous    = state;
        currentTime = step * dt;

        [candidate, ok, iters, residual] = newtonStep(fem, cfg, coeff, dt, previous, currentTime);
        result.newtonIterations = result.newtonIterations + iters;

        if (ok == false)
            error('%s: Newton solver failed at t = %.4f (residual = %.3e). Try reducing cfg.dt.', method, currentTime, residual);
        end

        state = candidate;
        result.energy(step+1) = kineticEnergy(fem, state);

        while(nextSnapshot <= numel(snapshotTimes) && abs(currentTime-snapshotTimes(nextSnapshot)) < 1e-8)
            result.uSnap(:,nextSnapshot) = state.u;
            result.vSnap(:,nextSnapshot) = state.v;
            nextSnapshot = nextSnapshot + 1;
        end

        if (mod(step, max(1, round(steps/8))) == 0 || step == steps)
            fprintf('  %-6s | Step %3d/%3d | Time: %5.2f s | Newton iterations: %2d | Residual: %9.2e\n', method, step, steps, currentTime, iters, residual);
        end
    end

    result.u = state.u;
    result.v = state.v;
    result.p = state.p;

    divergence = fem.Bx * state.u + fem.By * state.v;
    divergence(fem.pinnedPressure) = [];
    result.divergence = norm(divergence)/max(norm([state.u;state.v]), eps);

    if(result.divergence > 1e-5)
        warning('%s continuity residual = %.3g.', method, result.divergence);
    end
end

% -------
function coeff = timeCoefficients(method, rhoInfinity)

    r = rhoInfinity;
    coeff.alpha = 1/(1 + r);
    coeff.gamma = coeff.alpha;

    switch method
        case 'GA-2'
            beta0 = (3 - r)/(2*(1 + r));
            beta2 = 0;
            beta3 = 0;

        case 'GA-23'
            beta0 = (10 - 5*r + r^2)/(6*(1 + r));
            beta2 = -(1 - r)^2/(6*(1 + r));
            beta3 = 0;

        case 'GA-234'
            beta0 = (35 - 21*r + 7*r^2 - r^3)/(20*(1 + r));
            beta2 = -((1 - r)^2 * (5 - r))/(20*(1 + r));
            beta3 = -(1 - r)^3/(20 * (1 + r)^2);

        otherwise
            error('Unknown time-integration method: %s', method);
    end

    coeff.beta = [beta0, 1-beta0, beta2, beta3];
end

% -------
function a = amplitudeExact(t, cfg)
    a = exp(cfg.velocityDecay * t) .* cos(cfg.omega * t);
end

% -------
function gain = forcingGain(t, cfg)
    gain = -cfg.omega * exp(cfg.velocityDecay * t) .* sin(cfg.omega * t);
end

% -------
function E = exactEnergyAt(t, cfg)
    E = 0.25 * amplitudeExact(t, cfg).^2;
end

% -------
function E = kineticEnergy(fem, state)
    E = 0.5 * (state.u' * fem.M * state.u + state.v' * fem.M * state.v)/fem.area;
end


% ============================================================
% NEWTON SOLVER
% ============================================================

function [next, converged, iterations, residual] = newtonStep(fem, cfg, coeff, dt, previous, time)

    n = fem.nVelocity;
    m = fem.nPressure;

    free  = fem.freeDofs;
    beta  = coeff.beta;
    alpha = coeff.alpha;
    gamma = coeff.gamma;
    gradDiv = cfg.gradDiv;

    historyU = -(beta(1)/(gamma * dt)) * previous.u + (beta(2) - beta(1) * (1 - gamma)/gamma) * previous.du + beta(3) * dt * previous.ddu + beta(4) * dt^2 * previous.dddu;
    historyV = -(beta(1)/(gamma * dt)) * previous.v + (beta(2) - beta(1) * (1 - gamma)/gamma) * previous.dv + beta(3) * dt * previous.ddv + beta(4) * dt^2 * previous.dddv;

    linearU = cfg.nu * fem.K * previous.u + gradDiv * (fem.Gxx  * previous.u + fem.Gxy * previous.v) - fem.Bx' * previous.p;
    linearV = cfg.nu * fem.K * previous.v + gradDiv * (fem.Gxy' * previous.u + fem.Gyy * previous.v) - fem.By' * previous.p;

    timeAlpha = (time - dt) + alpha * dt;
    gain = forcingGain(timeAlpha, cfg);

    rhsU = -fem.M * historyU - (1-alpha) * linearU + gain * fem.forceU;
    rhsV = -fem.M * historyV - (1-alpha) * linearV + gain * fem.forceV;
    rhs = [rhsU; rhsV; zeros(m,1)];

    base = beta(1)/(gamma * dt) * fem.M + alpha * cfg.nu * fem.K;
    Axx = base + alpha * gradDiv * fem.Gxx;
    Axy = alpha * gradDiv * fem.Gxy;
    Ayx = Axy';
    Ayy = base + alpha * gradDiv * fem.Gyy;
    tolerance = cfg.newtonTolerance * max(1, norm(rhs(free),inf));

    boundaryScale = amplitudeExact(time, cfg);
    boundaryU = boundaryScale * fem.boundaryU0;
    boundaryV = boundaryScale * fem.boundaryV0;

    next         = previous;
    converged    = false;
    iterations   = 0;
    residual     = inf;
    totalUpdates = 0;

    for newtonAttempt = 1:2

        if (newtonAttempt == 1)
            u = previous.u + dt * previous.du;
            v = previous.v + dt * previous.dv;
        else
            u = previous.u;
            v = previous.v;
        end

        p = previous.p;

        u(fem.fixedVelocity)  = boundaryU;
        v(fem.fixedVelocity)  = boundaryV;
        p(fem.pinnedPressure) = cfg.pressureValue;

        for iter = 1:cfg.maxNewtonIterations

            uAlpha = alpha * u + (1 - alpha) * previous.u;
            vAlpha = alpha * v + (1 - alpha) * previous.v;
            [C, Duu, Duv, Dvu, Dvv] = convection(fem,uAlpha,vAlpha);

            F = gaMomentumResidual(fem, Axx, Axy, Ayx, Ayy, C, rhs, u, v, p, previous, alpha);
            f = F(free);

            residual = norm(f,inf);
            if (isfinite(residual) == false); break; end

            if (residual < tolerance)
                next.u = u;
                next.v = v;
                next.p = p;

                next.du = (u - previous.u)/(gamma * dt) - ((1 - gamma)/gamma) * previous.du;
                next.dv = (v - previous.v)/(gamma * dt) - ((1 - gamma)/gamma) * previous.dv;

                next.ddu = (next.du - previous.du)/(gamma * dt) - ((1 - gamma)/gamma) * previous.ddu;
                next.ddv = (next.dv - previous.dv)/(gamma * dt) - ((1 - gamma)/gamma) * previous.ddv;

                next.dddu = (next.ddu - previous.ddu)/(gamma * dt) - ((1 - gamma)/gamma) * previous.dddu;
                next.dddv = (next.ddv - previous.ddv)/(gamma * dt) - ((1 - gamma)/gamma) * previous.dddv;

                converged  = true;
                iterations = totalUpdates;
                return;
            end

            J = [Axx + alpha * (C + Duu),        Axy + alpha * Duv,   -alpha * fem.Bx'; ...
                       Ayx + alpha * Dvu,  Ayy + alpha * (C + Dvv),   -alpha * fem.By'; ...
                          alpha * fem.Bx,           alpha * fem.By,   sparse(m,m)];

            delta = zeros(2 * n + m,1);
            delta(free) = -J(free,free)\f;

            if any(~isfinite(delta)); break; end

            correction = 1;
            accepted   = false;

            merit = norm(f,2);

            for trial = 1:20

                uTrial = u + correction * delta(1:n);
                vTrial = v + correction * delta(n + 1:2 * n);
                pTrial = p + correction * delta(2 * n + 1:end);

                uTrialAlpha = alpha * uTrial + (1 - alpha) * previous.u;
                vTrialAlpha = alpha * vTrial + (1 - alpha) * previous.v;

                CTrial = convection(fem, uTrialAlpha, vTrialAlpha);

                FTrial = gaMomentumResidual(fem, Axx, Axy, Ayx, Ayy, ...
                            CTrial, rhs, uTrial, vTrial, pTrial, previous, alpha);

                if (norm(FTrial(free),2) <= (1 - 1e-4 * correction) * merit)
                    accepted = true;
                    break;
                end

                correction = correction/2;
            end

            if (accepted == false); break; end

            u = uTrial;
            v = vTrial;
            p = pTrial;

            totalUpdates = totalUpdates + 1;
            iterations   = totalUpdates;
        end
    end
end

% -------
function F = gaMomentumResidual(fem, Axx, Axy, Ayx, Ayy, C, rhs, u, v, p, previous, alpha)

    n = fem.nVelocity;

    uAlpha = alpha * u + (1 - alpha) * previous.u;
    vAlpha = alpha * v + (1 - alpha) * previous.v;

    F = [Axx * u + Axy * v + C * uAlpha - alpha * fem.Bx' * p - rhs(1:n); ...
         Ayx * u + Ayy * v + C * vAlpha - alpha * fem.By' * p - rhs(n + 1:2 * n); ...
         fem.Bx * uAlpha + fem.By * vAlpha];
end

% -------
function [C, Duu, Duv, Dvu, Dvv] = convection(fem, u, v)

    [points, weights] = triangleQuadrature();

    nElements = numel(fem.elements);
    n = fem.nVelocity;

    withTangent = nargout > 1;

    rows = zeros(36 * nElements, 1);
    cols = zeros(36 * nElements, 1);
    cValues = zeros(36 * nElements, 1);

    if(withTangent)
        uuValues = cValues; uvValues = cValues;
        vuValues = cValues; vvValues = cValues;
    end

    for e = 1:nElements
        
        element = fem.elements(e);

        ids = element.nodes;
        ue  = u(ids);
        ve  = v(ids);
        Ce  = zeros(6);

        if(withTangent)
            DuuE = zeros(6);
            DuvE = zeros(6);
            DvuE = zeros(6);
            DvvE = zeros(6);
        end

        for q = 1:numel(weights)
            [N, Nx, Ny] = quadraticShape(points(q,:), element.gradL);

            weight = weights(q) * element.area;

            direction = (N' * ue) * Nx + (N' * ve) * Ny;

            localConvection = N * direction';

            Ce = Ce + 0.5 * weight * (localConvection - localConvection');

            if(withTangent)
                DuuE = DuuE + 0.5 * weight * (N * (Nx' * ue) - Nx * (N' * ue)) * N';
                DuvE = DuvE + 0.5 * weight * (N * (Ny' * ue) - Ny * (N' * ue)) * N';
                DvuE = DvuE + 0.5 * weight * (N * (Nx' * ve) - Nx * (N' * ve)) * N';
                DvvE = DvvE + 0.5 * weight * (N * (Ny' * ve) - Ny * (N' * ve)) * N';
            end
        end

        positions = (e - 1) * 36 + (1:36);

        [rr, cc] = ndgrid(ids, ids);
        rows(positions)    = rr(:);
        cols(positions)    = cc(:);
        cValues(positions) = Ce(:);

        if(withTangent)
            uuValues(positions) = DuuE(:);
            uvValues(positions) = DuvE(:);
            vuValues(positions) = DvuE(:);
            vvValues(positions) = DvvE(:);
        end
    end

    C = sparse(rows, cols, cValues, n, n);

    if(withTangent)
        Duu = sparse(rows, cols, uuValues, n, n);
        Duv = sparse(rows, cols, uvValues, n, n);
        Dvu = sparse(rows, cols, vuValues, n, n);
        Dvv = sparse(rows, cols, vvValues, n, n);
    end
end


% ============================================================
% PLOTTING
% ============================================================

function plotEnergyComparison(cfg, solution, method)

    time        = solution.time;
    exactEnergy = exactEnergyAt(time, cfg);
    error       = abs(solution.energy - exactEnergy);

    fig = figure(1);
    set(fig, 'Color', 'w', 'Position', [150 250 1200 400]);
    ax = getappdata(fig, 'EnergyAxes');

    if(isempty(ax) || any(~isgraphics(ax)))

        clf(fig);
        layout = tiledlayout(fig, 1, 2);

        ax1 = nexttile(layout);
        ax2 = nexttile(layout);

        hold(ax1, 'on');
        hold(ax2, 'on');

        exactTime = linspace(0, cfg.finalTime, 500)';
        plot(ax1, exactTime, exactEnergyAt(exactTime, cfg), 'k--', 'LineWidth', 2, 'DisplayName', 'Exact solution');

        setappdata(fig, 'EnergyAxes', [ax1, ax2]);

    else
        ax1 = ax(1);
        ax2 = ax(2);
    end

    plot(ax1, time, solution.energy, 'LineWidth', 2, 'DisplayName', method);
    xlabel(ax1, 'Time (s)');
    ylabel(ax1, 'Kinetic energy per unit area');

    semilogy(ax2, time, max(error, eps), 'LineWidth', 2, 'DisplayName', method);
    xlabel(ax2, 'Time (s)');
    ylabel(ax2, 'Absolute energy error');
    
    title(ax1, 'Kinetic energy');
    title(ax2, 'Energy error');

    legend(ax1, 'show', 'Location', 'northeast');
    legend(ax2, 'show', 'Location', 'northeast');

    grid(ax1, 'on');
    grid(ax2, 'on');

    ax1.Toolbar.Visible = 'off';
    ax2.Toolbar.Visible = 'off';

    exportgraphics(fig, fullfile(cfg.outputFolder, 'energy_comparison.png'), 'Resolution', 300);
end

% -------
function plotVelocityComparison(fem, cfg, solution, method, firstPlot)

    [~, index] = min(abs(fem.nodes(:,2) - cfg.L/6));
    yLine = fem.nodes(index,2);
    
    ids = find(abs(fem.velocityNodes(:,2) - yLine) < 1e-9*cfg.L);
    [xNodes, order] = sort(fem.velocityNodes(ids,1));
    ids = ids(order);

    time = solution.time(end);
   decay = amplitudeExact(time, cfg);

    x = linspace(0, cfg.L, 800)';
    exactU =  sin(x) .* cos(yLine) * decay;
    exactV = -cos(x) .* sin(yLine) * decay;

    fig = figure(2);
    if(firstPlot); clf(fig); end
    set(fig, 'Color', 'w', 'Position', [150 250 1200 400]);

    ax1 = subplot(1,2,1);
    ax2 = subplot(1,2,2);

    hold(ax1, 'on');
    hold(ax2, 'on');

    if(firstPlot)
        plot(ax1, x, exactU, 'k--', 'LineWidth', 2, 'DisplayName', 'Exact solution');
        plot(ax2, x, exactV, 'k--', 'LineWidth', 2, 'DisplayName', 'Exact solution');
    end

    plot(ax1, xNodes, solution.u(ids), 'LineWidth', 2, 'DisplayName', method);
    plot(ax2, xNodes, solution.v(ids), 'LineWidth', 2, 'DisplayName', method);

    xlabel(ax1, 'x');
    xlabel(ax2, 'x');
    ylabel(ax1, 'Horizontal velocity u');
    ylabel(ax2, 'Vertical velocity v');

    title(ax1, 'Horizontal velocity');
    title(ax2, 'Vertical velocity');

    legend(ax1, 'show', 'Location', 'northeast');
    legend(ax2, 'show', 'Location', 'northeast');

    grid(ax1, 'on');
    grid(ax2, 'on');

    ax1.Toolbar.Visible = 'off';
    ax2.Toolbar.Visible = 'off';

   exportgraphics(fig, fullfile(cfg.outputFolder, 'velocity_comparison.png'), 'Resolution', 300);
end


% ============================================================
% VIDEO GENERATION
% ============================================================
function createFlowAnimation(fem, cfg, solution)

    t = solution.snapTimes(:);
    n = cfg.tracerGridSize;
    x = linspace(0, cfg.L, n);

    [X, Y] = meshgrid(x, x);
    [P, inside] = velocityDisplayInterpolation(fem, X, Y);
    
    nSnap = numel(t);
    U = reshape(P*solution.uSnap, n, n, nSnap);
    V = reshape(P*solution.vSnap, n, n, nSnap);

    nTrail = max(2, round(cfg.tracerTrailFrames));
    s = cfg.L * ((1:cfg.particlesPerSide) - 0.5)/cfg.particlesPerSide;

    [xp, yp] = meshgrid(s, s);
    xp = xp(:);
    yp = yp(:);

    xTrail = NaN(numel(xp), nTrail);
    yTrail = NaN(numel(yp), nTrail);
    
    xTrail(:,1) = xp;
    yTrail(:,1) = yp;

    sizes = repmat(linspace(3,  0.5, nTrail), numel(xp), 1);
    fade  = repmat(linspace(1, 0.05, nTrail), numel(xp), 1);

    redBlue = interp1([1 128 256], [0.1 0.3 0.8; 1 1 1; 0.8 0.1 0.1], 1:256);

    fig = figure('Color', 'w', ...
              'Position', [50 120 1344 432], ...
                'Resize', 'off', ...
               'MenuBar', 'none', ...
               'ToolBar', 'none', ...
               'DefaultAxesFontSize', 11);

    tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'loose', 'Padding', 'compact');

    [~,   speedImg] = fieldPanel(tl, x, zeros(n), [0 1.05], turbo, 'Speed |u|', inside);
    [ax2,  vortImg] = fieldPanel(tl, x, zeros(n), [-2.05 2.05], redBlue, 'Vorticity \omega and particle trajectories', inside);
    
    hold(ax2, 'on');
    trails = scatter(ax2, xTrail(:), yTrail(:), sizes(:), 'k', 'filled', 'MarkerFaceAlpha', 'flat', 'AlphaData', fade(:));

    ax3   = nexttile(tl);
    exact = 0.25 * exp(-4 * cfg.nu * solution.time) .* cos(cfg.omega * solution.time).^2;
    
    hold(ax3, 'on');
    exactLine = plot(ax3, solution.time, exact, 'k-', 'LineWidth', 2);
    femLine   = plot(ax3, solution.time, solution.energy, 'b--', 'LineWidth', 2);
    energyDot = plot(ax3, NaN, NaN, 'ro', 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
    
    grid(ax3, 'on');
    pbaspect(ax3, [1 1 1]);

    xlim(ax3, t([1 end]).');
    ylim(ax3, [0 1.08 * max(exact)]);

    xlabel(ax3, 'Time (s)');
    ylabel(ax3, 'Kinetic energy per unit area');
    title(ax3, 'Kinetic energy');

    legend(ax3, [femLine exactLine], {'FEM', 'Exact solution'}, 'Location', 'northeast');
    heading = sgtitle(tl, '', 'FontWeight', 'bold');

    writer = VideoWriter(fullfile(cfg.outputFolder, 'TaylorGreen_Combined.mp4'), 'MPEG-4');
    writer.FrameRate = cfg.videoFrameRate;
    writer.Quality   = 100;
    open(writer);
    printOpt = sprintf('-r%d', round(5/3*get(groot, 'ScreenPixelsPerInch')));

    nFrames    = max(2, round(cfg.tracerVideoSeconds*cfg.videoFrameRate));
    frameTimes = linspace(t(1), t(end), nFrames);

    for k = 1:nFrames

        tk = frameTimes(k);
        j = min(find(t <= tk, 1, 'last'), nSnap-1);
        w = (tk - t(j))/(t(j+1) - t(j));

        Unow = (1 - w) * U(:,:,j) + w * U(:,:,j+1);
        Vnow = (1 - w) * V(:,:,j) + w * V(:,:,j+1);

        [~, dUdy] = gradient(Unow, x, x);
        [dVdx, ~] = gradient(Vnow, x, x);

        if(k > 1)
            h  = tk - frameTimes(k - 1);

            up = interp2(x, x, Unow, xp, yp);
            vp = interp2(x, x, Vnow, xp, yp);

            xp = mod(xp + h * up, cfg.L);
            yp = mod(yp + h * vp, cfg.L);

            xTrail = [xp, xTrail(:,1:end-1)];
            yTrail = [yp, yTrail(:,1:end-1)];
        end

        set(speedImg, 'CData', hypot(Unow, Vnow));
        set(vortImg,  'CData', dVdx - dUdy);
        set(trails,   'XData', xTrail(:), 'YData', yTrail(:));
        set(energyDot,'XData', tk, 'YData', interp1(solution.time, solution.energy, tk));
        heading.String = sprintf('Taylor-Green vortex | t = %.2f s', tk);
        drawnow;

        frame = print(fig, '-RGBImage', printOpt);
        writeVideo(writer, frame(1:end-mod(end,16), 1:end-mod(end,16), :));
    end

    close(writer);
    close(fig);
end

% -------
function [P, inside] = velocityDisplayInterpolation(fem, X, Y)

    mesh    = triangulation(fem.triangles, fem.nodes);
    points  = [X(:), Y(:)];
    element = pointLocation(mesh, points);
    valid   = find(~isnan(element));
    inside  = reshape(~isnan(element), size(X));

    B = cartesianToBarycentric(mesh, element(valid), points(valid,:));
    a = B(:,1);
    b = B(:,2);
    c = B(:,3);
    shape = [a.*(2*a-1), b.*(2*b-1), c.*(2*c-1), 4*a.*b, 4*b.*c, 4*c.*a];

    cols = fem.velocityTriangles(element(valid), :);
    P = sparse(repelem(valid, 6), reshape(cols.', [], 1), reshape(shape.', [], 1), numel(X), fem.nVelocity);
end

% -------
function [ax, img] = fieldPanel(layout, x, data, limits, map, name, inside)
    ax  = nexttile(layout);
    img = imagesc(ax, x, x, data, 'AlphaData', double(inside));

    set( ax, 'YDir', 'normal');
    axis(ax, 'image');
    clim(ax, limits);

    colormap(ax, map);
    colorbar(ax);

    xlabel(ax, 'x');
    ylabel(ax, 'y');

    title(ax, name);
end