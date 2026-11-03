
clc; clear; close all;

% -------------------------------------------------
% 2D lid-driven cavity with P2/P1 finite
% -------------------------------------------------

cfg.method = 'GA-234';
cfg.rhoInfinity = 0;

cfg.N = 33;
cfg.stretch = 1.6;
cfg.ReList = [100 400 1000];

cfg.dt = 0.2;
cfg.finalTime = 40;

cfg.gradDiv = 0.03;

cfg.newtonTolerance = 1e-8;
cfg.maxNewtonIterations = 22;

cfg.pressureNode = 1;
cfg.pressureValue = 0;

cfg.snapshotEvery = 2;

cfg.saveVideo      = true;
cfg.videoFrames    = 300;
cfg.videoFrameRate = 30;
cfg.videoWidth     = 1600;

cfg.displayGridSize = 400;

cfg.outputFolder = 'LidDrivenCavity';

if (isfolder(cfg.outputFolder) == false)
    mkdir(cfg.outputFolder);
end

% ------------------------------------------
% Finite Element Mesh and Initial Conditions
% ------------------------------------------

fprintf('\n-----------------------------------------------------');
fprintf('\nLid-driven cavity | P2/P1 finite element simulation\n');
fprintf('-----------------------------------------------------\n');

fem = buildMesh(cfg);
fem = setupBoundaryDOFs(fem, cfg);

xu = fem.velocityNodes(:,1);
yu = fem.velocityNodes(:,2);

lid = (abs(yu - 1) < 1e-12);

initialU = zeros(fem.nVelocity, 1);
initialU(lid) = 1;

zeroVelocity = zeros(fem.nVelocity, 1);

fem.initialState.u = initialU;
fem.initialState.v = zeroVelocity;
fem.initialState.p = zeros(fem.nPressure, 1);

fem.initialState.du = zeroVelocity;
fem.initialState.dv = zeroVelocity;

fem.initialState.ddu = zeroVelocity;
fem.initialState.ddv = zeroVelocity;

fem.initialState.dddu = zeroVelocity;
fem.initialState.dddv = zeroVelocity;

fem.boundaryU0 = initialU(fem.fixedVelocity);
fem.boundaryV0 = zeros(numel(fem.fixedVelocity), 1);

fem = assembleFEMMatrices(fem);

% ------------------------------------------
% Reynolds Number Comparison
% ------------------------------------------

nRe = numel(cfg.ReList);
solutions  = cell(nRe, 1);

for k = 1:nRe

    runCfg = cfg;
    runCfg.Re = cfg.ReList(k);
    runCfg.nu = 1/runCfg.Re;

    fprintf('\nRe = %d | Method: %s (rho Infinity = %.2f) | dt = %.4f\n', runCfg.Re, cfg.method, cfg.rhoInfinity, cfg.dt);

    solutions{k}  = solveFlow(fem, runCfg, cfg.method, cfg.dt);
end

% ------------------------------------------
% Plotting
% ------------------------------------------

plotCenterlineComparison(fem, cfg, solutions);

if (cfg.saveVideo)
    createCavityAnimation(fem, cfg, solutions);
end

fprintf('\nSimulation completed successfully. Results saved to: %s\n', cfg.outputFolder);

% ============================================================
% FINITE ELEMENT SETUP
% ============================================================

function fem = buildMesh(cfg)

    n  = cfg.N;
    xi = linspace(0, 1, n);

    if (cfg.stretch > 0)
        xi = 0.5 * (1 + tanh(cfg.stretch * (xi - 0.5))/tanh(cfg.stretch/2));
    end

    [X, Y] = ndgrid(xi, xi);
    nodes = [X(:), Y(:)];

    triangles = zeros(2 * (n-1)^2, 3);
    row = 0;

    for j = 1:n-1
        for i = 1:n-1

            a = (j-1) * n + i;
            b = a + 1;
            d = a + n;
            e = d + 1;

            row = row + 1;
            triangles(row,:) = [a, b, e];

            row = row + 1;
            triangles(row,:) = [a, e, d];
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

    nv = fem.nVelocity;
    np = fem.nPressure;

    rows = zeros(36 * nElements, 1);
    cols = zeros(36 * nElements, 1);

    mValues   = zeros(36 * nElements, 1);
    kValues   = zeros(36 * nElements, 1);

    gxxValues = zeros(36 * nElements, 1);
    gxyValues = zeros(36 * nElements, 1);
    gyyValues = zeros(36 * nElements, 1);

    pRows = zeros(18 * nElements, 1);
    pCols = zeros(18 * nElements, 1);

    bxValues = zeros(18 * nElements, 1);
    byValues = zeros(18 * nElements, 1);

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

        position = (e - 1) * 36 + (1:36);

        [rr, cc] = ndgrid(uNodes, uNodes);
        rows(position) = rr(:);
        cols(position) = cc(:);

        mValues(position) = Me(:);
        kValues(position) = Ke(:);

        gxxValues(position) = Gxxe(:);
        gxyValues(position) = Gxye(:);
        gyyValues(position) = Gyye(:);

        position = (e - 1) * 18 + (1:18);

        [pr, pc] = ndgrid(pNodes, uNodes);
        pRows(position) = pr(:);
        pCols(position) = pc(:);

        bxValues(position) = Bxe(:);
        byValues(position) = Bye(:);
    end

    fem.M = sparse(rows, cols, mValues, nv, nv);
    fem.K = sparse(rows, cols, kValues, nv, nv);

    fem.Bx = sparse(pRows, pCols, bxValues, np, nv);
    fem.By = sparse(pRows, pCols, byValues, np, nv);

    fem.Gxx = sparse(rows, cols, gxxValues, nv, nv);
    fem.Gxy = sparse(rows, cols, gxyValues, nv, nv);
    fem.Gyy = sparse(rows, cols, gyyValues, nv, nv);

    fprintf('\nFinite element assembly completed\n');
    fprintf('  Total degrees of freedom = %d\n', 2*nv + np);
    fprintf('  Total mesh elements      = %d\n', nElements);
    fprintf('  P1 mesh vertices         = %d\n', np);
    fprintf('  Velocity nodes (P2)      = %d\n', nv);
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

function result = solveFlow(fem, cfg, method, dt)

    steps = round(cfg.finalTime/dt);
    coeff = timeCoefficients(method, cfg.rhoInfinity);
    state = fem.initialState;

    result.newtonIterations = 0;
    currentTime = 0;

    maxSnapshots = floor(steps/cfg.snapshotEvery) + 2;

    result.snapTimes = zeros(maxSnapshots, 1);
    result.uSnap = zeros(fem.nVelocity, maxSnapshots);
    result.vSnap = zeros(fem.nVelocity, maxSnapshots);

    nSnapshots = 1;
    result.uSnap(:,1) = state.u;
    result.vSnap(:,1) = state.v;

    for step = 1:steps

        previous    = state;
        currentTime = step * dt;

        [candidate, ok, iters, residual] = newtonStep(fem, cfg, coeff, dt, previous);
        result.newtonIterations = result.newtonIterations + iters;

        if (ok == false)
            error('%s: Newton solver failed at t = %.4f (residual = %.3e). Try reducing cfg.dt.', method, currentTime, residual);
        end

        state = candidate;

        if (mod(step, cfg.snapshotEvery) == 0)
            nSnapshots = nSnapshots + 1;
            result.snapTimes(nSnapshots) = currentTime;
            result.uSnap(:,nSnapshots) = state.u;
            result.vSnap(:,nSnapshots) = state.v;
        end

        if (mod(step, max(1, round(steps/8))) == 0 || step == 1 || step == steps)
            fprintf('  %-6s | Step %3d/%3d | Time: %6.2f s | Newton iterations: %2d | Residual: %9.2e\n', method, step, steps, currentTime, iters, residual);
        end
    end

    result.time = currentTime;

    if (abs(result.snapTimes(nSnapshots) - currentTime) > 1e-8)
        nSnapshots = nSnapshots + 1;
        result.snapTimes(nSnapshots) = currentTime;
        result.uSnap(:,nSnapshots) = state.u;
        result.vSnap(:,nSnapshots) = state.v;
    end

    result.snapTimes = result.snapTimes(1:nSnapshots);
    result.uSnap = result.uSnap(:,1:nSnapshots);
    result.vSnap = result.vSnap(:,1:nSnapshots);

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


% ============================================================
% NEWTON SOLVER
% ============================================================

function [next, converged, iterations, residual] = newtonStep(fem, cfg, coeff, dt, previous)

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

    rhsU = -fem.M * historyU - (1-alpha) * linearU;
    rhsV = -fem.M * historyV - (1-alpha) * linearV;
    rhs = [rhsU; rhsV; zeros(m,1)];

    base = beta(1)/(gamma * dt) * fem.M + alpha * cfg.nu * fem.K;
    Axx = base + alpha * gradDiv * fem.Gxx;
    Axy = alpha * gradDiv * fem.Gxy;
    Ayx = Axy';
    Ayy = base + alpha * gradDiv * fem.Gyy;
    tolerance = cfg.newtonTolerance * max(1, norm(rhs(free),inf));

    boundaryU = fem.boundaryU0;
    boundaryV = fem.boundaryV0;

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

                trialMerit = norm(FTrial(free),2);

                if (isfinite(trialMerit) && trialMerit <= (1 - 1e-4 * correction) * merit)
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
% GHIA ET AL. (1982) REFERENCE DATA
% ============================================================

function reference = ghiaReference(re)

    y = [1; .9766; .9688; .9609; .9531; .8516; .7344; .6172; .5; .4531; .2813; .1719; .1016; .0703; .0625; .0547; 0];
    x = [1; .9688; .9609; .9531; .9453; .9063; .8594; .8047; .5; .2344; .2266; .1563; .0938; .0781; .0703; .0625; 0];

    switch re
        case 100
            u = [1;  .84123;  .78871;  .73722;  .68717;  .23151;  .00332; -.13641; -.20581; -.21090; -.15662; -.10150; -.06434; -.04775; -.04192; -.03717; 0];
            v = [0; -.05906; -.07391; -.08864; -.10313; -.16914; -.22445; -.24533;  .05454;  .17527;  .17507;  .16077;  .12317;  .10890;  .10091;  .09233; 0];

        case 400
            u = [1;  .75837;  .68439;  .61756;  .55892;  .29093;  .16256;  .02135; -.11477; -.17119; -.32726; -.24299; -.14612; -.10338; -.09266; -.08186; 0];
            v = [0; -.12146; -.15663; -.19254; -.22847; -.23827; -.44993; -.38598;  .05186;  .30174;  .30203;  .28124;  .22965;  .20920;  .19713;  .18360; 0];

        case 1000
            u = [1;  .65928;  .57492;  .51117;  .46604;  .33304;  .18719;  .05702; -.06080; -.10648; -.27805; -.38289; -.29730; -.22220; -.20196; -.18109; 0];
            v = [0; -.21388; -.27669; -.33714; -.39188; -.51500; -.42665; -.31966;  .02526;  .32235;  .33075;  .37095;  .32627;  .30353;  .29012;  .27485; 0];

        otherwise
            error('Ghia reference implemented for Re = 100, 400, 1000 only.');
    end


    [reference.y, order] = sort(y);
    reference.u = u(order);

    [reference.x, order] = sort(x);
    reference.v = v(order);
end

% ============================================================
% PLOTTING
% ============================================================

function plotCenterlineComparison(fem, cfg, solutions)

    nRe = numel(cfg.ReList);

    xu = fem.velocityNodes(:,1);
    yu = fem.velocityNodes(:,2);

    verticalNodes   = find(abs(xu - 0.5) < 1e-10);
    horizontalNodes = find(abs(yu - 0.5) < 1e-10);

    [yLine, orderU] = sort(yu(verticalNodes));
    [xLine, orderV] = sort(xu(horizontalNodes));

    verticalNodes   = verticalNodes(orderU);
    horizontalNodes = horizontalNodes(orderV);

    uColor = [0.70 0.09 0.17];
    vColor = [0.13 0.40 0.67];

    uLimits = [-0.5, 1.0];
    vLimits = [-0.6, 0.4];


    figureWidth  = 1340;
    figureHeight = 480;

    leftMargin  = 85;
    rightMargin = 67;
    gap = 124;

    panelBottom = 123;
    panelHeight = 258;

    panelWidth = (figureWidth - leftMargin - rightMargin - (nRe - 1)*gap)/nRe;

    fig = figure(1);
    clf(fig);

    set(fig, 'Color'               , 'w', ...
             'Units'               ,   'pixels', ...
             'Position'            , [55 85 figureWidth figureHeight], ...
             'PaperPositionMode'   , 'auto', ...
             'DefaultAxesFontName' , 'Helvetica', ...
             'DefaultAxesFontSize' , 11, ...
             'DefaultTextFontName' , 'Helvetica');

    annotation(fig, 'textbox'            , [0.08 0.902 0.84 0.052], ...
                    'String'             , sprintf('Lid-driven cavity | %s vs Ghia et al. (1982)', cfg.method), ...
                    'EdgeColor'          , 'none', ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment'  , 'middle', ...
                    'FontWeight'         , 'bold', ...
                    'FontSize'           , 15, ...
                    'Interpreter', 'none');

    firstU = [];

    positions = zeros(nRe,4);
    axesV = gobjects(nRe,1);
    axesU = gobjects(nRe,1);

    for k = 1:nRe

        left = leftMargin + (k - 1)*(panelWidth + gap);

        position = [left panelBottom panelWidth panelHeight];

        positions(k,:) = position;

        axV = axes(fig, 'Units', 'pixels', 'Position', position, 'PositionConstraint', 'innerposition');
        
        hold(axV, 'on');
        
        axU = axes(fig, 'Units', 'pixels', 'Position', position, 'PositionConstraint', 'innerposition', 'Color', 'none');

        hold(axU, 'on');

        axV.Toolbar.Visible = 'off';
        axU.Toolbar.Visible = 'off';

        vLine = solutions{k}.v(horizontalNodes);
        uLine = solutions{k}.u(verticalNodes);

        reference = ghiaReference(cfg.ReList(k));

        plot(axV, xLine, vLine, '-', 'Color', vColor, 'LineWidth', 2.3);

        plot(axV, reference.x, reference.v, 'o', 'MarkerSize', 6.2, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', vColor, 'LineWidth', 1.45);
        
        plot(axU, uLine, yLine, '-', 'Color', uColor, 'LineWidth', 2.3);

        plot(axU, reference.u, reference.y, 'o', 'MarkerSize', 6.2, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', uColor, 'LineWidth', 1.45);

        set(axV, ...
            'XLim'     , [0 1], ...
            'YLim'     , vLimits, ...
            'XTick'    , 0:0.2:1, ...
            'YTick'    , -0.6:0.2:0.4, ...
            'XColor'   , vColor, ...
            'YColor'   , vColor, ...
            'Box'      , 'off', ...
            'TickDir'  , 'out', ...
            'LineWidth', 0.95, ...
            'XGrid'    , 'on', ...
            'YGrid'    , 'on', ...
            'GridColor', [0.82 0.84 0.87], ...
            'GridAlpha', 0.65, ...
            'Layer'    , 'bottom');

        set(axU, ...
            'XLim'         , uLimits, ...
            'YLim'         , [0 1], ...
            'XTick'        , -0.5:0.5:1, ...
            'YTick'        , 0:0.2:1, ...
            'XColor'       , uColor, ...
            'YColor'       , uColor, ...
            'XAxisLocation', 'top', ...
            'YAxisLocation', 'right', ...
            'Box'          , 'off', ...
            'TickDir'      , 'out', ...
            'LineWidth'    , 0.95);

        xlabel(axV, 'x', 'FontSize', 12);
        ylabel(axV, 'v(x, 0.5)', 'FontSize', 12, 'Color', vColor);

        xlabel(axU, 'u(0.5, y)', 'FontSize', 12, 'Color', uColor);
        ylabel(axU, 'y', 'FontSize', 12, 'Color', uColor);

        set(axV, 'Units', 'pixels', 'Position', position);
        set(axU, 'Units', 'pixels', 'Position', position);


        annotation(fig, 'textbox', [left/figureWidth, 0.045, panelWidth/figureWidth, 0.041], ...
                         'String', sprintf('Re = %d', cfg.ReList(k)), ...
                      'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', ...
              'VerticalAlignment', 'middle', ...
                     'FontWeight', 'bold', ...
                       'FontSize', 13, ...
                          'Color', [0.15 0.15 0.15], ...
                    'Interpreter', 'none');

        if k == 1
            firstU = axU;
        end

        axesV(k) = axV;
        axesU(k) = axU;

    end

    hU = plot(firstU, NaN, NaN, '-', 'Color', uColor, 'LineWidth', 2.3);
    hV = plot(firstU, NaN, NaN, '-', 'Color', vColor, 'LineWidth', 2.3);
    hG = plot(firstU, NaN, NaN, 'o', 'MarkerSize', 6.2, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', [0.20 0.20 0.20], 'LineWidth', 1.4);

    lgd = legend(firstU,  [hU hV hG], {sprintf('u(0.5, y)  %s', cfg.method), sprintf('v(x, 0.5)  %s', cfg.method), 'Ghia et al. (1982)'}, ...
             'Location', 'southeast', 'Orientation', 'vertical', 'Box', 'on', ...
             'FontSize', 9, 'Interpreter', 'none', 'AutoUpdate', 'off');

    lgd.Color = 'white';
    lgd.EdgeColor = [0.72 0.72 0.72];
    lgd.LineWidth = 0.65;

    drawnow;

    for k = 1:nRe
        set(axesV(k), 'Units', 'pixels', 'Position', positions(k,:));
        set(axesU(k), 'Units', 'pixels', 'Position', positions(k,:));
    end

    lgd.Units = 'pixels';
    legendPosition = lgd.Position;

    legendPosition(1) = positions(1,1) + positions(1,3) - legendPosition(3) - 10;
    legendPosition(2) = positions(1,2) + 10;
    
    lgd.Position = legendPosition;


    set(fig, 'PaperPositionMode', 'auto');

    print(fig, fullfile(cfg.outputFolder, 'centerline_comparison.png'), '-dpng', '-r300');
end


% ============================================================
% VIDEO GENERATION
% ============================================================

function createCavityAnimation(fem, cfg, solutions)

    nRe = numel(cfg.ReList);
    n = max(400, cfg.displayGridSize);
    x = linspace(0, 1, n);
    y = x;

    [X, Y] = meshgrid(x, y);
      mesh = triangulation(fem.triangles, fem.nodes);

    samplePoints = min(max([X(:), Y(:)], 1e-12), 1-1e-12);
    [element, bary] = pointLocation(mesh, samplePoints);

    a = bary(:,1);
    b = bary(:,2);
    c = bary(:,3);

    shape = [a.*(2*a-1), b.*(2*b-1), c.*(2*c-1), 4*a.*b, 4*b.*c, 4*c.*a];

    nodeIds = fem.velocityTriangles(element,:);

    P = sparse(repelem((1:numel(X))',6), reshape(nodeIds.',[],1), reshape(shape.',[],1), numel(X), fem.nVelocity);

    rowNames = {'Speed |u|', 'Horizontal velocity u', 'Vertical velocity v'};

    limits = [inf, -inf; inf, -inf; inf, -inf];

    for k = 1:nRe
        values = {hypot(solutions{k}.uSnap, solutions{k}.vSnap), solutions{k}.uSnap, solutions{k}.vSnap};

        for r = 1:3
            limits(r,1) = min(limits(r,1), min(values{r}(:)));
            limits(r,2) = max(limits(r,2), max(values{r}(:)));
        end
    end

    maps = {sequentialMap(256), divergingMap(limits(2,:),256), divergingMap(limits(3,:),256)};

    streamLevels = [-0.1175 -0.115 -0.11 -0.1 -0.09 -0.07 -0.05 -0.03 -0.01 -0.001 -0.0001 0.0001 0.0005 0.001 0.003];

    panel = 300;
    gap = 22;

    leftMargin = 90;
    rightMargin = 140;
    topMargin = 120;
    bottomMargin = 75;

    figureWidth  = leftMargin + nRe*panel + (nRe-1)*gap + rightMargin;
    figureHeight = bottomMargin + 3*panel + 2*gap + topMargin;

    fig = figure('Color'              , 'w', ...
                 'Visible'            , 'off', ...
                 'Units'              , 'pixels', ...
                 'Position'           , [40 40 figureWidth figureHeight], ...
                 'Resize'             , 'off', ...
                 'MenuBar'            , 'none', ...
                 'ToolBar'            , 'none', ...
                 'PaperPositionMode'  , 'auto', ...
                 'DefaultAxesFontName', 'Helvetica', ...
                 'DefaultAxesFontSize', 11, ...
                 'DefaultTextFontName', 'Helvetica');

    panels = gobjects(3,nRe);
    images = gobjects(3,nRe);

    for r = 1:3
        for c = 1:nRe
            position = [(leftMargin + (c-1) * (panel + gap)), (bottomMargin + (3 - r) * (panel + gap)), panel, panel];

             ax = axes(fig, 'Units', 'pixels', 'Position', position);
            img = imagesc(ax, x, y, zeros(n));
            hold(ax, 'on');

            rectangle(ax, 'Position', [0 0 1 1], 'EdgeColor', [0.25 0.25 0.28], 'LineWidth', 0.8);

            set(ax, ...
            'YDir', 'normal', 'DataAspectRatio',   [1 1 1], 'XLim',   [0 1], 'YLim', [0 1], ...
           'XTick',  0:0.5:1, 'YTick',  0:0.5:1, 'TickDir',  'out', 'Layer', 'top', ...
             'Box', 'off', 'LineWidth', 0.8);

            clim(ax, limits(r,:));
            colormap(ax, maps{r});

            if(r < 3)
                set(ax, 'XTickLabel', []);
            else
                xlabel(ax, 'x');
            end

            if(c > 1)
                set(ax, 'YTickLabel', []);
            else
                ylabel(ax, 'y');
            end

            if(r == 1)
                title(ax, sprintf('Re = %d', cfg.ReList(c)), 'FontWeight', 'bold', 'FontSize', 14, 'Color', [0.13 0.13 0.13], 'Interpreter', 'none');
            end

            panels(r,c) = ax;
            images(r,c) = img;

        end

        cb = colorbar(panels(r,nRe));

        set(cb, 'Units', 'pixels', 'Position', [(leftMargin + nRe * panel + (nRe - 1) * gap + 16), (bottomMargin + (3 - r) * (panel + gap)), 14, panel], ...
                'TickDirection', 'out', 'LineWidth', 0.8);
        ylabel(cb, rowNames{r});
    end


    for r = 1:3
        for c = 1:nRe
            position = [(leftMargin + (c - 1) * (panel + gap)), (bottomMargin + (3 - r) * (panel + gap)), panel, panel];
            set(panels(r,c), 'Units', 'pixels', 'Position', position);
        end
    end

    mainTitle = annotation(fig, 'textbox', [0.08 0.942 0.84 0.045], ...
                                'String' , sprintf('Lid-driven cavity | %s   |   t = 0.0', cfg.method), ...
                              'EdgeColor', 'none', ...
                    'HorizontalAlignment', 'center', ...
                      'VerticalAlignment', 'middle', ...
                               'FontName', 'Helvetica', ...
                             'FontWeight', 'bold', ...
                               'FontSize', 16, ...
                            'Interpreter', 'none');

    streamLines = gobjects(1,nRe);
      zeroLines = gobjects(2,nRe);

    printOpt = sprintf('-r%d', max(80, round(cfg.videoWidth/figureWidth * get(groot, 'ScreenPixelsPerInch'))));
    writer = VideoWriter(fullfile(cfg.outputFolder,'cavity_flow.mp4'), 'MPEG-4');
    writer.FrameRate = cfg.videoFrameRate;
    writer.Quality = 100;
    open(writer);

    for frame = 1:cfg.videoFrames
            t = cfg.finalTime * (frame - 1) /max(1,cfg.videoFrames - 1);
        for c = 1:nRe

            solution = solutions{c};
            u = reshape(P * interpolateSnapshot(solution.snapTimes, solution.uSnap, t), n, n);
            v = reshape(P * interpolateSnapshot(solution.snapTimes, solution.vSnap, t), n, n);

            set(images(1,c), 'CData', hypot(u,v));
            set(images(2,c), 'CData', u);
            set(images(3,c), 'CData', v);

            for r = 2:3

                if isgraphics(zeroLines(r-1,c)); delete(zeroLines(r-1,c)); end

                if (r == 2)
                    field = u;
                else
                    field = v;
                end

                if(min(field(:)) < 0 && max(field(:)) > 0)
                    [~,zeroLines(r-1,c)] = contour(panels(r,c), X, Y, field, [0 0], 'LineColor', [0.12 0.12 0.14],'LineWidth', 0.8);
                end

            end

            psi = cumtrapz(y, u, 1);

            if(isgraphics(streamLines(c))); delete(streamLines(c)); end

            if max(abs(psi(:))) > min(abs(streamLevels))
                [~,streamLines(c)] = contour(panels(1,c), X, Y, psi, streamLevels, 'LineColor', [0.97 0.97 0.97], 'LineWidth', 0.6);
            end

        end

        mainTitle.String = sprintf('Lid-driven cavity | %s   |   t = %.1f', cfg.method, t);
        drawnow;

        videoFrame = print(fig, '-RGBImage', printOpt);

        if(~isa(videoFrame, 'uint8'))

            if isfloat(videoFrame) && max(videoFrame(:)) <= 1
                videoFrame = uint8(255 * videoFrame);
            else
                videoFrame = uint8(videoFrame);
            end

        end

        frameHeight = size(videoFrame,1);
        frameWidth  = size(videoFrame,2);

        videoFrame = videoFrame(1:frameHeight-mod(frameHeight,2), 1:frameWidth-mod(frameWidth,2), :);

        writeVideo(writer, videoFrame);
    end
    close(writer);
    close(fig);
end

% -------
function field = interpolateSnapshot(times, snapshots, t)

    if(t <= times(1))
        field = snapshots(:,1);
        return;
    end

    if(t >= times(end))
        field = snapshots(:,end);
        return;
    end

    i = find(times <= t, 1, 'last');

    theta = (t - times(i))/(times(i + 1) - times(i));
    field = (1 - theta) * snapshots(:,i) + theta * snapshots(:,i + 1);
end

% -------
function map = sequentialMap(n)

    anchors = [0.267 0.005 0.329; 0.283 0.141 0.458; 0.254 0.265 0.530; ...
               0.207 0.372 0.553; 0.164 0.471 0.558; 0.128 0.567 0.551; ...
               0.135 0.659 0.518; 0.267 0.749 0.441; 0.478 0.821 0.318; ...
               0.741 0.873 0.150; 0.993 0.906 0.144];

    map = interp1(linspace(0,1,size(anchors,1)), anchors, linspace(0,1,n)');
end

% -------
function map = divergingMap(limits, n)

    blue = [1.000 1.000 1.000; ...
            0.573 0.773 0.871; ...
            0.263 0.576 0.765; ...
            0.129 0.400 0.675; ...
            0.020 0.188 0.380];

    red = [1.000 1.000 1.000; ...
           0.957 0.647 0.510; ...
           0.839 0.376 0.302; ...
           0.698 0.094 0.168; ...
           0.404 0.000 0.122];

    anchors = linspace(0,1,size(blue,1));

    value = linspace(limits(1),limits(2),n)';

    negative = value < 0;
    positive = value > 0;

    map = ones(n,3);

    for c = 1:3
        map(negative,c) = interp1(anchors, blue(:,c), min(1,value(negative)/limits(1)));
        map(positive,c) = interp1(anchors, red(:,c),  min(1,value(positive)/limits(2)));
    end

end
