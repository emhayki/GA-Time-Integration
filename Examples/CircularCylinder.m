
clc; clear; close all;

% -------------------------------------------------
% 2D flow past a circular cylinder with P2/P1 finite elements
% -------------------------------------------------

cfg.Lx = 2.2;
cfg.Ly = 0.41;
cfg.meshSize = 0.055;

cfg.cx = 0.20;
cfg.cy = 0.20;
cfg.R  = 0.05;
cfg.D  = 2 * cfg.R;

cfg.ringPoints = 48;
cfg.ringRadii = [cfg.R, 0.059, 0.071, 0.087, 0.109, 0.134, 0.160];
cfg.wakeEnd = 1.55;
cfg.wakeHalfWidth = 0.13;

cfg.Re = 100;
cfg.Umean = 1.0;
cfg.Umax = 1.5 * cfg.Umean;
cfg.nu = cfg.Umean * cfg.D / cfg.Re;

cfg.finalTime = 12;
cfg.dt = 0.04;
cfg.rampTime = 1.2;
cfg.rhoInfinity = 0;

cfg.gradDiv = 0.03;

cfg.newtonTolerance = 1e-9;
cfg.maxNewtonIterations = 22;

cfg.pressureNode = [];
cfg.pressureValue = 0;
cfg.backflowSmoothing = 0.01;

cfg.outputFolder = 'CircularCylinder5';

cfg.saveVideo = true;
cfg.videoFrameRate = 30;
cfg.videoSpeed = 0.5;
cfg.videoWidth = 1920;
cfg.displayHeight = 420;

cfg.vorticitySigma = 0.5;
cfg.vorticityLimit = [];

if(isfolder(cfg.outputFolder) == false)
    mkdir(cfg.outputFolder);
end

% ------------------------------------------
% Finite Element Mesh and Initial Conditions
% ------------------------------------------

fprintf('\n-------------------------------------------------------');
fprintf('\nCircular cylinder | P2/P1 finite element simulation\n');
fprintf('-------------------------------------------------------\n');

fem = buildMesh(cfg);
fem = setupBoundaryDOFs(fem,cfg);
fem = setupCylinder(fem,cfg);

xu = fem.velocityNodes(:,1);
yu = fem.velocityNodes(:,2);

inlet = (abs(xu) < 1e-8);
eta = yu(abs(xu) < 1e-8)/cfg.Ly;

inletU = zeros(fem.nVelocity,1);
inletU(inlet) = cfg.Umax * 4 * eta .* (1-eta);

fem.boundaryU0 = inletU(fem.fixedVelocity);
fem.boundaryV0 = zeros(numel(fem.fixedVelocity),1);

zeroVelocity = zeros(fem.nVelocity,1);

fem.initialState.u = zeroVelocity;
fem.initialState.v = zeroVelocity;

fem.initialState.p = zeros(fem.nPressure,1);

fem.initialState.du = zeroVelocity;
fem.initialState.dv = zeroVelocity;

fem.initialState.ddu = zeroVelocity;
fem.initialState.ddv = zeroVelocity;

fem.initialState.dddu = zeroVelocity;
fem.initialState.dddv = zeroVelocity;

fem = assembleFEMMatrices(fem);

% ------------------------------------------
% Reference Solution and Method Comparison
% ------------------------------------------
referenceDt = cfg.dt/4;
snapshotTimes = (0:referenceDt:cfg.finalTime)';

fprintf('\nReference solution (GA-234, rho Infinity = %.2f) | dt = %.4f\n', cfg.rhoInfinity, referenceDt);
reference = solveFlow(fem, cfg, 'GA-234', referenceDt, snapshotTimes);

methods = {'GA-2', 'GA-23', 'GA-234'};
solutions = cell(numel(methods),1);

for k = 1:numel(methods)
    fprintf('\nMethod: %s (rho Infinity = %.2f)| dt = %.4f\n', methods{k}, cfg.rhoInfinity, cfg.dt);
  
    solutions{k} = solveFlow(fem,cfg,methods{k},cfg.dt,[]);
end

% ------------------------------------------
% Plotting
% ------------------------------------------

plotResponseComparison(cfg, reference, solutions, methods);

if(cfg.saveVideo); createFlowAnimation(fem, cfg, reference); end

fprintf('\nSimulation completed. Results saved to: %s\n', cfg.outputFolder);

% ============================================================
% FINITE ELEMENT SETUP
% ============================================================

function fem = buildMesh(cfg)

    [nodes, triangles] = cylinderMesh(cfg);

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

    a = nodes(triangles(:,1),:);
    b = nodes(triangles(:,2),:);
    c = nodes(triangles(:,3),:);

    fem.nElements = nElements;
    fem.area = sum(abs((b(:,1)-a(:,1)).*(c(:,2)-a(:,2)) ...
                     - (b(:,2)-a(:,2)).*(c(:,1)-a(:,1))))/2;
end

% -------
function fem = setupBoundaryDOFs(fem, cfg)

    boundaryEdges = freeBoundary(triangulation(fem.triangles, fem.nodes));
    [~, boundaryEdgeIndex] = ismember(sort(boundaryEdges,2), fem.edges, 'rows');

    midpoints = (fem.nodes(boundaryEdges(:,1),:) + fem.nodes(boundaryEdges(:,2),:))/2;
    onOutlet  = abs(midpoints(:,1) - cfg.Lx) < 1e-8;

    fixedEdges = boundaryEdges(~onOutlet,:);
    fixedEdgeIndex = boundaryEdgeIndex(~onOutlet);

    fem.fixedVelocity  = unique([fixedEdges(:); fem.nPressure + fixedEdgeIndex]);
    fem.pinnedPressure = cfg.pressureNode;

    fixed = [fem.fixedVelocity; ...
             fem.nVelocity + fem.fixedVelocity; ...
             2 * fem.nVelocity + fem.pinnedPressure];

    fem.freeDofs = setdiff((1:2 * fem.nVelocity + fem.nPressure)', fixed);

    outletVertices = boundaryEdges(onOutlet,:);
    fem.outletEdges = [outletVertices, fem.nPressure + boundaryEdgeIndex(onOutlet)];
    fem.outletLength = abs(fem.nodes(outletVertices(:,2),2) ...
                         - fem.nodes(outletVertices(:,1),2));
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
    fprintf('  Total mesh nodes         = %d\n', np);
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

% -------
function fem = setupCylinder(fem, cfg)

    boundaryEdges = freeBoundary(triangulation(fem.triangles, fem.nodes));
    a = fem.nodes(boundaryEdges(:,1),:);
    b = fem.nodes(boundaryEdges(:,2),:);

    ra = hypot(a(:,1) - cfg.cx, a(:,2) - cfg.cy);
    rb = hypot(b(:,1) - cfg.cx, b(:,2) - cfg.cy);
    onCylinder = abs(ra - cfg.R) < 1e-8 & abs(rb - cfg.R) < 1e-8;

    fem.bodyEdges = boundaryEdges(onCylinder,:);

    fem.bodyElements = zeros(size(fem.bodyEdges,1),1);

    for k = 1:size(fem.bodyEdges,1)
        shared = sum(ismember(fem.triangles, fem.bodyEdges(k,:)),2) == 2;
        fem.bodyElements(k) = find(shared,1,'first');
    end
end

% -------
function [nodes, triangles] = cylinderMesh(cfg)

    h = cfg.meshSize;
    nRing = cfg.ringPoints;

    nx = max(2, ceil(2 * cfg.Lx/h));
    ny = max(2, ceil(2 * cfg.Ly/h));

    x = linspace(0, cfg.Lx, nx+1)';
    y = linspace(0, cfg.Ly, ny+1)';

    outer = [         x(1:end-1), zeros(nx,1); ...
             cfg.Lx * ones(ny,1), y(1:end-1); ...
                flipud(x(2:end)), cfg.Ly * ones(nx,1); ...
                     zeros(ny,1), flipud(y(2:end))];

    theta = 2 * pi * (0:nRing-1)'/nRing;
    radii = cfg.ringRadii(:);
    rings = zeros(nRing * numel(radii),2);

    for k = 1:numel(radii)
        ids = (k-1) * nRing + (1:nRing);
        rings(ids,:) = [cfg.cx + radii(k) * cos(theta), ...
                        cfg.cy + radii(k) * sin(theta)];
    end

    cylinder = rings(1:nRing,:);

    [X,Y] = meshgrid((h/2:h:cfg.Lx-h/2), (h/2:h:cfg.Ly-h/2));
    background = [X(:), Y(:)];

    r = hypot((background(:,1) - cfg.cx), (background(:,2) - cfg.cy));
    background = background(r > radii(end) + 0.45 * h,:);

    hw = 0.85 * h;
    [Xw,Yw] = meshgrid(cfg.cx + radii(end) + 0.75 * h:hw:min(cfg.wakeEnd, cfg.Lx - h), ...
              max(h/2, cfg.cy - cfg.wakeHalfWidth):hw:min(cfg.Ly - h/2, cfg.cy + cfg.wakeHalfWidth));

    wake = [Xw(:),Yw(:)];
    existing = [outer; rings; background];
    keepWake = false(size(wake,1), 1);

    for k = 1:size(wake,1)
        d2 = sum((existing - wake(k,:)).^2,2);
        keepWake(k) = min(d2) > (0.55 * h)^2;
    end

    points = unique(round([existing; wake(keepWake,:)], 11),'rows');
    [~, outerId] = ismember(round(outer, 11), points,'rows');
    [~,  bodyId] = ismember(round(cylinder, 11), points,'rows');

    outerNext = [(2:numel(outerId))'; 1];
    bodyNext  = [(2:nRing)'; 1];
    constraints = [outerId, outerId(outerNext); ...
                    bodyId,  bodyId(bodyNext)];

    DT = delaunayTriangulation(points,constraints);
    nodes = DT.Points;
    triangles = DT.ConnectivityList;

    C = (nodes(triangles(:,1),:) + ...
         nodes(triangles(:,2),:) + ...
         nodes(triangles(:,3),:))/3;

    inCylinder = inpolygon(C(:,1), C(:,2), cylinder(:,1), cylinder(:,2));

    inChannel = (C(:,1) > 0) & (C(:,1) < cfg.Lx) & ...
                (C(:,2) > 0) & (C(:,2) < cfg.Ly);

    triangles = triangles(inChannel & ~inCylinder,:);


        used = unique(triangles(:));
    oldToNew = zeros(size(nodes,1),1);
    oldToNew(used) = 1:numel(used);
    
    triangles = oldToNew(triangles);
    nodes = nodes(used,:);
end

% ============================================================
% TIME INTEGRATION
% ============================================================

function result = solveFlow(fem, cfg, method, dt, snapshotTimes)

    steps = round(cfg.finalTime/dt);
    coeff = timeCoefficients(method, cfg.rhoInfinity);
    state = fem.initialState;

    result.time = (0:steps)' * dt;

    result.response = zeros(steps+1,1);
    result.response(1) = liftCoefficient(fem, cfg, state);

    result.snapTimes = snapshotTimes(:);
    result.uSnap = zeros(fem.nVelocity, numel(snapshotTimes));
    result.vSnap = zeros(fem.nVelocity, numel(snapshotTimes));
    result.newtonIterations = 0;
    result.maxResidual = 0;

    if(isempty(snapshotTimes) == false)
        result.uSnap(:,1) = state.u;
        result.vSnap(:,1) = state.v;
    end

    nextSnapshot = 2;

    for step = 1:steps
        currentTime = step * dt;

        [candidate, ok, iters, residual] = newtonStep(fem, cfg, coeff, dt, state, currentTime);
        result.newtonIterations = result.newtonIterations + iters;
 
        if(ok == false)
            error('%s: Newton solver failed at t = %.4f (residual = %.3e). Try reducing cfg.dt.', method, currentTime, residual);
        end


        state = candidate;
        result.maxResidual = max(result.maxResidual, residual);
        result.response(step + 1) = liftCoefficient(fem, cfg,state);

        while(nextSnapshot <= numel(snapshotTimes) && abs(currentTime-snapshotTimes(nextSnapshot)) < 1e-8)
            result.uSnap(:,nextSnapshot) = state.u;
            result.vSnap(:,nextSnapshot) = state.v;
            nextSnapshot = nextSnapshot + 1;
        end

        if(mod(step, max(1, round(steps/8))) == 0 || step == steps)
            w = numel(num2str(steps));
            fprintf('  %-6s | Step %*d/%d | Time: %6.2f s | Newton iterations: %2d | Residual: %9.2e\n', method, w, step, steps, currentTime, iters, residual);
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
function r = inletRamp(t,cfg)
    s = max(0, min(1, t/cfg.rampTime));
    r = s^4 * (35 - 84*s + 70*s^2 - 20*s^3);
end

% -------
function CL = liftCoefficient(fem,cfg,state)

    integrationPoints = [(0.5 - 0.5/sqrt(3)), (0.5 + 0.5/sqrt(3))];
    lift = 0;

    for k = 1:size(fem.bodyEdges,1)

        edge = fem.bodyEdges(k,:);
        p1 = fem.nodes(edge(1),:);
        p2 = fem.nodes(edge(2),:);

        len = norm(p2-p1);
        normal = [-(p2(2) - p1(2)), (p2(1) - p1(1))]/len;

        if(dot(normal, [cfg.cx, cfg.cy] - (p1 + p2)/2) < 0)
            normal = -normal;
        end

        e = fem.bodyElements(k);
        element = fem.elements(e);
        corners = fem.triangles(e,:);

        first = find(corners == edge(1),1);
        second = find(corners == edge(2),1);

        for q = 1:2

            L = zeros(1,3);
            L(first) = 1 - integrationPoints(q);
            L(second) = integrationPoints(q);

            [~, Nx, Ny] = quadraticShape(L, element.gradL);

            ue = state.u(element.nodes);
            ve = state.v(element.nodes);

            gradU = [Nx' * ue, Ny' * ue; ...
                     Nx' * ve, Ny' * ve];

            pressure = L * state.p(corners);

            traction = pressure * normal' - cfg.nu * (gradU + gradU') * normal';

            lift = lift + 0.5 * len * traction(2);
        end
    end

    CL = lift/(0.5 * cfg.Umean^2 * cfg.D);
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
 
    rhsU = -fem.M * historyU - (1-alpha) * linearU;
    rhsV = -fem.M * historyV - (1-alpha) * linearV;
    rhs = [rhsU; rhsV; zeros(m,1)];
 
    base = beta(1)/(gamma * dt) * fem.M + alpha * cfg.nu * fem.K;
    Axx = base + alpha * gradDiv * fem.Gxx;
    Axy = alpha * gradDiv * fem.Gxy;
    Ayx = Axy';
    Ayy = base + alpha * gradDiv * fem.Gyy;
    tolerance = cfg.newtonTolerance * max(1, norm(rhs(free),inf));
 
    boundaryScale = inletRamp(time, cfg);
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

            [C, Duu, Duv, Dvu, Dvv] = convection(fem, uAlpha, vAlpha);
            [~, ~, Ouu, Ovu, Ovv]   = outletTerm(fem, cfg, uAlpha, vAlpha);
 
            F = gaMomentumResidual(fem, cfg, Axx, Axy, Ayx, Ayy, C, rhs, u, v, p, previous, alpha);
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
 
            J = [Axx + alpha * (C + Duu + Ouu),              Axy + alpha * Duv,   -alpha * fem.Bx'; ...
                     Ayx + alpha * (Dvu + Ovu),  Ayy + alpha * (C + Dvv + Ovv),   -alpha * fem.By'; ...
                                alpha * fem.Bx,                 alpha * fem.By,   sparse(m,m)];
 
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
 
                FTrial = gaMomentumResidual(fem, cfg, Axx, Axy, Ayx, Ayy, ...
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
function F = gaMomentumResidual(fem, cfg, Axx, Axy, Ayx, Ayy, C, rhs, u, v, p, previous, alpha)

    n = fem.nVelocity;

    uAlpha = alpha * u + (1 - alpha) * previous.u;
    vAlpha = alpha * v + (1 - alpha) * previous.v;
    [Ru,Rv] = outletTerm(fem, cfg, uAlpha, vAlpha);

    F = [Axx * u + Axy * v + C * uAlpha + Ru - alpha * fem.Bx' * p - rhs(1:n); ...
         Ayx * u + Ayy * v + C * vAlpha + Rv - alpha * fem.By' * p - rhs(n + 1:2 * n); ...
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

% -------
function [Ru, Rv, Juu, Jvu, Jvv] = outletTerm(fem, cfg, u, v)

    n = fem.nVelocity;
    Ru = zeros(n,1);
    Rv = zeros(n,1);

    withTangent = nargout > 2;
    Juu = sparse(n,n);
    Jvu = sparse(n,n);
    Jvv = sparse(n,n);

    smoothing = cfg.backflowSmoothing * cfg.Umean;
    s = 0.5 + 0.5 * sqrt(3/5) * [-1,0,1];
    w = [5/18, 8/18, 5/18];

    for k = 1:size(fem.outletEdges, 1)

        ids = fem.outletEdges(k,:);
        len = fem.outletLength(k);

        ue = u(ids);
        ve = v(ids);

        Ju = zeros(3);
        Jx = zeros(3);
        Jy = zeros(3);

        for q = 1:numel(w)

            N = [((1 - s(q)) * (1 - 2*s(q))); (s(q) * (2 * s(q)-1)); (4 * s(q) * (1 - s(q)))];

            un = N' * ue;
            vn = N' * ve;

            g = sqrt(un^2 + smoothing^2);

            weight = 0.5 * len * w(q);

            Ru(ids) = Ru(ids) + weight * g * un * N;
            Rv(ids) = Rv(ids) + weight * g * vn * N;

            if(withTangent)
                NN = N * N';
                Ju = Ju + weight * (g + un^2/g) * NN;
                Jx = Jx + weight * (un/g) * vn * NN;
                Jy = Jy + weight * g * NN;
            end
        end

        if(withTangent)
            Juu(ids, ids) = Juu(ids, ids) + Ju;
            Jvu(ids, ids) = Jvu(ids, ids) + Jx;
            Jvv(ids, ids) = Jvv(ids, ids) + Jy;
        end
    end
end

% ============================================================
% PLOTTING
% ============================================================

function plotResponseComparison(cfg, reference, solutions, methods)

    colours = [0.75 0.26 0.23; 0.15 0.48 0.73; 0.06 0.60 0.48];

    fig = figure('Color','w','Position',[150 250 1200 500]);
     ax = axes(fig);
    hold(ax,'on');

    plot(ax,reference.time, reference.response,'k-', 'LineWidth',2.5,'DisplayName','Reference Solution');

    for k = 1:numel(methods)
        plot(ax,solutions{k}.time, solutions{k}.response, 'Color', colours(k,:), 'LineWidth', 2, 'DisplayName', methods{k});
    end

    xlabel(ax,'Physical time t');
    ylabel(ax,'Lift coefficient C_L(t)');

    title(ax, sprintf('Circular cylinder | Re = %g | Lift response', cfg.Re));
    xlim(ax, [0 cfg.finalTime]);

    legend(ax, 'show', 'Location', 'northwest');
    grid(ax,'on');
    box(ax,'on');

    exportgraphics(fig,fullfile(cfg.outputFolder,'time_response.png'), 'Resolution',300);
end

% ============================================================
% VIDEO GENERATION
% ============================================================

function createFlowAnimation(fem, cfg, solution)

    snapTimes = solution.snapTimes(:);

    ny = cfg.displayHeight;
    nx = min(1900, max(480, round(ny*cfg.Lx/cfg.Ly)));

    x = linspace(0, cfg.Lx, nx);
    y = linspace(0, cfg.Ly, ny);

    [X, Y] = meshgrid(x, y);
    [P, inside] = displayInterp(fem, X, Y);

    vort = vorticityOperator(fem, cfg);
       L = cfg.vorticityLimit;

    if(isempty(L))
        [Xc, Yc] = meshgrid(vort.x, vort.y);

        wake = vort.inside & (hypot(Xc - cfg.cx, Yc - cfg.cy) > 2*cfg.D);
           n = numel(snapTimes);
        peak = 0;

        for k = unique(round(linspace(round(n/2), n, 6)))
            w = vorticityField(vort, solution.uSnap(:,k), solution.vSnap(:,k));
            m = sort(abs(w(wake)));

            if(isempty(m) == false)
                peak = max(peak, m(max(1, round(0.99*numel(m))))); 
            end
        end
        L = max(0.5, ceil(2 * peak)/2);
    end

    fh = 173;
    fw = (fh * cfg.Lx)/cfg.Ly;
    left = max(70, (1090 - fw)/2);

    yTrace = 54;
    yVort  = 238;
    ySpeed = 447;

    fig = figure('Color'    , 'w', ...
                 'Visible'  , 'on', ...
                 'Units'    , 'pixels', ...
                 'Position' , [40 40 1200 680], ...
                 'Resize'   , 'off', ...
                 'MenuBar'  , 'none', ...
                 'ToolBar'  , 'none', ...
     'DefaultAxesFontName'  , 'Helvetica', ...
     'DefaultAxesFontSize'  , 11, ...
     'DefaultTextFontName'  , 'Helvetica');

    [~,  speedImg] = fieldPanel(fig, [left ySpeed fw fh], x, y, inside, [0 2], sequentialMap(256), {'Velocity magnitude','|u| / U_{mean}'}, cfg);
    
    [axV, vortImg] = fieldPanel(fig, [left yVort fw fh], x, y, inside, [-L L], divergingMap(L, 256), {'Vorticity','\omega D / U_{mean}'}, cfg);
    
    xlabel(axV, 'x / D');

    axT = axes(fig, 'Units', 'pixels', 'Position', [left yTrace fw 130]);
    hold(axT, 'on');
  
    plot(axT, solution.time, solution.response, 'Color',[0.78 0.80 0.84], 'LineWidth',1.4);
    trace  = plot(axT, NaN, NaN, 'Color',[0.06 0.38 0.70], 'LineWidth',1.8);
    marker = plot(axT, NaN, NaN, 'o', 'MarkerFaceColor',[0.85 0.33 0.10], 'MarkerEdgeColor','w', 'MarkerSize',7);
    
    r = solution.response(isfinite(solution.response));
    b = 1.15 * max(max(abs(r)), 0.05);

    set(axT, 'XLim',[0 solution.time(end)], 'YLim',[-b b], 'TickDir','out', 'Box','off');
    
      grid(axT, 'on');
    xlabel(axT, 'Physical time t');
    ylabel(axT, 'Lift coefficient C_L');

   heading = annotation(fig, 'textbox', [0 0.915 1 0.06], ...
                    'String'              , ''       , ...
                    'EdgeColor'           , 'none'   , ...
                    'HorizontalAlignment' , 'center' , ...
                    'VerticalAlignment'   , 'middle' , ...
                    'FontWeight'          , 'bold'   , ...
                    'FontSize'            , 14       , ...
                    'Interpreter'         , 'none');

    videoTimes = 0:cfg.videoSpeed/cfg.videoFrameRate:snapTimes(end);

    if(videoTimes(end) < snapTimes(end)-1e-9)
        videoTimes(end+1) = snapTimes(end);
    end

    writer = VideoWriter(fullfile(cfg.outputFolder,'reference_flow.mp4'), 'MPEG-4');
    writer.FrameRate = cfg.videoFrameRate;
    writer.Quality = 100;
    open(writer);

    printOption = sprintf('-r%d', max(80, round(cfg.videoWidth/1200*get(groot,'ScreenPixelsPerInch'))));

    nF = numel(videoTimes);

    for f = 1:nF

        t = videoTimes(f);
        u = interpSnapshot(snapTimes, solution.uSnap, t);
        v = interpSnapshot(snapTimes, solution.vSnap, t);

        speed = hypot(P * u, P * v)/cfg.Umean;
        omega = interp2(vort.x, vort.y, vorticityField(vort,u,v), X, Y, 'linear');
            j = find(solution.time <= t + 1e-9, 1, 'last');

        set(speedImg, 'CData', reshape(speed, ny, nx));
        set(vortImg,  'CData', reshape(omega, ny, nx));

        set(trace,  'XData', solution.time(1:j), 'YData', solution.response(1:j));
        set(marker, 'XData', solution.time(j),   'YData', solution.response(j));
        
        heading.String = sprintf('Circular cylinder | Re = %g | t = %.2f (t U/D = %.1f)', cfg.Re, t, t * cfg.Umean/cfg.D);
        drawnow;

        img = print(fig, '-RGBImage', printOption);
        writeVideo(writer, img(1:2 * floor(end/2), 1:2 * floor(end/2), :));
    end

    close(writer);
    close(fig);
end

% -------
function [ax, img] = fieldPanel(fig, pos, x, y, inside, limits, map, barLabel, cfg)
     
     th = linspace(0, 2*pi, 361);
     ax = axes(fig, 'Units','pixels', 'Position',pos);
    img = imagesc(ax, x/cfg.D, y/cfg.D, zeros(numel(y),numel(x)), 'AlphaData', double(inside));
    hold(ax, 'on');

    fill(ax, (cfg.cx + cfg.R * cos(th))/cfg.D, (cfg.cy + cfg.R * sin(th))/cfg.D, [0.88 0.89 0.91], 'EdgeColor', [0.12 0.12 0.14], 'LineWidth',1.2);
    rectangle(ax, 'Position', [0 0 cfg.Lx/cfg.D cfg.Ly/cfg.D], 'EdgeColor', [0.25 0.25 0.28], 'LineWidth',0.8);

    set(ax, 'YDir', 'normal', ...
 'DataAspectRatio',  [1 1 1], ...
            'XLim',  [0 cfg.Lx/cfg.D], ...
            'YLim',  [0 cfg.Ly/cfg.D], ...
           'XTick', 0:2:floor(cfg.Lx/cfg.D), ...
           'YTick', 0:floor(cfg.Ly/cfg.D), ...
         'TickDir', 'out', ...
           'Layer', 'top', ...
             'Box', 'off', ...
       'LineWidth',  0.8);
    
    clim(ax, limits);
    colormap(ax, map);
    bar = colorbar(ax);

    set(bar, 'Units', 'pixels', 'Position', [pos(1) + (pos(3) + 16) pos(2) 14 pos(4)], 'TickDirection','out', 'LineWidth',0.8);
    ylabel(bar, barLabel);
    set(ax, 'Units','pixels', 'Position',pos);
    ylabel(ax, 'y / D');
end

% -------
function vort = vorticityOperator(fem, cfg)

    ny = max(60, round(cfg.displayHeight/3));
    nx = round(ny * cfg.Lx/cfg.Ly);

    vort.x = linspace(0, cfg.Lx, nx);
    vort.y = linspace(0, cfg.Ly, ny);

    vort.dx = vort.x(2) - vort.x(1);
    vort.dy = vort.y(2) - vort.y(1);

    [X, Y] = meshgrid(vort.x, vort.y);
    [vort.P, vort.inside] = displayInterp(fem, X, Y);

    vort.mask  = double(vort.inside);
    vort.size  = [ny, nx];
    vort.scale = cfg.D/cfg.Umean;

    sigma = max(cfg.vorticitySigma * cfg.meshSize/vort.dx, 0.5);
    k = exp(-0.5 * ((-ceil(3 * sigma):ceil(3 * sigma))/sigma).^2);

    vort.kernel = k/sum(k);
    vort.weight = max(conv2(vort.kernel', vort.kernel, vort.mask, 'same'), 1e-3);
end

% -------
function omega = vorticityField(vort, u, v)

    uGrid = reshape(vort.P * u, vort.size);
    vGrid = reshape(vort.P * v, vort.size);

    uSmooth = conv2(vort.kernel', vort.kernel, vort.mask .* uGrid, 'same') ./ vort.weight;
    vSmooth = conv2(vort.kernel', vort.kernel, vort.mask .* vGrid, 'same') ./ vort.weight;

    [dvdx, ~] = gradient(vSmooth, vort.dx, vort.dy);
    [~, dudy] = gradient(uSmooth, vort.dx, vort.dy);

    omega = (dvdx - dudy) * vort.scale;

end

% -------
function f = interpSnapshot(times, snaps, t)
    
    if(t >= times(end))
        f = snaps(:,end);
        return; 
    end

    i = find(times <= t, 1, 'last');
    a = (t - times(i))/(times(i+1) - times(i));
    f = (1 - a) * snaps(:,i) + a * snaps(:,i+1);
end

% -------
function [P, inside] = displayInterp(fem, X, Y)

    mesh = triangulation(fem.triangles, fem.nodes);
     pts = [X(:), Y(:)];
      el = pointLocation(mesh, pts);
      
     valid = find(~isnan(el));
    inside = reshape(~isnan(el), size(X));

    B = cartesianToBarycentric(mesh, el(valid), pts(valid,:));
    a = B(:,1);
    b = B(:,2);
    c = B(:,3);
    N = [a .* (2 * a-1), b .* (2 * b-1), c .* (2 * c-1), 4 * a .* b, 4 * b .* c, 4 * c .* a];

    cols = fem.velocityTriangles(el(valid),:);
       P = sparse(repelem(valid,6), reshape(cols.',[],1), reshape(N.',[],1), numel(X), fem.nVelocity);
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
function map = divergingMap(L, n)
    blue = [1 1 1; 0.573 0.773 0.871; 0.263 0.576 0.765; 0.129 0.400 0.675; 0.020 0.188 0.380];
    red  = [1 1 1; 0.957 0.647 0.510; 0.839 0.376 0.302; 0.698 0.094 0.168; 0.404 0 0.122];

    a = linspace(0, 1, 5);
    v = linspace(-L, L, n)';

    map = ones(n, 3);
    map(v < 0,:) = interp1(a, blue, min(1, -v(v < 0)/L));
    map(v > 0,:) = interp1(a, red,  min(1,  v(v > 0)/L));
end