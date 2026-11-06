
clc; clear; close all;

% -------------------------------------------------
% 2D flow over a backward-facing step with P2/P1 finite
% -------------------------------------------------

cfg.Lx = 40;
cfg.Ly = 3;

cfg.stepX = 4;
cfg.stepH = 1;
cfg.H = cfg.stepH;

cfg.nyBelow = 8;
cfg.nyAbove = 14;
cfg.nxUpstream = 8;
cfg.nxDownstream = 72;

cfg.gradeY = 2.0;
cfg.gradeUpstream = 1.5;
cfg.gradeDownstream = 2.5;

cfg.Re = 229;
cfg.Umean = 1.0;
cfg.nu = cfg.Umean * cfg.H / cfg.Re;

cfg.dt = 0.5;
cfg.rampTime = 5;
cfg.finalTime = 150;
cfg.rhoInfinity = 0;

cfg.gradDiv = 0.03;

cfg.newtonTolerance = 1e-9;
cfg.maxNewtonIterations = 22;

cfg.pressureNode = [];
cfg.pressureValue = 0;

cfg.backflowSmoothing = 0.01;

cfg.stations = [0 4 4.88 6.11 8.17 14.29];
cfg.nearWall = 0.025;

cfg.outputFolder = 'BackwardFacingStep';

cfg.saveVideo = true;
cfg.videoFrameRate = 30;
cfg.videoSpeed = 10;
cfg.videoTime = 120;
cfg.videoWidth = 1920;
cfg.displayHeight = 240;

if(isfolder(cfg.outputFolder) == false)
    mkdir(cfg.outputFolder);
end

% ------------------------------------------
% Finite Element Mesh and Initial Conditions
% ------------------------------------------

fprintf('\n-------------------------------------------------------');
fprintf('\nBackward-facing step | P2/P1 finite element simulation\n');
fprintf('--------------------------------------------------------\n');

fem = buildMesh(cfg);
fem = setupBoundaryDOFs(fem,cfg);

data = benchmarkData();

xu = fem.velocityNodes(:,1);
yu = fem.velocityNodes(:,2);

inlet = (abs(xu) < 1e-8);

inletU = zeros(fem.nVelocity,1);
inletU(inlet) = interp1(data.inletY, data.inletU, yu(inlet), 'pchip');

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
fem = setupConvectionData(fem);

% ------------------------------------------
% Solution
% ------------------------------------------
snapshotTimes = (0:cfg.dt:cfg.finalTime)';

fprintf('\nSolution (GA-234, rho Infinity = %.2f) | dt = %.4f\n', cfg.rhoInfinity, cfg.dt);
solution = solveFlow(fem, cfg, 'GA-234', cfg.dt, snapshotTimes);

% ------------------------------------------
% Plotting
% ------------------------------------------

plotProfiles(fem, cfg, solution, data);

if(cfg.saveVideo); createFlowAnimation(fem, cfg, solution); end

fprintf('\nSimulation completed. Results saved to: %s\n', cfg.outputFolder);

% ============================================================
% FINITE ELEMENT SETUP
% ============================================================

function fem = buildMesh(cfg)

    [nodes, triangles] = stepMesh(cfg);

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
function fem = setupConvectionData(fem)

    [points, weights] = triangleQuadrature();

    nq = numel(weights);
    ne = fem.nElements;

    gradients = cat(3, fem.elements.gradL);
    gradX = squeeze(gradients(1,:,:))';
    gradY = squeeze(gradients(2,:,:))';
    area  = [fem.elements.area]';

    shapeValue = zeros(nq, 6);
    shapeX     = zeros(ne, 6, nq);
    shapeY     = zeros(ne, 6, nq);

    for q = 1:nq
        
        [N, d1] = quadraticShape(points(q,:), [1 0 0; 0 0 0]);
        [~, d2] = quadraticShape(points(q,:), [0 1 0; 0 0 0]);
        [~, d3] = quadraticShape(points(q,:), [0 0 1; 0 0 0]);

        dNdL = [d1, d2, d3];

        shapeValue(q,:) = N';
        shapeX(:,:,q)   = gradX * dNdL';
        shapeY(:,:,q)   = gradY * dNdL';
    end

    ids = fem.velocityTriangles;

    rowIndex = repmat(1:6, 1, 6);
    colIndex = repelem(1:6, 6);

    fem.conv.ids   = ids;
    fem.conv.shape = shapeValue;

    fem.conv.shapeX = shapeX;
    fem.conv.shapeY = shapeY;
    fem.conv.weight = area * weights';

    fem.conv.rowIndex = rowIndex;
    fem.conv.colIndex = colIndex;

    fem.conv.rows = ids(:, rowIndex);
    fem.conv.cols = ids(:, colIndex);
end

% -------
function [nodes, triangles] = stepMesh(cfg)

    xi = linspace(0, 1, cfg.nxUpstream + 1);
    xUpstream = sort(cfg.stepX - cfg.stepX * ((exp(cfg.gradeUpstream * xi) - 1)/(exp(cfg.gradeUpstream) - 1)));

    xi = linspace(0, 1, cfg.nxDownstream + 1);
    xDownstream = cfg.stepX + (cfg.Lx - cfg.stepX) * ((exp(cfg.gradeDownstream * xi) - 1)/(exp(cfg.gradeDownstream) - 1));

    xx = [xUpstream(1:end-1), xDownstream];

    eta = linspace(0, 1, cfg.nyBelow + 1);
    yBelow = cfg.stepH * 0.5 * (1 + tanh(cfg.gradeY * (eta - 0.5))/tanh(cfg.gradeY/2));

    eta = linspace(0, 1, cfg.nyAbove + 1);
    yAbove = cfg.stepH + (cfg.Ly - cfg.stepH) * 0.5 * (1 + tanh(cfg.gradeY * (eta - 0.5))/tanh(cfg.gradeY/2));

    yy = [yBelow, yAbove(2:end)];

    [X, Y] = meshgrid(xx, yy);
    nodes = [X(:), Y(:)];

    nx = numel(xx) - 1;
    ny = numel(yy) - 1;

    triangles = zeros(2 * nx * ny, 3);
    row = 0;

    for i = 1:nx
        for j = 1:ny

            xm = (xx(i) + xx(i+1))/2;
            ym = (yy(j) + yy(j+1))/2;

            if(xm < cfg.stepX && ym < cfg.stepH); continue; end

            a = (i-1) * (ny+1) + j;
            b = i * (ny+1) + j;
            c = b + 1;
            d = a + 1;

            row = row + 1;
            triangles(row,:) = [a, b, c];

            row = row + 1;
            triangles(row,:) = [a, c, d];
        end
    end

    triangles = triangles(1:row,:);

    used = unique(triangles(:));
    oldToNew = zeros(size(nodes,1),1);
    oldToNew(used) = 1:numel(used);

    triangles = oldToNew(triangles);
    nodes = nodes(used,:);
end


% ============================================================
% BENCHMARK DATA
% ============================================================

function data = benchmarkData()

    data.inletY = [ ...
        1.0000000 1.0346534 1.0769792 1.1286761 1.1918188 1.2689414 ...
        1.3631392 1.4781927 1.6187193 1.7903589 2.0000000 2.2096411 ...
        2.3812807 2.5218073 2.6368608 2.7310586 2.8081812 2.8713239 ...
        2.9230208 2.9653466 3.0000000]';

    data.inletU = [ ...
        0.000000000 0.159516590 0.343900880 0.545084650 0.755216410 0.961046240 ...
        1.143818200 1.282067700 1.359530700 1.377773500 1.361100000 1.325760600 ...
        1.248506500 1.113912800 0.935216310 0.737315800 0.542620330 0.365704240 ...
        0.213463690 0.087295732 0.000000000]';

    % ------------------------------------------
    % Denham & Patrick (1974) Experimental Data
    % ------------------------------------------

    % Station 1: x/H = 0
    data.denham(1).u = [ ...
        0.00000000 0.10271492 0.29251403 0.56716419 0.77482713 1.11646520 ...
        1.26607150 1.35762200 1.38218400 1.37995080 1.35538890 1.35092330 ...
        1.30403180 1.21918000 1.06510800 0.87084263 0.56046537 0.37289942 ...
        0.16970275 0.08038541 0.00000000]';

    data.denham(1).y = [ ...
        1.0000000 1.0211465 1.0694832 1.1299042 1.2054291 1.3262702 ...
        1.4863849 1.6102464 1.7431719 1.8881808 2.0241268 2.1479893 ...
        2.2778925 2.4229023 2.5437426 2.6796895 2.8035509 2.8730341 ...
        2.9274133 2.9636656 2.9999179]';

    % Station 2: x/H = 4
    data.denham(2).u = [ ...
        0.00000000 0.19724855 0.57844739 0.79564233 1.13473180 1.29208740 ...
        1.37630550 1.38295440 1.38517090 1.36965670 1.32754800 1.27214100 ...
        1.18792230 1.04386490 0.82445344 0.50974277 0.35238774 0.14627387 ...
        0.00443232]';

    data.denham(2).y = [ ...
        1.0000000 1.0479760 1.1499248 1.2188896 1.3538214 1.4887532 ...
        1.6236859 1.7496223 1.8755586 2.0134891 2.1514197 2.2833527 ...
        2.4182846 2.5532163 2.6881489 2.8140853 2.8770535 2.9460183 ...
        2.9999910]';

    % Station 3: x/H = 4.88
    data.denham(3).u = [ ...
        -0.01542870 -0.03967493 -0.03747055 -0.05510364 -0.03526617  0.02204185 ...
         0.17853786  0.45185429  0.84639990  1.11089950  1.35556220  1.35776600 ...
         1.34013290  1.30266240  1.20127070  1.00950780  0.49593735  0.32401263 ...
         0.14988353  0.05069619  0.00220438  0.00000000]';

    data.denham(3).y = [ ...
        0.06262416 0.14314094 0.21471141 0.42047650 0.68588428 0.81113257 ...
        0.94830963 1.08250430 1.21968130 1.35089440 1.55665960 1.74453200 ...
        1.95924340 2.16799100 2.35884570 2.53479050 2.80914350 2.87176780 ...
        2.93141040 2.96719580 2.99999850 0.00000000]';

    % Station 4: x/H = 6.11
    data.denham(4).u = [ ...
         0.00000013 -0.03326069 -0.05099951 -0.09312998 -0.06430374 -0.04434715 ...
         0.05986942  0.31930460  0.61643423  0.91134752  1.17078170  1.33265170 ...
         1.37256440  1.27721660  1.14417330  0.90026011  0.68517365  0.38804306 ...
         0.25056519  0.09756483  0.00221790]';

    data.denham(4).y = [ ...
        0.000000000 0.051000509 0.129000050 0.335999700 0.525000220 0.602999810 ...
        0.744000200 0.944999220 1.082999200 1.220999200 1.355999000 1.547999000 ...
        1.778999200 2.150998800 2.345998300 2.546998000 2.678998500 2.813998200 ...
        2.876997900 2.933998300 2.999998600]';

    % Station 5: x/H = 8.17
    data.denham(5).u = [ ...
         0.00000008 -0.06860153 -0.09737103 -0.10622303 -0.10843555 -0.11064808 ...
        -0.07966702  0.04868515  0.12613876  0.30538853  0.52668484  0.80330570 ...
         1.05115800  1.22155580  1.33441740  1.25253780  1.21491730  0.96485241 ...
         0.72806461  0.52225979  0.27883347  0.16818531  0.06196220  0.00221260]';

    data.denham(5).y = [ ...
        0.000000000 0.071856163 0.131736010 0.200598760 0.272454920 0.353293120 ...
        0.398203220 0.610778330 0.679640110 0.808382750 0.952095090 1.104789400 ...
        1.233531900 1.344310500 1.748501500 1.946106900 2.146705600 2.344310000 ...
        2.535926800 2.676644800 2.808381800 2.883232100 2.946106400 2.999998600]';

    % Station 6: x/H = 14.29
    data.denham(6).u = [ ...
        0.0000000 0.0728826 0.1258878 0.4019575 0.6603574 1.0247697 ...
        1.2036619 1.2036619 1.0976514 0.8790044 0.3334913 0.2076034 ...
        0.0927591 0.0463800 0.0176677 0.0022096]';

    data.denham(6).y = [ ...
        0.00000000 0.18525925 0.25398397 0.53187287 0.74103618 0.93227035 ...
        1.11752950 1.52091620 1.92131380 2.12748980 2.53685140 2.67131400 ...
        2.78784780 2.85956070 2.94621370 2.99999860]';

end


% ============================================================
% TIME INTEGRATION
% ============================================================

function result = solveFlow(fem, cfg, method, dt, snapshotTimes)

    steps = round(cfg.finalTime/dt);
    coeff = timeCoefficients(method, cfg.rhoInfinity);
    state = fem.initialState;

    result.time = (0:steps)' * dt;

    result.snapTimes = snapshotTimes(:);
    result.uSnap = zeros(fem.nVelocity, numel(snapshotTimes));
    result.vSnap = zeros(fem.nVelocity, numel(snapshotTimes));
    result.pSnap = zeros(fem.nPressure, numel(snapshotTimes));

    result.newtonIterations = 0;
    result.maxResidual = 0;

    if(isempty(snapshotTimes) == false)
        result.uSnap(:,1) = state.u;
        result.vSnap(:,1) = state.v;
        result.pSnap(:,1) = state.p;
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

        while(nextSnapshot <= numel(snapshotTimes) && abs(currentTime-snapshotTimes(nextSnapshot)) < 1e-8)
            result.uSnap(:,nextSnapshot) = state.u;
            result.vSnap(:,nextSnapshot) = state.v;
            result.pSnap(:,nextSnapshot) = state.p;
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

function plotProfiles(fem, cfg, solution, data)

    blue = [0.00 0.45 0.74];
    red  = [0.85 0.20 0.20];
    gray = [0.55 0.55 0.55];

    fig = figure('Color', 'w', 'Position',[80 80 1300 600]);
    ax = axes(fig, 'Position', [0.09 0.13 0.86 0.73]);
    hold(ax, 'on');

    plot(ax, [cfg.stepX cfg.stepX], [0 cfg.stepH], '-', 'Color', gray, 'LineWidth', 1.5);
    plot(ax ,[cfg.stepX 17], [0 0], '-', 'Color', gray, 'LineWidth', 1.2);

    xNodes = unique(fem.nodes(:,1));

    for k = 1:numel(cfg.stations)

        s = cfg.stations(k);

        [~,i] = min(abs(xNodes - s));
        x = xNodes(i);

        ids = find(abs(fem.velocityNodes(:, 1) - x) < 1e-8);

        if(s <= cfg.stepX)
            ids = ids(fem.velocityNodes(ids, 2) >= cfg.stepH);
        end

        [y,order] = sort(fem.velocityNodes(ids, 2));
        uFEM = solution.u(ids(order));

        hf = plot(ax,s + uFEM,y, '-', 'Color', blue, 'LineWidth', 2);
        hd = plot(ax,s + data.denham(k).u(:), data.denham(k).y(:), 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', red,'LineWidth', 1.2);

    end

    xlabel(ax,'x / H', 'FontSize', 12);
    ylabel(ax,'y / H', 'FontSize', 12);

    xlim(ax,[0 17]);
    ylim(ax,[0 cfg.Ly]);

    set(ax,'XTick', 0:2:16, 'YTick', 0:0.5:cfg.Ly);
    set(ax,'FontName', 'Arial','FontSize', 11, 'LineWidth', 1);
    set(ax,'TickDir', 'out', 'GridAlpha', 0.12);
    grid(ax,'on');

    title(ax,sprintf('Backward-facing step | Re = %g | GA-234',cfg.Re),'FontSize',14,'FontWeight','bold');

    ax.Title.Units = 'normalized';
    ax.Title.Position = [0.5 1.03 0];

    lgd = legend(ax, [hf hd], {'GA-234 FEM (P2/P1)','Denham & Patrick (1974)'}, 'Location', 'southwest');
    set(lgd,'FontName','Arial','FontSize', 11, 'Color', 'white', 'Box', 'on');
    lgd.AutoUpdate = 'off';

    exportgraphics(fig,fullfile(cfg.outputFolder,'profiles.png'),'Resolution',300);
end


% ============================================================
% VIDEO GENERATION
% ============================================================

function createFlowAnimation(fem, cfg, solution)

    times = solution.snapTimes(:);

    ny = cfg.displayHeight;
    nx = min(1900, max(480, round(ny * cfg.Lx/cfg.Ly)));

    x = linspace(0, cfg.Lx, nx);
    y = linspace(0, cfg.Ly, ny);

    [X, Y] = meshgrid(x, y);
    [P, inside] = displayInterp(fem, X, Y);

    Q = pressureInterp(fem, X, Y);

    pLimits = [min(solution.pSnap(:, end)) max(solution.pSnap(:, end))];
    uLimits = [min(solution.uSnap(:, end)) max(solution.uSnap(:, end))];
    vLimit  = max(0.5 * max(abs(solution.vSnap(:))), 1e-6);

    if(diff(uLimits) < 1e-8); uLimits = uLimits + [-0.5 0.5]; end
    if(diff(pLimits) < 1e-8); pLimits = pLimits + [-0.5 0.5]; end
    pLevels = linspace(pLimits(1), pLimits(2), 26);

    figW = 1200;
    figH = 580;
    left = 75;

    fw = 1005;
    fh = 72;
    yP = 75;
    yV = 240;
    yU = 405;

    fig = figure('Color','w','Visible','on','Units','pixels','Position',[40 40 figW figH], 'Resize','off','MenuBar','none','ToolBar','none');

    [axU, uImg] = fieldPanel(fig,[left yU fw fh],x,y,inside,uLimits,jet(256),cfg);
    [axV, vImg] = fieldPanel(fig,[left yV fw fh],x,y,inside,[-vLimit vLimit],jet(256),cfg);
    [axP, pImg] = fieldPanel(fig,[left yP fw fh],x,y,inside,pLimits,jet(256),cfg);

    set(axU,'XTick',[]);
    set(axV,'XTick',[]);
    xlabel(axP,'x / H','FontSize',11);

    panelTitle(fig,'Horizontal velocity', left, (yU + fh + 12), fw, figW, figH);
    panelTitle(fig,'Vertical velocity', left, (yV + fh + 12), fw, figW, figH);
    panelTitle(fig,'Pressure', left, (yP + fh + 12), fw, figW, figH);

    heading = annotation(fig,'textbox',  [0.06 0.922 0.88 0.055], ...
                              'String',  '', ...
                           'EdgeColor',  'none', ...
                 'HorizontalAlignment',  'center', ...
                   'VerticalAlignment',  'middle', ...
                            'FontName',  'Arial', ...
                          'FontWeight',  'bold', ...
                            'FontSize',   14, ...
                         'Interpreter',  'none');

    tEnd = min(times(end), 100);
    videoTimes = 0:cfg.videoSpeed/cfg.videoFrameRate:tEnd;
    if(videoTimes(end) < tEnd - 1e-9); videoTimes(end + 1) = tEnd; end

    writer = VideoWriter(fullfile(cfg.outputFolder,'BackwardFacingStep.mp4'),'MPEG-4');
    writer.FrameRate = cfg.videoFrameRate;
    writer.Quality = 100;
    open(writer);

    dpi = max(150, round(cfg.videoWidth/figW * get(groot,'ScreenPixelsPerInch')));
    printOption  = sprintf('-r%d',dpi);
    contourLines = gobjects(0);

    for f = 1:numel(videoTimes)

        t = videoTimes(f);
        u = interpSnapshot(times,solution.uSnap, t);
        v = interpSnapshot(times,solution.vSnap, t);
        p = interpSnapshot(times,solution.pSnap, t);

        uGrid = reshape(P*u, ny, nx);
        vGrid = reshape(P*v, ny, nx);
        pGrid = reshape(Q*p, ny, nx);

        set(uImg,'CData', uGrid);
        set(vImg,'CData', vGrid);
        set(pImg,'CData', pGrid);

        if(~isempty(contourLines))
            delete(contourLines(isgraphics(contourLines)));
        end

        contourLines = gobjects(0);
        maskedP = pGrid;
        maskedP(~inside) = NaN;
        if(any(isfinite(maskedP(:))) && max(maskedP(:),[],'omitnan') - min(maskedP(:),[],'omitnan') > 1e-10)
            [~,contourLines] = contour(axP, X/cfg.H, Y/cfg.H, maskedP, pLevels, 'LineColor', 'w', 'LineWidth', 0.65);
        end

        heading.String = sprintf('Backward-facing step | Re = %g | GA-234 | t = %.1f', cfg.Re, t);
        drawnow;
        frame = print(fig, '-RGBImage', printOption);
        frame = frame(1:2*floor(size(frame,1)/2), 1:2*floor(size(frame,2)/2), :);
        writeVideo(writer,frame);
    end

    close(writer);
    close(fig);
end

% -------
function panelTitle(fig, label, left, bottom, width, figW, figH)
    
annotation(fig,'textbox', [left/figW bottom/figH width/figW 23/figH], ...
                'String', label, ...
             'EdgeColor', 'none', ...
   'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'middle', ...
              'FontSize', 12, ...
              'FontName', 'Arial', ...
            'FontWeight', 'bold', ...
           'Interpreter', 'none');
end

% -------
function [ax,img] = fieldPanel(fig, pos, x, y, inside, limits, map, cfg)

     ax = axes(fig, 'Units', 'pixels', 'Position', pos);
    img = imagesc(ax, x/cfg.H, y/cfg.H, zeros(numel(y),numel(x)), 'AlphaData', double(inside));
    hold(ax,'on');

    fill(ax,[0 cfg.stepX cfg.stepX 0]/cfg.H, [0 0 cfg.stepH cfg.stepH]/cfg.H, [0.87 0.89 0.92], 'EdgeColor', 'none'); 
    plot(ax,[0 cfg.stepX cfg.stepX cfg.Lx cfg.Lx 0 0]/cfg.H, [cfg.stepH cfg.stepH 0 0 cfg.Ly cfg.Ly cfg.stepH]/cfg.H, '-', 'Color', [0.23 0.23 0.23], 'LineWidth', 0.85);

    set(ax, 'YDir', 'normal', 'XLim', [0 cfg.Lx/cfg.H], 'YLim', [0 cfg.Ly/cfg.H], 'XTick', 0:4:floor(cfg.Lx/cfg.H), 'YTick', 0:floor(cfg.Ly/cfg.H), 'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', 'Box', 'on', 'LineWidth', 0.8);
    ylabel(ax,'y / H');
    
    colormap(ax, map);
    clim(ax, limits);
    
    cb = colorbar(ax);
    set(cb, 'Units', 'pixels', 'Position', [(pos(1) + pos(3) + 16) pos(2) 14 pos(4)], 'TickDirection', 'out', 'LineWidth', 0.8, 'FontSize', 9);
    set(ax, 'Units', 'pixels', 'Position', pos);
end

% -------
function value = interpSnapshot(times, snapshots,t)

    if(t >= times(end))
        value = snapshots(:,end);
        return;
    end

    k = find(times <= t, 1, 'last');
    a = (t - times(k))/(times(k + 1) - times(k));
    value = (1-a) * snapshots(:,k) + a * snapshots(:,k + 1);
end

% -------
function [P,inside] = displayInterp(fem, X, Y)

    mesh = triangulation(fem.triangles, fem.nodes);
    pts = [X(:) Y(:)];
    el  = pointLocation(mesh, pts);

    valid  = find(~isnan(el));
    inside = reshape(~isnan(el), size(X));

    B = cartesianToBarycentric(mesh, el(valid), pts(valid,:));
    a = B(:,1);
    b = B(:,2);
    c = B(:,3);
    N = [(a .* (2 * a-1)) (b .* (2 * b-1)) (c .* (2 * c-1)) (4 * a .* b) (4 * b .* c) (4 * c .* a)];

    ids = fem.velocityTriangles(el(valid), :);
      P = sparse(repelem(valid, 6), reshape(ids.',[], 1), reshape(N.', [], 1), numel(X), fem.nVelocity);
end

% -------
function Q = pressureInterp(fem, X, Y)

    mesh = triangulation(fem.triangles, fem.nodes);
    pts = [X(:) Y(:)];

    el = pointLocation(mesh, pts);
    valid = find(~isnan(el));

    B = cartesianToBarycentric(mesh, el(valid), pts(valid,:));
  ids = fem.triangles(el(valid), :);
    Q = sparse(repelem(valid,3), reshape(ids.',[],1), reshape(B.',[],1), numel(X), fem.nPressure);
end