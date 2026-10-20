
clc; clear; close all;

% -------------------------------------------------
% 2D Taylor-Green vortex with P2/P1 finite elements
% -------------------------------------------------

% Domain and mesh
cfg.L = 2 * pi;
cfg.meshSize = 0.1;

% Physical parameters
cfg.nu = 0.1;

% Time integration
cfg.finalTime = 10;
cfg.dt = 0.25;
cfg.rhoInfinity = 0;

% Numerical stabilisation
cfg.gradDiv = 0.03;

% Newton solver
cfg.newtonTolerance = 1e-9;
cfg.maxNewtonIterations = 22;

% Taylor-Green decay
cfg.velocityDecay = -2 * cfg.nu;

% Pressure reference
cfg.pressureNode = 1;
cfg.pressureValue = 0;

% Output settings
cfg.outputFolder = 'ClassicalTaylorGreen';

% Video settings
cfg.saveCombinedVideo = true;
cfg.videoFrameRate = 24;
cfg.tracerVideoSeconds = 12;

% Tracer settings
cfg.particlesPerSide = 50;
cfg.tracerGridSize = 220;
cfg.tracerTrailFrames = 30;

% Create output folder
if (isfolder(cfg.outputFolder) == false)
    mkdir(cfg.outputFolder);
end
