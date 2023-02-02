%{ 
...
Copyright Nikita
A fully automated CR3BP library

REFERENCES FOR ALGORITHMS(Any one will do as long as you have no doubt)
--------------------------------------------------------------------------
  1) W S koon ,MW Lo, JE Marsden, SD Ross - Dynamical Systems, the Three Body Problem and Space Mission Design 
    (Main refrence for almost everything)
  2) Daniel J Grebow - Generating Periodic Orbits in the circular
     Restricted Three Body Problem With Applications to Lunar South Pole
     Coverage - MS Thesis (2006)(Insights for differential Correction and
     continuation Halo orbits)
  3) Srianish Vutukuri - Spacecraft Trajectory Design Techniques using
     Resonant Orbits (for DC and continuation)
  4) Thomas Pavlak - Mission Design In the Earth Moon System Transfer
     Trajectories and station keeping
  5) TM Vaquero Escribano - Poincare Sections and resonant orbits in the
     Restricted Three Body Problem.

Library Functions
  1) calcEigen
  2) calcStabilityIdx
  3) continuationBifurcation
  4) CR3BP
  5) DFmatrix3D
  6) diffCorrec
  7) equilPts
  8) event_yCrossing
  9) findBifurcationPt
  10) getGlobalVariable
  12) haloFamily
  13) haloOrbit
  14) initGuess
  15) integrate
  16) jacobiValue3D
  17) lyapunovFamily
  18) lyapunovOrbit
  19) orbitInvManifoldIC
  20) plotBifucationSln
  21) plotFamily
  22) plotOrbit
  23) pdeudoArcDiffCorrec
  24) pseudoArcLengthCont
  25) richardson3Halo
  26) stm_X
  27) twoLevelDiffCorrec
  28) varEq_initialization
  29) varEq_stmDot

See Function descriptions for more details
%}

%% --------------------------------- Mainline Program ----------------------------------------------------------------
clearvars;clc;close all

%% User Input
userInput.Dimension         = 3;   % 2/3
userInput.mu                = 0.0121505856; % system Parameter
userInput.lagrangePt        = 1; % lagrange Point
userInput.orbitCount        = 150; % No. of Orbits in the family
userInput.plotDiffCorrec    = 1;   % Is Differential Correction Plotted (0/1)
userInput.orbit             = 'halo';  % which Orbit
userInput.type              = 'northern';  % 'northern'/'southern' for 'halo' else 'none'
userInput.tolerance         = 1e-6;    % solution tolerance
userInput.jacobianConst     = 3.1; % Required jacobian Constant
userInput.manifoldType      = 'unstable';
userInput.diffCorrecPlot    = figure('Name','Differential Correction','NumberTitle','off');
userInput.orbitPlot         = figure('Name','Orbit Plot','NumberTitle','off');
userInput.manifoldPlot      = figure('Name','Manifold Plot','NumberTitle','off');
tic

globalVar                   = getGlobalVariable(userInput);
fprintf('System /mu value %f\n',globalVar.userInput.mu);

%% Input Orbit Code here

familyPar = haloFamily(globalVar);
figure()
plot(familyPar.IC(:,1), familyPar.IC(:,3));

%plotFamily(globalVar);
%plotBifurcationSln(globalVar);
toc