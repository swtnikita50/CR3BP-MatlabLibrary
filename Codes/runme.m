%{ 
...
Copyright Nikita

A fresh Modification to reduce the number of varibles in workspace is started on 20
Feb,2020

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



File dependencies (Explicit -E)
  1) GlobalData

File dependencies (Implicit - I)
  1) diffCorrec             - will call E(1) I(3) I(4) I(5)
  2) integrate             - will call I(5) 
  3) stm_X         - willcall I(6),I(3) and I(3) will call (7)
  4) CR3BP
  5) varEq_initialization 
  6) varEq_stmDot         - will call E(1) I(8) I(5) 
  7) Dfmatrix3D
  8) initGuess

See Function descriptions for more details

****************************************************************
                    Read before proceeding
****************************************************************

NOTE:
-----
1) See the descriptions of each functions inside respective function files(if required).
2) Some figures might be small, try to zoom and see 
3) Eq.Point Manifolds might require different time of integration depending on mu
   value. Change it and see 
4) Except the "LyapunovPlotter" all the plots requires only global
   variables in workspace
5) "LyapunovPlotter" requires saved parameters of the orbits (Changing the
    file name might require you to change the variable name inside the plotter.
    This is made for convience because once the parameters are saved for a
    particular primary and secondary, one can extract and plot the required
    data any number of times without having to worry to run whole script again.
6) ENTER ONLY "userInput" asked for . This is fully sufficient.
7) To know the flow see the file depencies above(if required).

...
%}

%% --------------------------------- Mainline Program ----------------------------------------------------------------
clearvars;clc;close all

%% User Input
userInput.Dimension      = 3;   % 2/3
userInput.mu        = 0.0121505856; % system Parameter
userInput.lagrangePt       = 1; % lagrange Point
userInput.orbitCount  = 10; % No. of Orbits in the family
userInput.plotDiffCorrec = 1;   % Is Differential Correction Plotted (0/1)
userInput.orbit          = 'halo';  % which Orbit
userInput.type           = 'northern';  % 'northern'/'southern' for 'halo' else 'none'
userInput.tolerance      = 1e-6;    % solution tolerance
userInput.jacobianConst  = 3.1; % Required jacobian Constant

tic

globalVar                  = GlobalData(userInput);
fprintf('System /mu value %f\n',globalVar.userInput.mu);

%% Input Orbit Code here

% OrbFamPar = familyPar;
% 
% figure()
% viscircles([0,0],1), grid on; hold on;
% for i = 1:userInput.NoOfFam
%     scatter(real(OrbFamPar.Eigens(i).C_Val(1)),imag(OrbFamPar.Eigens(i).C_Val(1)),'ob','filled');
%     scatter(real(OrbFamPar.Eigens(i).C_Val(2)),imag(OrbFamPar.Eigens(i).C_Val(2)),'ob','filled');
% end
% xlabel('Real')
% ylabel('Imaginary')
% title('Eigen Structure Associated with Center Subspace');
% 
% 
% plotBifurcationSln(globalVar,OrbFamPar);
% toc
% takes 21.493200 seconds for earth moon system(for 30- orbits) with i7 - 4 core and 16g Ram ,without
% parallel pool.
