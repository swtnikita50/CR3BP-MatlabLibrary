%{ 
...
Modified and Updated on 7 Aug 2019 -Karthi
Modified on 30/6/2022 to generate a full cr3bp library- Nikita
(last updated : (1) 21:26, 10/8/2019  (2) 16:06 13/08/2019)  
(3) Including Manifolds 14/08/2019 - 15/08/2019 (4) Eigen Value Plot and
others 20/08/2019
4) Halo Orbit: 03/07/22

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
  2) LyapOrbitParameters     - will call I(8) and I(1)
  3) Plotter                 - will call I(10),I(11),I(12)
 
 

File dependencies (Implicit - I)
  1) DiffCorrec             - will call E(1) I(3) I(4) I(5)
  2) Integrator             - will call I(5) 
  3) StateTransAndX         - willcall I(6),I(3) and I(3) will call (7)
  4) CRes3BP_EOM
  5) VarEq_Init 
  6) VarEqAndSTMDOT         - will call E(1) I(8) I(5) 
  7) Dfmatrix3D
  8) InitialGuess 
  9) SunPlanetMoonParameters
  10) PlotContourEquilPoints  - will call E(1)
  11) LyapunovPlotter         - Plots the Lyapunov Orbits
  12) plotEqPointManifold     - plots eq.points manifolds

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
6) ENTER ONLY "UserDat" asked for . This is fully sufficient.
7) To know the flow see the file depencies above(if required).


Description for entering Primary Secondary
-------------------------------------------
1) Primary can be 'Sun' or any of the nine
    Planets('Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus',
    'Neptune' and 'Pluto')
2) Secondary any of nine planets(if primary is sun) or natural satellites for the planets
3) Case does not matter while entering any of these two 
4) Enter as string inside double or single quotes.

Natural Satellites included are 
Jupiter - 'Callisto','Ganymede','Europa' and 'Io'.
Saturn - 'Titan'
Neptune - 'Triton'.
Earth - 'moon'

...
%}

%% --------------------------------- Mainline Program ----------------------------------------------------------------
clearvars;clc;close all

%% Get the UserData
% ------------------
% The program can handle both 2 and 3 dimension
UserDat.Dimension      = 3;% input('Give the problem dimension: ');
% Input the Primary/Secondary
% See the description above before entering Primary and Secondary
UserDat.mu        = 0.0121505856;%input('Enter the primary body (string input): ');
% Input for which equilibrium points you want the data for?
UserDat.PointLoc       = 1;%input('For how many  equilibrium points do you want orbit parameters? : ');
UserDat.NoOfFam        = 20;%input('Number of Family? : ');
% Check if you want to see correction plot(Enter 0 or 1).
UserDat.CorrectionPlot = 1;%input('Do You want to see correction plot of differential correction?:');


tic

G_var                  = GlobalData(UserDat);
fprintf('mu value %f\n',G_var.Constants.mu);

[LyapOrbFam] = LyapOrbitFamilyParameters(UserDat,G_var);
figure()
%plot(LyapOrbFam.Energy, LyapOrbFam.StabilityIdx(:,1).Saddle);hold on;
plot(LyapOrbFam.Energy, LyapOrbFam.time,'-r','LineWidth',2);
% Ax = 0.1025;
% % [HaloOrbFam]              = HaloOrbitFamilyParameters(UserDat,G_var, 'southern',3);
% % ans1 = PeriodicOrbitInvariantMfdsIC(G_var,HaloOrb,80, 'stable',1);
% % figure()
% % set(0,'DefaultAxesColorOrder',jet(length(HaloOrb.Energy)));
% % for i = 1:length(HaloOrb.Energy)
% %     [t,x] = Integrator(G_var,G_var.IntFunc.EOM,HaloOrb.IC(i,:),[0 HaloOrb.time(i)],'forward');
% %     plot3(x(:,1),x(:,2),x(:,3),'LineWidth',2);hold on; grid on;
% %     scatter3(1-UserDat.mu,0,0,'p','filled');
% % end
% % for i = 1:length(ans1(:,1))
% %     [t,x] = Integrator(G_var,G_var.IntFunc.EOM,ans1(i,1:6),[0 3*HaloOrb.time],'backward');
% %     plot3(x(:,1),x(:,2),x(:,3),'g');hold on; grid on;
% %     plot3(1-UserDat.mu,0,0,'p');
% %end


% [LyapOrb]              = LyapOrbitParameters(UserDat,G_var, 2.8);
% ans = PeriodicOrbitInvariantMfdsIC(G_var,LyapOrb,80, 'stable',1);
% figure()
% for i = 1:length(ans(:,1))
%     [t,x] = Integrator(G_var,G_var.IntFunc.EOM,ans(i,1:6),[0 2*LyapOrb.time],'backward');
%     plot(x(:,1),x(:,2),'g');hold on; grid on;
%     plot(1-UserDat.mu,0,'p');
% end

toc
% takes 21.493200 seconds for earth moon system(for 30- orbits) with i7 - 4 core and 16g Ram ,without
% parallel pool.






























%% Not Required (This is natural Parameter Continuation)
% ================================================================================================
% [XGuess] = InitialGuess(UserDat.PointLoc);
% [Corrected.time,Corrected.InitialX,Corrected.DF] =...
% GrebowContinuation(XGuess(1),UserDat.NoofFam,UserDat.CorrectionPlot);
% =================================================================================================