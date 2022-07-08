%{
...
Created on 22/2/2020 18:40
This files store some constant values which might be required else where
and some frequently used functions.

This file calculates below mentioned values 

1) mu - Mass paprameter
2) Equilibrium Points Location (x,y values)
3) Energy at each equilibrium points 
4) Eigen Values and Eigen Vectors of each equlibrium point.
5) Ax1  and Ax2 - initial seed amplitudes to calulate two seed initial
   guesses.
6) Energy to plot the contour plots

stores the below mentioned functions

7) Options for ODE
8) Equation of motion function
9) Variational equation function

Inputs
-------
1) Primary and Secondary bodies (user defined in "MAIN_LyapOrbit")

Outputs
--------
1) G_var - structure with fields 
                  - Constants - structure with fields
                        - mu
                        - gam
                        - Ax1
                        - Ax2
                  - LagPts - structure with fields
                        - Gamma (distance fom the Eq.pts to primaries)
                        - L1,L2,L3,L4,L5(location of eq.points)
                        - Energy(a structure with fields L1 to L5)
                        - EigeVec(a structure with fields L1 to L5)
                        - EigeVal(a structure with fields L1 to L5)
                  - IntFunc - structure with fields
                        - ODEoptions
                        - EOM
                        - VarEqAndSTMdot

Dependencies
------------
mu -  calculated through "SunPlanetMoonParameters"(see description of that
      inside the file)
Equilibrium points and related values(2,3 and 4) are calculated through
"equil_pts_position"( see description)

Reference
----------
1) Wang Sang Koon, Martin W Lo,JE Marsden, Shane D Ross - "Dynamical
   Systems,the Three Body Problem and Space Mission Design", 2011

2) Code form Prof.Shane Ross

Initial seed amplitude value calculations are fron Ref(2)
...
%}
function G_var = GlobalData(UserDat)

G_var.Constants.mu              = UserDat.mu;
G_var.Constants.gam             = (G_var.Constants.mu/3)^(1/3) ;
G_var.Constants.Ax1             = 2e-2*G_var.Constants.gam;% initial amplitude 1 for seed orbit
G_var.Constants.Ax2             = 2*G_var.Constants.Ax1; % initial amplitude 2 for second seed orbit

G_var.LagPts                    = equil_pts_position(G_var.Constants.mu); % see "equil_pts_position.m"

G_var.Constants.ReqEnergy       = G_var.LagPts.Energy.L2;

G_var.IntFunc.ODEoptions        = odeset('Reltol',3.e-14,'Abstol',1.e-16);% added 24/2/2020 16:10
G_var.IntFunc.EOM               = @(t,x) CRes3BP_EOM(t,x,G_var.Constants.mu);% added 24/2/2020 18:55
G_var.IntFunc.VarEqAndSTMdot    = @(t,x) VarEqAndSTMDOT(t,x,G_var.Constants.mu);% added 24/2/2020 18:55
G_var.UserDat = UserDat;

                                                                                   

