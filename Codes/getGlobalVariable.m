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
1) globalVar - structure with fields 
                  - const - structure with fields
                        - mu
                        - gam
                        - Ax1
                        - Ax2
                  - lagPts - structure with fields
                        - Gamma (distance fom the Eq.pts to primaries)
                        - L1,L2,L3,L4,L5(location of eq.points)
                        - Energy(a structure with fields L1 to L5)
                        - EigeVec(a structure with fields L1 to L5)
                        - EigeVal(a structure with fields L1 to L5)
                  - internalFunc - structure with fields
                        - odeOptions
                        - dynamics
                        - varEq_STMdot

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
function globalVar = getGlobalVariable(userInput)

globalVar.userInput             = userInput;
globalVar.const.gam             = (userInput.mu/3)^(1/3) ;
globalVar.const.Ax1             = 2e-2*globalVar.const.gam;% initial amplitude 1 for seed orbit
globalVar.const.Ax2             = 2*globalVar.const.Ax1; % initial amplitude 2 for second seed orbit

globalVar.lagPts                    = equilPts(userInput.mu); % see "equil_pts_position.m"

globalVar.const.jacobianMax       = globalVar.lagPts.jacobianConst(userInput.lagrangePt);

globalVar.functions.odeOptions      = odeset('Reltol',3.e-14,'Abstol',1.e-16);% added 24/2/2020 16:10
globalVar.functions.systemDynamics  = @(t,x) CR3BP(t,x,userInput.mu);% added 24/2/2020 18:55
globalVar.functions.varEq_stmDot    = @(t,x) varEq_stmDot(t,x,userInput.mu);% added 24/2/2020 18:55


                                                                                   

