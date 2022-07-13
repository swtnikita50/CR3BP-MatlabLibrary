%{
...
Created on Feb 20 2020 by 16:05 

This file calculates
 
Inputs
------
1) mu - Mass Parameter

Outputs
-------
1) LagPts a structure file conating the following
        * LagrangePoints
        * Gamma Values for L1,L2 and L3
        * Energy values at each Lagrange point
        * Eigen Values and Eigen Vectors of each lagrange points

Since all these are Constants I have moved this file to "GlobalData.m"

=========================================================================
Source File  : Prof.Shane Ross (revised 2.19.04) ,till "xPos,yPos"
WebSite      : http://www.dept.aoe.vt.edu/~sdross/
Modification : calculate energy , eigen values,eigen vectors and store in structure 
=========================================================================

...
%}
function[LagPts]= equil_pts_position(mu)
% calculates the gamma values and possitions of equilibrium points


mu2=1-mu;


poly1 = [1   -1*(3-mu)  (3-2*mu)  -mu   2*mu  -mu ];
poly2 = [1      (3-mu)  (3-2*mu)  -mu  -2*mu  -mu ];
poly3 = [1      (2+mu)  (1+2*mu)  -mu2 -2*mu2 -mu2];

rt1 = roots(poly1);
rt2 = roots(poly2);
rt3 = roots(poly3);

for k=1:5
    if isreal(rt1(k)) 
        gamma1=rt1(k); 
    end
    if isreal(rt2(k)) 
        gamma2=rt2(k); 
    end
    if isreal(rt3(k)) 
        gamma3=rt3(k); 
    end
   
end
LagPts.Gamma = [gamma1;gamma2;gamma3];

xL1=1-mu-gamma1;
xL2=1-mu+gamma2;
xL3=-mu-gamma3;
xL4=0.5-mu;
xL5=xL4;
yL4=sqrt(3)/2;
yL5=-yL4;
yL1=0;

xPos=[xL1,xL2,xL3,xL4,xL5];
yPos=[yL1,yL1,yL1,yL4,yL5];

% Lagrange Points 
LagPts.L1 = [xPos(1),yPos(1)];
LagPts.L2 = [xPos(2),yPos(2)];
LagPts.L3 = [xPos(3),yPos(3)];
LagPts.L4 = [xPos(4),yPos(4)];
LagPts.L5 = [xPos(5),yPos(5)];


% Calculate The Energies at each Lagrange Point


LagPts.Energy.L1 = jacobiValue3D([LagPts.L1,0,0,0,0],mu);
LagPts.Energy.L2 = jacobiValue3D([LagPts.L2,0,0,0,0],mu);
LagPts.Energy.L3 = jacobiValue3D([LagPts.L3,0,0,0,0],mu);
LagPts.Energy.L4 = jacobiValue3D([LagPts.L4,0,0,0,0],mu);
LagPts.Energy.L5 = jacobiValue3D([LagPts.L5,0,0,0,0],mu);

% Find the Eigen Values and Eigen Vectors at each Lagrange Point
 [LagPts.SEigVal.L1,LagPts.USEigVal.L1,LagPts.CEigVal.L1,LagPts.SEigVec.L1,LagPts.USEigVec.L1,LagPts.CEigVec.L1,~,~] =...
     CalcEigenValVec(Dfmatrix3D([LagPts.L1,0,0,0,0],mu),0);
 [LagPts.SEigVal.L2,LagPts.USEigVal.L2,LagPts.CEigVal.L2,LagPts.SEigVec.L2,LagPts.USEigVec.L2,LagPts.CEigVec.L2,~,~] =...
     CalcEigenValVec(Dfmatrix3D([LagPts.L2,0,0,0,0],mu),0);
  [LagPts.SEigVal.L3,LagPts.USEigVal.L3,LagPts.CEigVal.L3,LagPts.SEigVec.L3,LagPts.USEigVec.L3,LagPts.CEigVec.L3,~,~] =...
     CalcEigenValVec(Dfmatrix3D([LagPts.L3,0,0,0,0],mu),0);
  [LagPts.SEigVal.L4,LagPts.USEigVal.L4,LagPts.CEigVal.L4,LagPts.SEigVec.L4,LagPts.USEigVec.L4,LagPts.CEigVec.L4,~,~] =...
     CalcEigenValVec(Dfmatrix3D([LagPts.L4,0,0,0,0],mu),0);
  [LagPts.SEigVal.L5,LagPts.USEigVal.L5,LagPts.CEigVal.L5,LagPts.SEigVec.L5,LagPts.USEigVec.L5,LagPts.CEigVec.L5,~,~] =...
     CalcEigenValVec(Dfmatrix3D([LagPts.L5,0,0,0,0],mu),0);

 % Lx         - x can be 1:5
% SEigVal.Lx  - Stable Eigen Value
% USEigVal.Lx - UnStable Eigen Value
% CEigVal.Lx  - Center Eigen Value
% SEigVec.Lx  - Stable Eigen Vector
% USEigVec.Lx - UnStable Eigen Vector
% CEigVec.Lx  - Center Eigen Vector

end