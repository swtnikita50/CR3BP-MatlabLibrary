%{
...
Created on Feb 20 2020 by 16:05 

This file calculates
 
Inputs
------
1) mu - Mass Parameter

Outputs
-------
1) lagPts a structure file conating the following
        * LagrangePoints
        * gamma Values for L1,L2 and L3
        * jacobianConst values at each Lagrange point
        * Eigen Values and Eigen Vectors of each lagrange points

Since all these are Constants I have moved this file to "GlobalData.m"

=========================================================================
Source File  : Prof.Shane Ross (revised 2.19.04) ,till "xPos,yPos"
WebSite      : http://www.dept.aoe.vt.edu/~sdross/
Modification : calculate jacobianConst , eigen values,eigen vectors and store in structure 
=========================================================================

...
%}
function[lagPts]= equilPts(mu)
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
lagPts.gamma = [gamma1;gamma2;gamma3];

xL1 = 1-mu-gamma1;
xL2 = 1-mu+gamma2;
xL3 = -mu-gamma3;
xL4 = 0.5-mu;
xL5 = xL4;
yL4 = sqrt(3)/2;
yL5 = -yL4;
yL1 = 0;
yL2 = 0;
yL3 = 0;

% Lagrange Points 
lagPts.pos(1,:) = [xL1,yL1];
lagPts.pos(2,:) = [xL2,yL2];
lagPts.pos(3,:) = [xL3,yL3];
lagPts.pos(4,:) = [xL4,yL4];
lagPts.pos(5,:) = [xL5,yL5];


% Calculate The Energies at each Lagrange Point


lagPts.jacobianConst(1) = jacobiValue3D([lagPts.pos(1,1),lagPts.pos(1,2),0,0,0,0],mu);
lagPts.jacobianConst(2) = jacobiValue3D([lagPts.pos(2,1),lagPts.pos(2,2),0,0,0,0],mu);
lagPts.jacobianConst(3) = jacobiValue3D([lagPts.pos(3,1),lagPts.pos(3,2),0,0,0,0],mu);
lagPts.jacobianConst(4) = jacobiValue3D([lagPts.pos(4,1),lagPts.pos(4,2),0,0,0,0],mu);
lagPts.jacobianConst(5) = jacobiValue3D([lagPts.pos(5,1),lagPts.pos(5,2),0,0,0,0],mu);

% Find the Eigen Values and Eigen Vectors at each Lagrange Point
 [lagPts.eig.Val.stable.L1,lagPts.eig.Val.unstable.stable.L1,lagPts.eig.Val.center.L1,lagPts.eig.Vec.stable.L1,lagPts.eig.Vec.unstable.L1,lagPts.eig.Vec.center.L1,~,~] =...
     calcEigen(Dfmatrix3D([lagPts.pos(1,1),lagPts.pos(1,2),0,0,0,0],mu),0);
 [lagPts.eig.Val.stable.L2,lagPts.eig.Val.unstable.stable.L2,lagPts.eig.Val.center.L2,lagPts.eig.Vec.stable.L2,lagPts.eig.Vec.unstable.L2,lagPts.eig.Vec.center.L2,~,~] =...
     calcEigen(Dfmatrix3D([lagPts.pos(2,1),lagPts.pos(2,2),0,0,0,0],mu),0);
  [lagPts.eig.Val.stable.L3,lagPts.eig.Val.unstable.stable.L3,lagPts.eig.Val.center.L3,lagPts.eig.Vec.stable.L3,lagPts.eig.Vec.unstable.L3,lagPts.eig.Vec.center.L3,~,~] =...
     calcEigen(Dfmatrix3D([lagPts.pos(3,1),lagPts.pos(3,2),0,0,0,0],mu),0);
  [lagPts.eig.Val.stable.L4,lagPts.eig.Val.unstable.stable.L4,lagPts.eig.Val.center.L4,lagPts.eig.Vec.stable.L4,lagPts.eig.Vec.unstable.L4,lagPts.eig.Vec.center.L4,~,~] =...
     calcEigen(Dfmatrix3D([lagPts.pos(4,1),lagPts.pos(4,2),0,0,0,0],mu),0);
  [lagPts.eig.Val.stable.L5,lagPts.eig.Val.unstable.stable.L5,lagPts.eig.Val.center.L5,lagPts.eig.Vec.stable.L5,lagPts.eig.Vec.unstable.L5,lagPts.eig.Vec.center.L5,~,~] =...
     calcEigen(Dfmatrix3D([lagPts.pos(1,1),lagPts.pos(1,2),0,0,0,0],mu),0);

% Lx         - x can be 1:5
% eig.Val.stable.Lx  - Stable Eigen Value
% eig.Val.unstable.stable.Lx - UnStable Eigen Value
% eig.Val.center.Lx  - Center Eigen Value
% eig.Vec.stable.Lx  - Stable Eigen Vector
% eig.Vec.unstable.Lx - UnStable Eigen Vector
% eig.Vec.center.Lx  - Center Eigen Vector

end