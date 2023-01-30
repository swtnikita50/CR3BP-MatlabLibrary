
%{
...
Created on  05/6/2022 17:23

This File plots the invariant manifolds along the lyapunov orbit.

Inputs
------
1) G_var - Global Data
2) LyapOrb - Lyapunov orbit parameters we get from LyapOrbitParameters.m
3) n - The number of nodes along the lyapunov orbit to plot manifolds


Outputs
--------


Dependencies
------------
1) stm_X
2) CR3BP
3) varEq
4) Integrator(fun,x0,[0 tspan]);

...
%}

% n = no. of points you need the Initital conditions for

function [Xn] = orbitInvManifoldIC(globalVar,orbPar,n,type, dir)

funVarEq = globalVar.functions.varEq_stmDot;
epsilon = 10^-6;

switch type 
    case 'stable'
        eigVal = orbPar.eigens.val.stable;
        eigVec = orbPar.eigens.vec.stable;
    case 'unstable'
        eigVal = orbPar.eigens.val.unstable;
        eigVec = orbPar.eigens.vec.unstable;
end
tf = orbPar.period;
X0 = orbPar.IC';

k = 1;
Y0_int = eigVec(:,k);


[~,PHItf,~,~,PHI] = stm_X(globalVar,X0,funVarEq,tf);
i = 1;

for m = 1:floor(length(PHI(:,1))/(n-1)):length(PHI(:,1))
    phi = reshape(PHI(m,1:36),6,6);
    X0 = reshape(PHI(m,37:42),6,1);
    Y0 = phi*Y0_int;
    %Normalize Ys0 to 1
    normY0 = norm(Y0(1:3));
    Y0 = Y0/normY0;
    XN(i,:) = X0 + dir*epsilon*Y0;
    i = i+1;
end
Xn(:,:,k) = XN(:,:);
end









