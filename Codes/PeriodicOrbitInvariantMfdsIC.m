
%{
...
Created on  05/6/2022 17:23

This File plots the invariant manifolds along the lyapunov orbit.

Inputs
------
1) G_var - Global Data
2) LyapOrb - Lyapunov orbit parameters we get from LyapOrbitParameters.m
3) n - The number of nodes along the lyapunov orbit to plot manifolds
4) system - defines external or internal system
5) libNo - Number of liberation point used
6) orbitNo - number of orbit used from the family

Outputs
--------


Dependencies
------------
1) StateTransAndX
2) EOM
3) VarEq
4) Integrator(fun,x0,[0 tspan]);

...
%}

function [Xn] = PeriodicOrbitInvariantMfdsIC(G_var,PeriodicOrbPar,n,type, dir)


funVarEq = G_var.IntFunc.VarEqAndSTMdot;
epsilon = 10^-6;

switch type 
    case 'stable'
        eigVec = PeriodicOrbPar.Eigens.S_EigVec;
        eigVal = PeriodicOrbPar.Eigens.S_EigVal;
    case 'unstable'
        eigVec = PeriodicOrbPar.Eigens.US_EigVec;
        eigVal = PeriodicOrbPar.Eigens.US_EigVal;
end
tf = PeriodicOrbPar.time;
X0 = PeriodicOrbPar.IC';


k = 1;
        Y0int = eigVec(:,k);


    [~,PHItf,~,~,PHI] = StateTransAndX(G_var,X0,funVarEq,tf);
    i = 1;
    for m = 1:floor(length(PHI(:,1))/(n-1)):length(PHI(:,1))
        phi = reshape(PHI(m,1:36),6,6);
        X0 = reshape(PHI(m,37:42),6,1);
        Y0 = phi*Y0int;
        %Normalize Ys0 to 1
        normY0 = norm(Y0(1:3));
        Y0 = Y0/normY0;
        XN(i,:) = X0 + dir*epsilon*Y0;
        i = i+1;
    end
    Xn(:,:,k) = XN(:,:);
end


                






