%{
...
Created on  12/07/2022 17:12

This File does the continuation and gets all the Lyapunov orbit parameters.

Inputs
------
1) UserDat - Supplied UserData in main file.
2) G_var   - Global data.
3) e       - specified jacobianContant

Outputs
--------
1) LyapOrb - A structure containing 
                - time      - full time period for lyapunov orbit.
                - IC        - Initial Conditions for Lyapunov Orbit.
                - Energy    - Energy Value for Lyapunov Orbit.
                - Monodromy - Monodromy matrix of Lyapunov Orbit.
                - Eigens    - Eigen Values and Eigen Vectors for Lyapunov Orbit.
                                    - Stable,Unstable and Center Eigen Values
                                    - Stable,Unstable and Center Eigen Vectors

Note that the size of eigen values and eigen vectors might change

Dependencies
------------
1) InitialGuess(PointLoc,G_var)
2) DiffCorrec(X_Guess,Plot,G_var)
                        

Reference for continuation 
--------------------------
1) Wang Sang Koon, Martin W Lo,JE Marsden, Shane D Ross - "Dynamical
   Systems,the Three Body Problem and Space Mission Design", 2011
...
%}
function [LyapOrbFam] = LyapOrbitFamilyParameters(UserDat,G_var)
% Extract the parameters
CorrecPlot = UserDat.CorrectionPlot;
NoOfFam = UserDat.NoOfFam;
EqPoint = UserDat.PointLoc;
mu = G_var.Constants.mu;
funVarEq = G_var.IntFunc.VarEqAndSTMdot;

[XGuessL] = InitialGuess(EqPoint,G_var);
for guess = 1:2
    fprintf('\n===============================================\n')
    fprintf('Obtaining the Corrected Values for guess:- %d\n',guess)
    fprintf('===============================================\n')
    [tCorrec(guess,1),xCorrec(guess,:),~] = DiffCorrec(XGuessL(guess,:),CorrecPlot,G_var);
    Energy(guess,1) = jacobiValue3D(xCorrec(guess,:),mu);
    [~,Monodromy(:,:,guess),~,~] = StateTransAndX(G_var,xCorrec(guess,:),funVarEq,tCorrec(guess,:));
    [Eigens(guess).S_EigVal,Eigens(guess).US_EigVal,Eigens(guess).C_Val,Eigens(guess).S_EigVec,...
        Eigens(guess).US_EigVec,Eigens(guess).C_EigVec] = CalcEigenValVec(Monodromy(:,:,guess),1) ;
    StabilityIdx(guess,1) = CalcStabilityIdx(Eigens(guess));
end

tol = 1e-6;

for family = 3:NoOfFam
    fprintf('\n===============================================\n')
    fprintf('Obtaining the Corrected Values for guess:- %d\n',family)
    fprintf('===============================================\n')
    delta = (xCorrec(family-1,:) - xCorrec(family-2,:));
    GuessX = xCorrec(family-1,:) + delta;
    [~,~,~,isMaxIterReached] = DiffCorrec(GuessX,CorrecPlot,G_var);
    while isMaxIterReached
        delta = delta/2;
        GuessX = xCorrec(family-1,:) + delta;
        [~,~,~,isMaxIterReached] = DiffCorrec(GuessX,CorrecPlot,G_var);
    end
    [tCorrec(family,1),xCorrec(family,:),~,~] = DiffCorrec(GuessX,CorrecPlot,G_var);
    Energy(family,1) = jacobiValue3D(xCorrec(family,:),mu);
    [~,Monodromy(:,:,family),~,~] = StateTransAndX(G_var,xCorrec(family,:),funVarEq,tCorrec(family,:));
    [Eigens(family).S_EigVal,Eigens(family).US_EigVal,Eigens(family).C_Val,Eigens(family).S_EigVec,...
        Eigens(family).US_EigVec,Eigens(family).C_EigVec] = CalcEigenValVec(Monodromy(:,:,family),1) ;
    StabilityIdx(family,1) = CalcStabilityIdx(Eigens(family));

end



LyapOrbFam.time      = tCorrec; %(NoofFam x 1) - Full Orbit Time
LyapOrbFam.IC        = xCorrec; %(NoofFam x UserDat.Dimension)
LyapOrbFam.Energy    = Energy;
LyapOrbFam.Monodromy = Monodromy;
LyapOrbFam.Eigens    = Eigens;
LyapOrbFam.StabilityIdx = StabilityIdx;
end
