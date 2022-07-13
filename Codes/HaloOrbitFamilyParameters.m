% Created on 12-07-22 (17:02)


function [HaloOrbFam] = HaloOrbitFamilyParameters(UserDat,G_var, type)
% Extract the parameters
CorrecPlot = UserDat.CorrectionPlot;
NoOfFam = UserDat.NoOfFam;
EqPoint = UserDat.PointLoc;
mu = G_var.Constants.mu;
funVarEq = G_var.IntFunc.VarEqAndSTMdot;

switch type
    case 'northern'
        m = 1;
        dir = 1;
    case 'southern'
        m = 3;
        dir = -1;
end

for guess = 1:2
    XGuessL(guess,:) = InitialGuess(EqPoint,G_var, 'halo', m);
    fprintf('\n===============================================\n')
    fprintf('Obtaining the Corrected Values for guess:- %d\n',guess)
    fprintf('===============================================\n')
    [tCorrec(guess,1),xCorrec(guess,:),~] = DiffCorrec(XGuessL(guess,:),CorrecPlot,G_var,'halo');
    Energy(guess,1) = jacobiValue3D(xCorrec(guess,:),mu);

    [~,Monodromy(:,:,guess),~,~] = StateTransAndX(G_var,xCorrec(guess,:),funVarEq,tCorrec(guess,:));
    [Eigens(guess).S_EigVal,Eigens(guess).US_EigVal,Eigens(guess).C_Val,Eigens(guess).S_EigVec,...
        Eigens(guess).US_EigVec,Eigens(guess).C_EigVec] = CalcEigenValVec(Monodromy(:,:,guess),1) ;
    StabilityIdx(guess,1) = CalcStabilityIdx(Eigens(guess));

end

tol = 1e-6;
delX0 = [0 0 1 0 0 0];
S = 0.01*dir;

for family = 3:NoOfFam+1
    fprintf('\n===============================================\n')
    fprintf('Obtaining the Corrected Values for guess:- %d\n',family)
    fprintf('===============================================\n')
    delta = delX0*S;
    GuessX = xCorrec(family-1,:) + delta;
    [~,~,~,isMaxIterReached] = DiffCorrec(GuessX,CorrecPlot,G_var,'halo');
    while isMaxIterReached
        delta = delta/2;
        GuessX = xCorrec(family-1,:) + delta;
        [~,~,~,isMaxIterReached] = DiffCorrec(GuessX,CorrecPlot,G_var,'halo');
    end
    [tCorrec(family,1),xCorrec(family,:),~,~] = DiffCorrec(GuessX,CorrecPlot,G_var,'halo');
    Energy(family,1) = jacobiValue3D(xCorrec(family,:),mu);



        [~,Monodromy(:,:,family),~,~] = StateTransAndX(G_var,xCorrec(family,:),funVarEq,tCorrec(family));
        [Eigens(family).S_EigVal,Eigens(family).US_EigVal,Eigens(family).C_Val,Eigens(family).S_EigVec,...
        Eigens(family).US_EigVec,Eigens(family).C_EigVec] = CalcEigenValVec(Monodromy(:,:,family),1) ;
    StabilityIdx(family,1) = CalcStabilityIdx(Eigens(family));
    S = (xCorrec(family-1,3) - xCorrec(family-2,3));
end

        
        HaloOrbFam.time      = tCorrec; %(NoofFam x 1) - Full Orbit Time
        HaloOrbFam.IC        = xCorrec; %(NoofFam x UserDat.Dimension)
        HaloOrbFam.Energy    = Energy;
        HaloOrbFam.Monodromy = Monodromy;
        HaloOrbFam.Eigens    = Eigens;
        HaloOrbFam.StabilityIdx = StabilityIdx;
    

end