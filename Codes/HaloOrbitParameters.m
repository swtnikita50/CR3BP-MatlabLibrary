%{
...
Created on  02/07/22 22:10 - Nikita

This File does the continuation and gets all the Halo orbit parameters.

Inputs
------
1) UserDat - Supplied UserData in main file.
2) G_var   - Global data.
3) e       - specified jacobianContant

Outputs
--------
1) HaloOrb - A structure containing 
                - time      - full time period for Halo orbit.
                - IC        - Initial Conditions for Halo Orbit.
                - Energy    - Energy Value for Halo Orbit.
                - Monodromy - Monodromy matrix of Halo Orbit.
                - Eigens    - Eigen Values and Eigen Vectors for Halo Orbit.
                                    - Stable,Unstable and Center Eigen Values
                                    - Stable,Unstable and Center Eigen Vectors

Note that the size of eigen values and eigen vectors might change

Dependencies
------------
1) InitialGuess(PointLoc,G_var)
2) DiffCorrec(X_Guess,Plot,G_var, isMaxIterReached)
                        

Reference for continuation 
--------------------------
1) Wang Sang Koon, Martin W Lo,JE Marsden, Shane D Ross - "Dynamical
   Systems,the Three Body Problem and Space Mission Design", 2011
...
%}
function [HaloOrb] = HaloOrbitParameters(UserDat,G_var, type,e)
% Extract the parameters
CorrecPlot = UserDat.CorrectionPlot;
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

for guess = 1
    XGuessL(guess,:) = InitialGuess(EqPoint,G_var, 'halo', m);
    fprintf('\n===============================================\n')
    fprintf('Obtaining the Corrected Values for guess:- %d\n',guess)
    fprintf('===============================================\n')
    [tCorrec(guess,1),xCorrec(guess,:),~] = DiffCorrec(XGuessL(guess,:),CorrecPlot,G_var,'halo');
    Energy(guess,1) = jacobiValue3D(xCorrec(guess,:),mu);

    [~,Monodromy(:,:,guess),~,~] = StateTransAndX(G_var,xCorrec(guess,:),funVarEq,tCorrec(guess,:));
    [Eigens(guess).S_EigVal,Eigens(guess).US_EigVal,Eigens(guess).C_Val,Eigens(guess).S_EigVec,...
        Eigens(guess).US_EigVec,Eigens(guess).C_EigVec] = CalcEigenValVec(Monodromy(:,:,guess),1) ;


end

eMax = Energy(1,1);
isOrbitFound = 1;

if e > eMax
    fprintf(']\n===============================================\n')
    fprintf(['Error: The jacobian(energy) value enetered is more \n' ...
        'than the maximum limit for the Halo Orbit of this Libration Point\n'])
    fprintf('===============================================\n');
    isOrbitFound = 0;
else
end

tol = 1e-6;
family = 2;
delX0 = [0 0 1 0 0 0];
S = 0.01*dir;

while Energy(family-1) > e
    fprintf('\n Obtaining the Corrected Values for guess:- %d\n',family);
%     if e < eMin
%         fprintf('\nError: Energy input is less than the %d.\n',eMin);
%         isOrbitFound = 0;
%         break;
%     end
    
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
    if abs(Energy(family)-Energy(family-1)) < tol
        disp('\n Error: The energy difference is negligible. \n');
        isOrbitFound = 0;
        break;
    end
    family = family+1;
    S = (xCorrec(family-1,3) - xCorrec(family-2,3));
end
xCorrecLowerAmp(1,:) = xCorrec(end-1,:); tCorrecLowerAmp(1,1) = tCorrec(end-1,1); EnergyLowerAmp(1,1) = Energy(end-1,1);
xCorrecUpperAmp(1,:) = xCorrec(end,:); tCorrecUpperAmp(1,1) = tCorrec(end,1); EnergyUpperAmp(1,1) = Energy(end,1);


switch isOrbitFound
    case 1
        fprintf('\n===============================================\n')
        fprintf('Found IC for orbit with energy lesser than and greater than %d\n',e)
        fprintf('===============================================\n')


        tol = 1e-10;
        eNew = EnergyUpperAmp(1,1);
        while abs(eNew-e)>tol
            %S = (xCorrecUpperAmp(1,3) - xCorrecLowerAmp(1,3))/2;
            %delta = delX0*S;
            delta = (xCorrecUpperAmp(1,:) - xCorrecLowerAmp(1,:))/2;
            GuessX = xCorrecLowerAmp(1,:) + delta;

            [tNew,xNew,~] = DiffCorrec(GuessX,CorrecPlot,G_var);
            eNew = jacobiValue3D(xNew,mu);
            % CHECK eNew with e and update in xcorrec and everywhere
            if eNew>e
                xCorrecLowerAmp(1,:) = xNew;
                tCorrecLowerAmp(1,1) = tNew;
                EnergyLowerAmp(1,1) = eNew;
            else
                xCorrecUpperAmp(1,:) = xNew;
                tCorrecUpperAmp(1,1) = tNew;
                EnergyUpperAmp(1,1) = eNew;
            end
        end

        [~,Monodromy,~,~] = StateTransAndX(G_var,xNew,funVarEq,tNew);
        [Eigens.S_EigVal,Eigens.US_EigVal,Eigens.C_Val,Eigens.S_EigVec,...
            Eigens.US_EigVec,Eigens.C_EigVec] = CalcEigenValVec(Monodromy,1);
        StabilityIdx = CalcStabilityIdx(Eigens);
        
        HaloOrb.time      = tNew; %(NoofFam x 1) - Full Orbit Time
        HaloOrb.IC        = xNew; %(NoofFam x UserDat.Dimension)
        HaloOrb.Energy    = eNew;
        HaloOrb.Monodromy = Monodromy;
        HaloOrb.Eigens    = Eigens;
        HaloOrb.StabilityIdx = StabilityIdx;
    case 0

        if e > eMax
            fprintf('\n===============================================\n')
            fprintf('Error: The jacobian(energy) value enetered is greater than the maximum limit for the Libration Point\n')
            fprintf('===============================================\n');
        else
            fprintf('\nError: Energy input is less than the %d.\n',Energy(end));
        end
        clear HaloOrb;
        HaloOrb = -1;
end

end