function [familyPar] = haloFamilyPseudoArcLengthCont(globalVar)
%% Extract the parameters
orbitCount = globalVar.userInput.orbitCount;
mu = globalVar.userInput.mu;
funVarEq = globalVar.functions.varEq_stmDot;
type = globalVar.userInput.type;
tol = globalVar.userInput.tolerance;

%% Initial Guess and Orbit Solution
switch type
    case 'northern'
        m = 1;
        dir = 1;
    case 'southern'
        m = 3;
        dir = -1;
end
fprintf('\n===============================================\n')
for i = 1
    xGuess(i,:) = initGuess(globalVar);
    fprintf('Starting differential correction for orbit no.: %d\n',i)
    
    [t(i,1),x(i,:),~] = diffCorrec(xGuess(i,:),globalVar);
    jacobianConst(i,1) = jacobiValue3D(x(i,:),mu);
    [~,monodromy(:,:,i),~,~] = stm_X(globalVar,x(i,:),funVarEq,t(i,:));
    [eigens(i).val.stable,eigens(i).val.unstable,eigens(i).val.center,eigens(i).val.p,eigens(i).vec.stable,...
        eigens(i).vec.unstable,eigens(i).vec.center,eigens(i).vec.p] = calcEigen(monodromy(:,:,i),1) ;
    stabilityIdx(i,1) = calcStabilityIdx(eigens(i));

end


delX0 = [0 0 1 0 0 0];
S = 0.01*dir;
for i = 2
    fprintf('Starting differential correction for orbit no.: %d\n',i)
    
    delta = delX0*S;        % Continutation Step
    xGuess = x(i-1,:) + delta;
    [~,~,~,isMaxIterReached] = diffCorrec(xGuess,globalVar);
    while isMaxIterReached
        delta = delta/2;
        xGuess = x(i-1,:) + delta;
        [~,~,~,isMaxIterReached] = diffCorrec(xGuess,globalVar);
    end
    [t(i,1),x(i,:),~,~] = diffCorrec(xGuess,globalVar);
    jacobianConst(i,1) = jacobiValue3D(x(i,:),mu);
    [~,monodromy(:,:,i),~,~] = stm_X(globalVar,x(i,:),funVarEq,t(i,:));
    [eigens(i).val.stable,eigens(i).val.unstable,eigens(i).val.center,eigens(i).val.p,eigens(i).vec.stable,...
        eigens(i).vec.unstable,eigens(i).vec.center,eigens(i).vec.p] = calcEigen(monodromy(:,:,i),1) ;
    stabilityIdx(i,1) = calcStabilityIdx(eigens(i));
end

for i = 3:orbitCount
    fprintf('Starting differential correction for orbit no.: %d\n',i)
    
    delta = (x(i-1,:) - x(i-2,:));  % Continuation Step
    xGuess = x(i-1,:) + delta;
    [~,~,~,isMaxIterReached] = diffCorrec(xGuess,globalVar);
    while isMaxIterReached
        delta = delta/2;
        xGuess = x(i-1,:) + delta;
        [~,~,~,isMaxIterReached] = diffCorrec(xGuess,globalVar);
    end
    
    [t(i,1),x(i,:),~] = diffCorrec(xGuess,globalVar);
    jacobianConst(i,1) = jacobiValue3D(x(i,:),mu);
    [~,monodromy(:,:,i),~,~] = stm_X(globalVar,x(i,:),funVarEq,t(i,:));
    [eigens(i).val.stable,eigens(i).val.unstable,eigens(i).val.center,eigens(i).val.p,eigens(i).vec.stable,...
        eigens(i).vec.unstable,eigens(i).vec.center,eigens(i).vec.p] = calcEigen(monodromy(:,:,i),1) ;
    stabilityIdx(i,1) = calcStabilityIdx(eigens(i));

end

%% Output Data
familyPar.period      = t; %(NoofFam x 1) - Full Orbit Time
familyPar.IC        = x; %(NoofFam x UserDat.Dimension)
familyPar.jacobianConst    = jacobianConst;
familyPar.monodromy = monodromy;
familyPar.eigens    = eigens;
familyPar.stabilityIdx = stabilityIdx;
    
end