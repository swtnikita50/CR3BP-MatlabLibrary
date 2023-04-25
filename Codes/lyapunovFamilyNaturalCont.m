function [familyPar] = lyapunovFamilyNaturalCont(globalVar)
%% Extract the parameters
orbitCount = globalVar.userInput.orbitCount;
mu = globalVar.userInput.mu;
funVarEq = globalVar.functions.varEq_stmDot;
tol = globalVar.userInput.tolerance;

%% Obtain Initial Guess
[guess] = initGuess(globalVar);
delX0 = [-1 0 0 0 0 0];
dels = 0.0001;

%% Differential Correction and Natural Parameter Continuation
fprintf('\n===============================================\n')
for i = 1:orbitCount
    fprintf('Starting differential correction for orbit no.: %d\n',i)
    
    if i > 2
        delta = dels*delX0;
        xGuess = x(i-1,:) + delta;
        [~,~,~,isMaxIterReached] = diffCorrec(xGuess,globalVar);
        while isMaxIterReached
            delta = 4*delta/5;
            xGuess = x(i-1,:) + delta;
            [~,~,~,isMaxIterReached] = diffCorrec(xGuess,globalVar);
        end
    else
        xGuess = guess(i,:);
    end
    [t(i,1),x(i,:),~] = diffCorrec(xGuess,globalVar);
    jacobianConst(i,1) = jacobiValue3D(x(i,:),mu);
    [~,monodromy(:,:,i),~,~] = stm_X(globalVar,x(i,:),funVarEq,t(i,:));
    [eigens(i).val.stable,eigens(i).val.unstable,eigens(i).val.center,eigens(i).val.p,eigens(i).vec.stable,...
        eigens(i).vec.unstable,eigens(i).vec.center,eigens(i).vec.p] = calcEigen(monodromy(:,:,i),1) ;
    stabilityIdx(i,1) = calcStabilityIdx(eigens(i));

end
fprintf('\n===============================================\n')

%% Output Data
familyPar.period      = t; %(orbitCount x 1) - Full Orbit Time
familyPar.IC        = x; %(orbitCount x userInput.Dimension)
familyPar.jacobianConst    = jacobianConst;
familyPar.monodromy = monodromy;
familyPar.eigens    = eigens;
familyPar.stabilityIdx = stabilityIdx;
end