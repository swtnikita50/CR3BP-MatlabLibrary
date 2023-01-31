%{
...
Created on  12/07/2022 17:12

This File does the continuation and gets all the Lyapunov orbit parameters.

Inputs
------
1) userInput - Supplied userInputa in main file.
2) globalVar   - Global data.
3) e       - specified jacobianContant

Outputs
--------
1) LyapOrb - A structure containing 
                - time      - full time period for lyapunov orbit.
                - IC        - Initial Conditions for Lyapunov Orbit.
                - jacobianConst    - jacobianConst Value for Lyapunov Orbit.
                - Monodromy - Monodromy matrix of Lyapunov Orbit.
                - eigens    - Eigen Values and Eigen Vectors for Lyapunov Orbit.
                                    - Stable,Unstable and Center Eigen Values
                                    - Stable,Unstable and Center Eigen Vectors

Note that the size of eigen values and eigen vectors might change

Dependencies
------------
1) InitialGuess(PointLoc,globalVar)
2) diffCorrec(X_Guess,Plot,globalVar)
                        

Reference for continuation 
--------------------------
1) Wang Sang Koon, Martin W Lo,JE Marsden, Shane D Ross - "Dynamical
   Systems,the Three Body Problem and Space Mission Design", 2011
...
%}
function [familyPar] = lyapunovFamily(globalVar)
%% Extract the parameters
orbitCount = globalVar.userInput.orbitCount;
mu = globalVar.userInput.mu;
funVarEq = globalVar.functions.varEq_stmDot;
tol = globalVar.userInput.tolerance;

%% Obtain Initial Guess
[guess] = initGuess(globalVar);

%% Differential Correction and Natural Parameter Continuation
fprintf('\n===============================================\n')
for i = 1:orbitCount
    fprintf('Starting differential correction for orbit no.: %d\n',i)
    
    if i > 2
        delta = (x(i-1,:) - x(i-2,:));  % Continuation Step
        xGuess = x(i-1,:) + delta;
        [~,~,~,isMaxIterReached] = diffCorrec(xGuess,globalVar);
        while isMaxIterReached
            delta = delta/2;
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
