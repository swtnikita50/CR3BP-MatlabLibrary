%{
...
Created on  30/6/2022 18:47
added jacobianConst sensitivity to the function on 30/6/2022

This File does the continuation and gets all the Lyapunov orbit parameters.

Inputs
------
1) UserDat - Supplied UserData in main file.
2) globalVar   - Global data.
3) jacobianConst      - specified jacobianContant

Outputs
--------
1) orbPar - A structure containing 
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
function [orbPar] = lyapunovOrbit(globalVar)
% Extract the parameters
jacobianConstReq = globalVar.userInput.jacobianConst;
lagrangePt = globalVar.userInput.lagrangePt;
mu = globalVar.userInput.mu;
funVarEq = globalVar.functions.varEq_stmDot;
tol = globalVar.userInput.tolerance;

[xGuess] = initGuess(globalVar);

jacobianMax = globalVar.const.jacobianMax;

isOrbitFound = 1;

if jacobianConstReq> jacobianMax
    fprintf(']\n===============================================\n')
    fprintf('Error: The jacobian(energy) value enetered is more than the maximum limit for the Libration Point\n')
    fprintf('===============================================\n');
    isOrbitFound = 0;
else
end
fprintf('\n===============================================\n')
    
for i = 1:2
    fprintf('Starting differential correction for orbit no.: %d\n',i)
    [t(i,1),x(i,:),~] = diffCorrec(xGuess(i,:),globalVar);
    jacobianConst(i,1) = jacobiValue3D(x(i,:),mu);
end

i = 3;
jacobianConstMin = 2.6;         % Check how to set this value?

if jacobianConstReq> jacobianConst(1)
    xSubs(1,:) = x(1,:); 
    tSubs(1,1) = t(1,1); 
    jacobianConstSubs(1,1) = jacobianConst(1,1);
    
    x(1,:) = x(2,:);
    t(1,1) = t(2,1);
    jacobianConst(1,1) = jacobianConst(2,1);
    
    x(2,:) = xSubs(1,:);
    t(2,1) = tSubs(1,1);
    jacobianConst(2,1) = jacobianConstSubs(1,1);
    while jacobianConst(i-1) < jacobianConstReq
    fprintf('Starting differential correction for orbit no.: %d\n',i)
    
    delta = (x(i-1,:) - x(i-2,:))/2;
    xGuess = x(i-1,:) + delta;
    [~,~,~,isMaxIterReached] = diffCorrec(xGuess,globalVar);
    while isMaxIterReached
        delta = delta/2;
        xGuess = x(i-1,:) + delta;
        [~,~,~,isMaxIterReached] = diffCorrec(xGuess,globalVar);
    end
    [t(i,1),x(i,:),~,~] = diffCorrec(xGuess,globalVar);
    jacobianConst(i,1) = jacobiValue3D(x(i,:),mu);
    if abs(jacobianConst(i)-jacobianConst(i-1)) < tol
        disp('\n Error: The Jacobian Const difference is negligible. \n');
        isOrbitFound = 0;
        break;
    end
    i = i+1;
    end
    % Updating Vaues for the next iteration
    xLowerAmp(1,:) = x(end,:); tLowerAmp(1,1) = t(end,1); jacobianConstLowerAmp(1,1) = jacobianConst(end,1);
    xUpperAmp(1,:) = x(end-1,:); tUpperAmp(1,1) = t(end-1,1); jacobianConstUpperAmp(1,1) = jacobianConst(end-1,1);

elseif jacobianConstReq< jacobianConst(2)
% Check if jacobianConst is between these two energies
%% Start Continuation
while jacobianConst(i-1) > jacobianConstReq
    fprintf('Starting differential correction for orbit no.: %d\n',i)
    if jacobianConstReq < jacobianConstMin
        fprintf('\nError: jacobianConst input is less than the %d.\n',jacobianConstMin);
        isOrbitFound = 0;
        break;
    end

    delta = (x(i-1,:) - x(i-2,:));
    xGuess = x(i-1,:) + delta;
    [~,~,~,isMaxIterReached] = diffCorrec(xGuess,globalVar);
    while isMaxIterReached
        delta = delta/2;
        xGuess = x(i-1,:) + delta;
        [~,~,~,isMaxIterReached] = diffCorrec(xGuess,globalVar);
    end
    [t(i,1),x(i,:),~,~] = diffCorrec(xGuess,globalVar);
    jacobianConst(i,1) = jacobiValue3D(x(i,:),mu);
    if abs(jacobianConst(i)-jacobianConst(i-1)) < tol
        disp('\n Error: The jacobianConst difference is negligible. \n');
        isOrbitFound = 0;
        break;
    end
    i = i+1;
end
xLowerAmp(1,:) = x(end-1,:); tLowerAmp(1,1) = t(end-1,1); jacobianConstLowerAmp(1,1) = jacobianConst(end-1,1);
xUpperAmp(1,:) = x(end,:); tUpperAmp(1,1) = t(end,1); jacobianConstUpperAmp(1,1) = jacobianConst(end,1);
end

switch isOrbitFound
    case 1
        fprintf('\n===============================================\n')
        fprintf('Found IC for orbit with jacobianConst lesser than and greater than %d\n',globalVar.userInput.jacobianConst)
        fprintf('===============================================\n')


        tol = 1e-10;        % Why tolerance is inreased??
        jacobianConstNew = jacobianConstUpperAmp(1,1);
        while abs(jacobianConstNew-jacobianConstReq)>tol
            delta = (xUpperAmp(1,:) - xLowerAmp(1,:))/2;
            xGuess = xLowerAmp(1,:) + delta;

            [t,x,~] = diffCorrec(xGuess,globalVar);
            jacobianConstNew = jacobiValue3D(x,mu);
            % CHECK eNew with jacobianConstand update in xcorrec and everywhere
            if jacobianConstNew>jacobianConstReq
                xLowerAmp(1,:) = x;
                tLowerAmp(1,1) = t;
                jacobianConstLowerAmp(1,1) = jacobianConstNew;
                %idx = (length(tCorrec))-1;
            else
                xUpperAmp(1,:) = x;
                tUpperAmp(1,1) = t;
                jacobianConstUpperAmp(1,1) = jacobianConstNew;
                %idx = (length(tCorrec));
            end
        end

        [~,monodromy,~,~] = stm_X(globalVar,x,funVarEq,t);
        [eigens.val.stable,eigens.val.unstable,eigens.val.center,eigens.val.p,eigens.vec.stable,...
            eigens.vec.unstable,eigens.vec.center,eigens.vec.p] = calcEigen(monodromy,1) ;
        stabilityIdx = calcStabilityIdx(eigens);
        orbPar.period      = t; %(NoofFam x 1) - Full Orbit Time
        orbPar.IC          = x; %(NoofFam x UserDat.Dimension)
        orbPar.jacobianConst    = jacobianConstNew;
        orbPar.monodromy = monodromy;
        orbPar.eigens    = eigens;
        orbPar.stabilityIdx = stabilityIdx;
    case 0

        if jacobianConstReq> jacobianMax
            fprintf('\n===============================================\n')
            fprintf('Error: The jacobian(jacobianConst) value enetered is greater than the maximum limit for the Libration Point\n')
            fprintf('===============================================\n');
        else
            fprintf('\nError: jacobianConst input is less than the %d.\n',jacobianConst(end));
        end
        clear orbPar;
        orbPar = -1;
end

