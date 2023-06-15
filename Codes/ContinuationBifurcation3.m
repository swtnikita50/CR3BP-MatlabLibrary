function [xAns] = ContinuationBifurcation3(xLower, xUpper, globalVar)

tol = globalVar.userInput.tolerance*10^-1;
mu = globalVar.userInput.mu;
funVarEq = globalVar.functions.varEq_stmDot;
switch globalVar.userInput.orbit
    case 'lyapunov'
        for i = 1:length(xLower.period)
            if ~isempty(xLower.stabilityIdx(i).saddle)
                bifurcationType = 3;
                reqX = xLower;
            else
                bifurcationType = 4;
                reqX = xUpper;
            end
            while abs(reqX.stabilityIdx(i).saddle-1) > tol
                delta = 1/2*(xUpper.IC(i,:)-xLower.IC(i,:));
                xGuess = xLower.IC(i,:) + delta;
                [~,~,~,isMaxIterReached] = diffCorrec(xGuess,globalVar);
                while isMaxIterReached
                    delta = delta/2;
                    xGuess = xLower.IC(i,:) + delta;
                    [~,~,~,isMaxIterReached] = diffCorrec(xGuess,globalVar);
                end
                [tMid(i),xMid(i,:),~,~] = diffCorrec(xGuess,globalVar);
                [~,monodromy,~,~] = stm_X(globalVar,xMid(i,:),funVarEq,tMid(i));
                [eigens.val.stable,eigens.val.unstable,eigens.val.center,eigens.val.p,eigens.vec.stable,...
            eigens.vec.unstable,eigens.vec.center,eigens.vec.p] = calcEigen(monodromy,1) ;
                stabilityIdx = calcStabilityIdx(eigens);
                if ~isempty(stabilityIdx.saddle) && bifurcationType == 3 || isempty(stabilityIdx.saddle) && bifurcationType == 4
                    xLower.IC(i,:) = xMid(i,:);
                    xLower.period(i) = tMid(i);
                    xLower.jacobianConst(i) = jacobiValue3D(xMid(i,:),mu);
                    xLower.monodromy(:,:,i) = monodromy;
                    xLower.eigens(i) = eigens;
                    xLower.stabilityIdx(i) = stabilityIdx;

                else
                    xUpper.IC(i,:) = xMid(i,:);
                    xUpper.period(i) = tMid(i);
                    xUpper.jacobianConst(i) = jacobiValue3D(xMid(i,:),mu);
                    xUpper.monodromy(:,:,i) = monodromy;
                    xUpper.eigens(i) = eigens;
                    xUpper.stabilityIdx(i) = stabilityIdx;
                end
                if bifurcationType == 3
                    reqX = xLower;
                else
                    reqX = xUpper;
                end
            end
            xAns = reqX;
        end
    case 'halo'
        delX0 = [0 0 1 0 0 0];
        S = 0.005;
        delta = delX0*S;
        for i = 1:length(xLower.period)
            if ~isempty(xLower.stabilityIdx(i).saddle)
                bifurcationType = 3;
                reqX.IC = xLower.IC(i,:);
                reqX.period = xLower.period(i);
                reqX.jacobianConst = xLower.jacobianConst(i);
                reqX.monodromy = xLower.monodromy(:,:,i);
                reqX.eigens = xLower.eigens(i);
                reqX.stabilityIdx = xLower.stabilityIdx(i);
                
            else
                bifurcationType = 4;
                reqX.IC = xUpper.IC(i,:);
                reqX.period = xUpper.period(i);
                reqX.jacobianConst = xUpper.jacobianConst(i);
                reqX.monodromy = xUpper.monodromy(:,:,i);
                reqX.eigens = xUpper.eigens(i);
                reqX.stabilityIdx = xUpper.stabilityIdx(i);
            end
            while abs(reqX.stabilityIdx(i).saddle-1) > tol
                delta = 1/2*(xUpper.IC(i,:)-xLower.IC(i,:));
                xGuess = xLower.IC(i,:) + delta;
                [~,~,~,isMaxIterReached] = diffCorrec(xGuess,globalVar);
                while isMaxIterReached
                    delta = delta/2;
                    xGuess = xLower(i,:) + delta;
                    [~,~,~,isMaxIterReached] = diffCorrec(xGuess,globalVar);
                end
                [tMid(i,1),xMid(i,:),~,~] = diffCorrec(xGuess,globalVar);
                S = (xMid(i,3) - xLower.IC(i,3));
                [~,monodromy,~,~] = stm_X(globalVar,xMid,funVarEq,tMid);
                [eigens.val.stable,eigens.val.unstable,eigens.val.center,eigens.val.p,eigens.vec.stable,...
            eigens.vec.unstable,eigens.vec.center,eigens.vec.p] = calcEigen(monodromy,1) ;
                stabilityIdx = calcStabilityIdx(eigens);

                if ~isempty(stabilityIdx.saddle) && bifurcationType == 3 || isempty(stabilityIdx.saddle) && bifurcationType == 4
                    xLower.IC(i,:) = xMid(i,:);
                    xLower.period(i) = tMid(i);
                    xLower.jacobianConst(i) = jacobiValue3D(xMid(i,:),mu);
                    xLower.monodromy(:,:,i) = monodromy;
                    xLower.eigens(i) = eigens;
                    xLower.stabilityIdx(i) = stabilityIdx;

                else
                    xUpper.IC(i,:) = xMid(i,:);
                    xUpper.period(i) = tMid(i);
                    xUpper.jacobianConst(i) = jacobiValue3D(xMid(i,:),mu);
                    xUpper.monodromy(:,:,i) = monodromy;
                    xUpper.eigens(i) = eigens;
                    xUpper.stabilityIdx(i) = stabilityIdx;
                end
                if bifurcationType == 3
                    reqX.IC = xLower.IC(i,:);
                    reqX.period = xLower.period(i);
                    reqX.jacobianConst = xLower.jacobianConst(i);
                    reqX.monodromy = xLower.monodromy(:,:,i);
                    reqX.eigens = xLower.eigens(i);
                    reqX.stabilityIdx = xLower.stabilityIdx(i);
                else
                    reqX.IC = xUpper.IC(i,:);
                    reqX.period = xUpper.period(i);
                    reqX.jacobianConst = xUpper.jacobianConst(i);
                    reqX.monodromy = xUpper.monodromy(:,:,i);
                    reqX.eigens = xUpper.eigens(i);
                    reqX.stabilityIdx = xUpper.stabilityIdx(i);
                end
            end
            xAns.IC(i,:) =               reqX.IC(i,:);
            xAns.period(i) =           reqX.period(i);
            xAns.jacobianConst(i) =    reqX.jacobianConst(i);
            xAns.monodromy(:,:,i) =        reqX.monodromy(:,:,i);
            xAns.eigens(i) =           reqX.eigens(i);
            xAns.stabilityIdx(i) =     reqX.stabilityIdx(i);
        end

end