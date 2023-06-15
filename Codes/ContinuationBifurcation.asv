function [xAns] = ContinuationBifurcation(xLower, xUpper, globalVar)

tol = globalVar.userInput.tolerance*10^-1;
mu = globalVar.userInput.mu;
funVarEq = globalVar.functions.varEq_stmDot;
switch globalVar.userInput.orbit
    case 'lyapunov'
        for i = 1:length(xLower.period)
            if xLower.stabilityIdx(i).center < 1
                bifurcationType = 1;    %Halo
            else
                bifurcationType = 2;    %Axial
            end
            while abs(xUpper.stabilityIdx(i).center-xLower.stabilityIdx(i).center) > tol
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
                if stabilityIdx.center < 1 && bifurcationType == 1 || stabilityIdx.center > 1 && bifurcationType == 2
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
            end
            delta = 1/2*(xUpper.IC(i,:)-xLower.IC(i,:));
            xGuess = xLower.IC(i,:) + delta;
            [~,~,~,isMaxIterReached] = diffCorrec(xGuess,globalVar);
            while isMaxIterReached
                delta = delta/2;
                xGuess = xLower.IC(i,:) + delta;
                [~,~,~,isMaxIterReached] = diffCorrec(xGuess,globalVar);
            end
            [tMid(i),xMid(i,:),~,~] = diffCorrec(xGuess,globalVar);
            [~,monodromy,~,~] = stm_X(globalVar,xMid(i,:),funVarEq,tMid);
            [eigens.val.stable,eigens.val.unstable,eigens.val.center,eigens.val.p,eigens.vec.stable,...
            eigens.vec.unstable,eigens.vec.center,eigens.vec.p] = calcEigen(monodromy,1) ;
                stabilityIdx = calcStabilityIdx(eigens);
            xAns.eigens(i) = eigens;
            xAns.stabilityIdx(i,1) = stabilityIdx;
            xAns.period(i) = tMid(i);
            xAns.monodromy(:,:,i) = monodromy;
            xAns.IC(i,:) = xMid(i,:);
            xAns.jacobianConst(i) = jacobiValue3D(xMid(i,:),mu);
        end
    case 'halo'
        delX0 = [0 0 1 0 0 0];
        S = 0.01;
        for i = 1:length(xLower.period)
            if xLower.stabilityIdx(i).center < 1
                bifurcationType = 1;
            else
                bifurcationType = 2;
            end
            while abs(xUpper.stabilityIdx(i).center-xLower.stabilityIdx(i).center) > tol
                delta = delX0*S;
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

                if stabilityIdx.center < 1 && bifurcationType == 1 || stabilityIdx.center > 1 && bifurcationType == 2
                    xLower.IC(i,:) = xMid(1,:);
                    xLower.period(i) = tMid;
                    xLower.jacobianConst(i) = jacobiValue3D(xMid(i,:),mu);
                    xLower.monodromy(:,:,i) = monodromy;
                    xLower.eigens(i) = eigens;
                    xLower.stabilityIdx(i) = stabilityIdx;

                else
                    xUpper.IC(i,:) = xMid(1,:);
                    xUpper.period(i) = tMid;
                    xUpper.jacobianConst(i) = jacobiValue3D(xMid(i,:),mu);
                    xUpper.monodromy(:,:,i) = monodromy;
                    xUpper.eigens(i) = eigens;
                    xUpper.stabilityIdx(i) = stabilityIdx;
                end
            end
            delta = 1/2*(xUpper.IC(i,:)-xLower.IC(i,:));
            xGuess = xLower.IC(i,:) + delta;
            [~,~,~,isMaxIterReached] = diffCorrec(xGuess,globalVar);
            while isMaxIterReached
                delta = delta/2;
                xGuess = xLower.IC(i,:) + delta;
                [~,~,~,isMaxIterReached] = diffCorrec(xGuess,globalVar);
            end
            [tMid(i),xMid(i,:),~,~] = diffCorrec(xGuess,globalVar);
            [~,monodromy,~,~] = StateTransAndX(globalVar,xMid(i,:),funVarEq,tMid);
            [eigens.val.stable,eigens.val.unstable,eigens.val.center,eigens.val.p,eigens.vec.stable,...
            eigens.vec.unstable,eigens.vec.center,eigens.vec.p] = calcEigen(monodromy,1) ;
                stabilityIdx = calcStabilityIdx(eigens);
            xAns.eigens(i) = eigens;
            xAns.stabilityIdx(i,1) = stabilityIdx;
            xAns.period(i) = tMid(i);
            xAns.monodromy(:,:,i) = monodromy;
            xAns.IC(i,:) = xMid(i,:);
            xAns.jacobianConst(i) = jacobiValue3D(xMid(i,:),mu);
        end

end