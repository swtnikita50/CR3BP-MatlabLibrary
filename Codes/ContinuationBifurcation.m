function [xAns] = ContinuationBifurcation(xLower, xUpper, G_var)

tol = 1e-4;
mu = G_var.Constants.mu;
funVarEq = G_var.IntFunc.VarEqAndSTMdot;
switch G_var.Constants.orbit
    case 'lyapunov'
        for i = 1:length(xLower.time)
            if xLower.StabilityIdx(i).C < 1
                bifurcationType = 1;
            else
                bifurcationType = 2;
            end
            while abs(xUpper.StabilityIdx(i).C-xLower.StabilityIdx(i).C) > tol
                delta = 1/2*(xUpper.IC(i,:)-xLower.IC(i,:));
                GuessX = xLower.IC(i,:) + delta;
                [~,~,~,isMaxIterReached] = DiffCorrec(GuessX,0,G_var);
                while isMaxIterReached
                    delta = delta/2;
                    GuessX = xLower.IC(i,:) + delta;
                    [~,~,~,isMaxIterReached] = DiffCorrec(GuessX,0,G_var);
                end
                [tMid(1),xMid(1,:),~,~] = DiffCorrec(GuessX,0,G_var);
                [~,Monodromy,~,~] = StateTransAndX(G_var,xMid,funVarEq,tMid);
                [Eigens.S_EigVal,Eigens.US_EigVal,Eigens.C_Val,Eigens.P_EigVal,Eigens.S_EigVec,...
                    Eigens.US_EigVec,Eigens.C_EigVec,Eigens.P_EigVec] = CalcEigenValVec(Monodromy,1) ;
                StabilityIdx = CalcStabilityIdx(Eigens);
                if StabilityIdx.C < 1 && bifurcationType == 1 || StabilityIdx.C > 1 && bifurcationType == 2
                    xLower.IC(i,:) = xMid(1,:);
                    xLower.time(i) = tMid;
                    xLower.Energy(i) = jacobiValue3D(xMid(i,:),mu);
                    xLower.Monodromy(:,:,i) = Monodromy;
                    xLower.Eigens(i) = Eigens;
                    xLower.StabilityIdx(i) = StabilityIdx;

                else
                    xUpper.IC(i,:) = xMid(1,:);
                    xUpper.time(i) = tMid;
                    xUpper.Energy(i) = jacobiValue3D(xMid(i,:),mu);
                    xUpper.Monodromy(:,:,i) = Monodromy;
                    xUpper.Eigens(i) = Eigens;
                    xUpper.StabilityIdx(i) = StabilityIdx;
                end
            end
            delta = 1/2*(xUpper.IC(i,:)-xLower.IC(i,:));
            GuessX = xLower.IC(i,:) + delta;
            [~,~,~,isMaxIterReached] = DiffCorrec(GuessX,0,G_var);
            while isMaxIterReached
                delta = delta/2;
                GuessX = xLower.IC(i,:) + delta;
                [~,~,~,isMaxIterReached] = DiffCorrec(GuessX,0,G_var);
            end
            [tMid(i),xMid(i,:),~,~] = DiffCorrec(GuessX,0,G_var);
            [~,Monodromy,~,~] = StateTransAndX(G_var,xMid(i,:),funVarEq,tMid);
            [Eigens.S_EigVal,Eigens.US_EigVal,Eigens.C_Val,Eigens.P_EigVal,Eigens.S_EigVec,...
                Eigens.US_EigVec,Eigens.C_EigVec,Eigens.P_EigVec] = CalcEigenValVec(Monodromy,1) ;
            StabilityIdx = CalcStabilityIdx(Eigens);
            xAns.Eigens(i) = Eigens;
            xAns.StabilityIdx(i,1) = StabilityIdx;
            xAns.time(i) = tMid(i);
            xAns.Monodromy(:,:,i) = Monodromy;
            xAns.IC(i,:) = xMid(i,:);
            xAns.Energy(i) = jacobiValue3D(xMid(i,:),mu);
        end
    case 'halo'
        delX0 = [0 0 1 0 0 0];
        S = 0.01;
        for i = 1:length(xLower.time)
            if xLower.StabilityIdx(i).C < 1
                bifurcationType = 1;
            else
                bifurcationType = 2;
            end
            while abs(xUpper.StabilityIdx(i).C-xLower.StabilityIdx(i).C) > tol
                delta = delX0*S;
                GuessX = xLower.IC(i,:) + delta;
                [~,~,~,isMaxIterReached] = DiffCorrec(GuessX,0,G_var,'halo');
                while isMaxIterReached
                    delta = delta/2;
                    GuessX = xLower(i,:) + delta;
                    [~,~,~,isMaxIterReached] = DiffCorrec(GuessX,0,G_var,'halo');
                end
                [tMid(i,1),xMid(i,:),~,~] = DiffCorrec(GuessX,0,G_var,'halo');
                Energy(i,1) = jacobiValue3D(xMid(i,:),mu);
                S = (xMid(i,3) - xLower.IC(i,3));
                [~,Monodromy,~,~] = StateTransAndX(G_var,xMid,funVarEq,tMid);
                [Eigens.S_EigVal,Eigens.US_EigVal,Eigens.C_Val,Eigens.P_EigVal,Eigens.S_EigVec,...
            Eigens.US_EigVec,Eigens.C_EigVec,Eigens.P_EigVec] = CalcEigenValVec(Monodromy,1) ;
        
                StabilityIdx = CalcStabilityIdx(Eigens);

                if StabilityIdx.C < 1 && bifurcationType == 1 || StabilityIdx.C > 1 && bifurcationType == 2
                    xLower.IC(i,:) = xMid(1,:);
                    xLower.time(i) = tMid;
                    xLower.Energy(i) = jacobiValue3D(xMid(i,:),mu);
                    xLower.Monodromy(:,:,i) = Monodromy;
                    xLower.Eigens(i) = Eigens;
                    xLower.StabilityIdx(i) = StabilityIdx;

                else
                    xUpper.IC(i,:) = xMid(1,:);
                    xUpper.time(i) = tMid;
                    xUpper.Energy(i) = jacobiValue3D(xMid(i,:),mu);
                    xUpper.Monodromy(:,:,i) = Monodromy;
                    xUpper.Eigens(i) = Eigens;
                    xUpper.StabilityIdx(i) = StabilityIdx;
                end
            end
            delta = 1/2*(xUpper.IC(i,:)-xLower.IC(i,:));
            GuessX = xLower.IC(i,:) + delta;
            [~,~,~,isMaxIterReached] = DiffCorrec(GuessX,0,G_var);
            while isMaxIterReached
                delta = delta/2;
                GuessX = xLower.IC(i,:) + delta;
                [~,~,~,isMaxIterReached] = DiffCorrec(GuessX,0,G_var);
            end
            [tMid(i),xMid(i,:),~,~] = DiffCorrec(GuessX,0,G_var);
            [~,Monodromy,~,~] = StateTransAndX(G_var,xMid(i,:),funVarEq,tMid);
            [Eigens.S_EigVal,Eigens.US_EigVal,Eigens.C_Val,Eigens.P_EigVal,Eigens.S_EigVec,...
                Eigens.US_EigVec,Eigens.C_EigVec,Eigens.P_EigVec] = CalcEigenValVec(Monodromy,1) ;
            StabilityIdx = CalcStabilityIdx(Eigens);
            xAns.Eigens(i) = Eigens;
            xAns.StabilityIdx(i,1) = StabilityIdx;
            xAns.time(i) = tMid(i);
            xAns.Monodromy(:,:,i) = Monodromy;
            xAns.IC(i,:) = xMid(i,:);
            xAns.Energy(i) = jacobiValue3D(xMid(i,:),mu);
        end

end