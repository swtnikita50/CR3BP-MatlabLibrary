% Incomplete
% Pseudo Arc-Length Continuation Method
% Created on: 11-07-22 (13:13)

function X_new = pseudoArcLengthCont(X_prev,dels,globalVar)
f1 = globalVar.functions.systemDynamics;
f2 = globalVar.functions.varEq_stmDot;
tol = globalVar.userInput.tolerance^2;
mu = globalVar.userInput.mu;
tspan = [0 10];

[tb,xb] = integrate(globalVar,f1,X_prev,tspan,'crossing');
[~,~,~,~, PHI] = stm_X(globalVar, X_prev,f2,tb(end));
X_DotDotf = CR3BP([],xb(end,:),globalVar.userInput.mu);


T_prev = tb(end);

switch globalVar.userInput.orbit
    case 'lyapunov'
        PHItf = reshape(PHI(end,1:36),6,6);
        DF = [PHItf(4,1), PHItf(4,5)];
        var = null(DF);
        if var(1,1) > 0
            var = -var;
        end
        delX_prev = [var(1, 1); 0; 0; 0; var(2,1); 0]';
        X_new = X_prev + dels*delX_prev;
            
    case 'halo'
        PHItf = reshape(PHI(end,1:36),6,6);
        DF = [PHItf(2,:); PHItf(4,:); PHItf(6,:)];
        
        var = null(DF);
        switch globalVar.userInput.type
            case 'northern'
                for i = 1:6
                    if var(i,2) > 0
                        delX_prev = [ var(i, 1); 0; var(i,2); 0; var(i,3); 0]';
                        X_new = X_prev + dels*delX_prev;
                        break
                    end
                end
            case 'southern'
                for i = 1:6
                    if var(i,2) < 0
                        delX_prev = [ var(i, 1); 0; var(i,2); 0; var(i,3); 0]';
                        X_new = X_prev + dels*delX_prev;
                        break
                    end
                end
        end
end

[tb,xb] = integrate(globalVar,f1,X_new,tspan,'crossing');
T_new = tb(end);
[tCorrec,xCorrec,~,isMaxIterReached] = pseudoArcDiffCorrec(X_new,T_new, globalVar,dels, X_prev,T_prev);

while isMaxIterReached
    dels = dels/2;
    X_new = X_prev + dels*delX_prev;
    [tCorrec,xCorrec,~,isMaxIterReached] = pseudoArcDiffCorrec(X_new,T_new, globalVar,dels, X_prev,T_prev);
end
[~,x] = integrate(globalVar,globalVar.functions.systemDynamics,xCorrec,[0 tCorrec],'forward');
figure(5)
plot3(x(:,1),x(:,2),x(:,3));hold on; grid on;
scatter3(globalVar.lagPts.pos(globalVar.userInput.lagrangePt,1),globalVar.lagPts.pos(globalVar.userInput.lagrangePt,2),0,'p','filled');


X_new = xCorrec;
end


