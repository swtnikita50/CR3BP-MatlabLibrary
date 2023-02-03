% Incomplete
% Created on: 11-07-22 (13:21)

function [tCorrec,xCorrec,DF,isMaxIterReached] = pseudoArcDiffCorrec(xGuess, tGuess,globalVar,dels, X_prev,T_prev)
% Extract from Global Data
orbit   = globalVar.userInput.orbit;
f1      = globalVar.functions.systemDynamics;
f2      = globalVar.functions.varEq_stmDot;
tol 	= globalVar.userInput.tolerance^2;
fig1    = globalVar.userInput.diffCorrecPlot;

% arbitrary time span, should be long enough to reach the y crossing event
tspan = [0 10];

isMaxIterReached = 0;
maxIter = 20;
iter = 1; % set the iter value

if xGuess(3) == 0 && strcmp(orbit,'halo')
    orbit = 'lyapunov';
end


del_xDotf = 1;
del_zDotf = 1;
Fx = [del_xDotf;del_zDotf];
Gx = 1;
while  norm(Fx)>tol
    
    switch orbit
        case 'lyapunov'
            xFree = [xGuess(1), xGuess(5), tGuess]';
            xFreeG = [xGuess(1), xGuess(5)]';
            xFree_prev = [X_prev(1), X_prev(5), T_prev]';
            xFreeG_prev = [X_prev(1), X_prev(5)]';
            delX_prev = (xFreeG - xFreeG_prev)/dels;
        case 'halo'
            xFree = [xGuess(1);xGuess(5)];
    end
    % Check the max iter and stop
    if iter > maxIter
        fprintf('\nMaximum number of iterations exceded! \n');
        isMaxIterReached = 1;
        break;
    end

    % Obtain a Baseline Solution (This gives the variational end point which is
    % the constarint to be satisfied)
    [tb,xb] = integrate(globalVar,f1,xGuess,tspan,'crossing');

    % Constraint Vector
    switch orbit
        case 'lyapunov'
            del_yf = xb(end,2) ;
            del_xDotf =xb(end,4);
            del_yDotf = xb(end,5);
            Fx = [del_yf; del_xDotf];
            Gx = [Fx; (xFreeG - xFreeG_prev)'*delX_prev-dels];
        case 'halo'
            del_yf = xb(end,2) ;
            del_xDotf =xb(end,4) ;
            del_zDotf = xb(end,6);
            Fx = [del_xDotf;del_zDotf]; % Constraint vector
    end


    %[tbnew,xbnew] = integrate(globalVar,f1,xGuess,[0 2*tb(end)]);
    
    %% Plotting Correction Plot
    if globalVar.userInput.plotDiffCorrec == 1
        figure(fig1)
        ax = gca; ax.ColorOrderIndex = iter;
        plot(xb(:,1),xb(:,2))
        ax = gca; ax.ColorOrderIndex = iter;
        hold on
        plot(xb(1,1),xb(1,2),'*')
        ax = gca; ax.ColorOrderIndex = iter;
        plot(xb(end,1),xb(end,2),'o')
        ax = gca; ax.ColorOrderIndex = iter;
        axis tight
    end
    %%


    % Get final STM(State Transistion Matrix)
    [t,PHItf,x,xf] = stm_X(globalVar,xGuess,f2,tb(end));
    
    % Get the final Vector Field from the EOM
    X_DotDotf = CR3BP([],xb(end,:),globalVar.userInput.mu);

    switch orbit
        case 'lyapunov'
            x_Dotf = X_DotDotf(1);
            y_Dotf = X_DotDotf(2);
            x_DotDotf = X_DotDotf(4);
            y_DotDotf = X_DotDotf(5);
            
            %DF = [PHItf(4,1), PHItf(4,5)];
            DF = [PHItf(2,1), PHItf(2,5), y_Dotf;PHItf(4,1), PHItf(4,5), x_DotDotf];
            DG = [DF; delX_prev', 0];
            DGn = DG'/(DG*DG');
            %DFn = DF'/(DF*DF');
            xFree = xFree - (DGn*Gx);
            %xFree = xFree - DF\Fx;
            xGuess(1) = xFree(1); % Corrected yDot value
            xGuess(5) = xFree(2);
            tGuess = xFree(3);
        case 'halo'
            y_Dotf = X_DotDotf(2);
            x_DotDotf = X_DotDotf(4);
            z_DotDotf = X_DotDotf(6);
            DF = ([PHItf(4,1) PHItf(4,5);PHItf(6,1) PHItf(6,5)] - (1/y_Dotf)*[x_DotDotf;z_DotDotf]*[PHItf(2,1) PHItf(2,5)]);
            DFn = -DF\Fx;
            xFree = xFree + DFn;
            % Guess Update
            xGuess(1) = xFree(1);
            xGuess(5) = xFree(2);
    end
    iter = iter + 1;  % Update iter

end
xCorrec = xGuess;
tCorrec = 2*tb(end);
end
