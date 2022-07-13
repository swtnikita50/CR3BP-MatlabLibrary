% Incomplete
% Created on: 11-07-22 (13:21)

function [tCorrec,xCorrec,DF,isMaxIterReached] = pseudoArcDiffCorrec(X_Guess,dels,Plot,G_var,orbitType)
% Extract from Global Data

if nargin<5
    orbitType = 'lyapunov';
end

f1 = G_var.IntFunc.EOM;
f2 = G_var.IntFunc.VarEqAndSTMdot;
tspan = [0 10];

isMaxIterReached = 0;
MaxIteration = 10;
iteration = 1; % set the iteration value

switch orbitType
    case 'lyapunov'
        Xfree = X_Guess(5);
    case 'halo'
        Xfree = [X_Guess(1);X_Guess(5)];
end
del_xDotf = 1;
del_zDotf = 1;
Fx = [del_xDotf;del_zDotf];
while  norm(Fx)>1.e-10
    ax = gca; ax.ColorOrderIndex = iteration;
    
% Check the max iteration and stop
if iteration > MaxIteration
    fprintf('\nMaximum number of iterations exceded: Check for errors in script. \n');
    isMaxIterReached = 1;
    break;
end

fprintf('\nDiffCorrecIter# %d',iteration)

% Obtain a Baseline Solution (This gives the variational end point which is
% the constarint to be satisfied)
[tb,xb] = Integrator(G_var,f1,X_Guess,tspan,'Events');

% Develop the constraint vector

switch orbitType
    case 'lyapunov'
        del_yf = xb(end,2) ;
        del_xDotf =xb(end,4);
        del_yDotf = xb(end,5);
        Fx = del_xDotf;
    case 'halo'
        del_yf = xb(end,2) ;
        del_xDotf =xb(end,4) ;
        del_zDotf = xb(end,6);
        Fx = [del_xDotf;del_zDotf]; % Constraint vector
end


% Plot the correction(Just to see how differential correction works) If
% Required - adopted from Prof SD Ross code 

[tbnew,xbnew] = Integrator(G_var,f1,X_Guess,[0 2*tb(end)]);
if Plot == 1
    switch orbitType
    case 'lyapunov'
        plot(xb(:,1),xb(:,2))
        ax = gca; ax.ColorOrderIndex = iteration;
        hold on
        plot(xb(1,1),xb(1,2),'*')
        ax = gca; ax.ColorOrderIndex = iteration;
        plot(xb(end,1),xb(end,2),'o')
        ax = gca; ax.ColorOrderIndex = iteration;
        axis tight
    case 'halo'
        plot3(xb(:,1),xb(:,2),xb(:,3))
        ax = gca; ax.ColorOrderIndex = iteration;

        hold on; grid on
        %plot3(xbnew(:,1),xbnew(:,2),xbnew(:,3))
        %ax = gca; ax.ColorOrderIndex = iteration;
        plot3(xb(1,1),xb(1,2),xb(1,3),'*')
        ax = gca; ax.ColorOrderIndex = iteration;
        plot3(xb(end,1),xb(end,2),xb(end,3),'o')
        plot3(1-G_var.Constants.mu,0,0,'p');
        ax = gca; ax.ColorOrderIndex = iteration;
        axis tight
end

end


% Now get the STM(State Transistion Matrix)

[t,PHItf,x,xf] = StateTransAndX(G_var,X_Guess,f2,tb(end));
% Get the final Vector Field from the EOM



X_DotDotf = CRes3BP_EOM([],xb(end,:),G_var.Constants.mu);

switch orbitType
    case 'lyapunov'
        y_Dotf = X_DotDotf(2);
        x_DotDotf = X_DotDotf(4);
        DF = (PHItf(4,5)*y_Dotf - PHItf(2,5)*x_DotDotf)/y_Dotf;
        Xfree = Xfree - (1/DF)*Fx;
        X_Guess(5) = Xfree; % Corrected yDot value
    case 'halo'
        y_Dotf = X_DotDotf(2);
        x_DotDotf = X_DotDotf(4);
        z_DotDotf = X_DotDotf(6);
        DF = ([PHItf(4,1) PHItf(4,5);PHItf(6,1) PHItf(6,5)] - (1/y_Dotf)*[x_DotDotf;z_DotDotf]*[PHItf(2,1) PHItf(2,5)]);
        DFn = -DF\Fx;
        Xfree = Xfree + DFn;
        % Guess Update
        X_Guess(1) = Xfree(1);
        X_Guess(5) = Xfree(2);
end


iteration = iteration + 1;  % Update Iteration   

end
xCorrec = X_Guess;
tCorrec = 2*tb(end);
end