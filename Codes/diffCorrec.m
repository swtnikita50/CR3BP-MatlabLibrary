%{
...
Created on  30/6/2022 17:42
taken from Karthi and then modified.

Applies differential correction on guess solutions

Inputs
------
1) xGuess - initial guess Calculated by the 'InitialGuess' function file.
2) globalVar - Global Variable

Outputs
--------
1) tCorrec - Corrected half time value for perpendicular XZ plane crossing
2) xCorrec - Corrected Initial Value to obtain the periodic orbit.
3) DF - Final DF matrix(this is just for the user to see)
4) isMaxIterReached - 1 if No of iters = maxIters is reached, 0 otheriwse

References
----------
1) Wang Sang Koon, Martin W Lo,JE Marsden, Shane D Ross - "Dynamical
   Systems,the Three Body Problem and Space Mission Design", 2011
2) Tatiana Mar Vaquero Escribano - "Poincare Sections And Resonant Orbits In 
   the Restricted Three Body Problem",MS Thesis,Purdue University,2010.
   (Chapter 2 - Sec 2.3 and 2.4).
3) Thomas A. Pavlak - "Mission Design Applications In the Earth Moon
   System:Transfer Trajectories And StationKeeping",MS Thesis, Purdue
   University, 2010 May
   (Chapter 2 - Sec 2.3 and 2.4)

Note: These references DONOT have different methods , underlying concept is Newton's
method, refrences are just for convinience, one can follow any source. 
...
%}

function [tCorrec,xCorrec,DF,isMaxIterReached] = diffCorrec(xGuess,globalVar)
% Extract from Global Data
orbit = globalVar.userInput.orbit;
f1 = globalVar.functions.systemDynamics;
f2 = globalVar.functions.varEq_stmDot;
tol = globalVar.userInput.tolerance^2;
fig1 = globalVar.userInput.diffCorrecPlot;

% arbitrary time span, should be long enough to reach the y crossing event
tspan = [0 10];

isMaxIterReached = 0;
maxIter = 10;
iter = 1; % set the iter value

if xGuess(3) == 0 && strcmp(orbit,'halo')
    orbit = 'lyapunov';
elseif xGuess(6) == 0 && strcmp(orbit,'axial')
    orbit = 'lyapunov';
end
del_xDotf = 1;
del_zDotf = 1;
del_zf = 1;
switch orbit
    case 'lyapunov'
        xFree = xGuess(5);
        Fx = [del_xDotf;del_zDotf];
    case 'halo'
        xFree = [xGuess(1);xGuess(5)];
        Fx = [del_xDotf;del_zDotf];
    case 'halo2'
        xFree = [xGuess(3);xGuess(5)];
        Fx = [del_xDotf;del_zDotf];
    case 'axial'
        xFree = [xGuess(1); xGuess(5)];
        Fx = [del_zf; del_xDotf]; % Constraint vector
        maxIter = 20;
        tol = globalVar.userInput.tolerance*10^-2;
    case 'axial2'
        xFree = [xGuess(1); xGuess(6)];
        Fx = [del_zf; del_xDotf]; % Constraint vector
        maxIter = 20;
        tol = globalVar.userInput.tolerance*10^-2;
end


while  norm(Fx)>tol
    % Check the max iter and stop
    if iter > maxIter
        fprintf('\nMaximum number of iterations exceded! \n');
        isMaxIterReached = 1;
        break;
    end

    % Obtain a Baseline Solution (This gives the variational end point which is
    % the constarint to be satisfied)
    [tb,xb] = integrate(globalVar,f1,xGuess,tspan,'crossing');

    % Develop the constraint vector
    switch orbit
        case 'lyapunov'
            del_yf = xb(end,2) ;
            del_xDotf =xb(end,4);
            del_yDotf = xb(end,5);
            Fx = del_xDotf;
        case {'halo','halo2'}
            del_yf = xb(end,2) ;
            del_xDotf =xb(end,4) ;
            del_zDotf = xb(end,6);
            Fx = [del_xDotf;del_zDotf]; % Constraint vector
        case {'axial', 'axial2'}
            del_yf = xb(end,2) ;
            del_xDotf =xb(end,4) ;
            del_zf = xb(end,3);
            Fx = [del_zf; del_xDotf]; % Constraint vector
    end


    % Plot the correction(Just to see how differential correction works)
    %[tbnew,xbnew] = integrate(globalVar,f1,xGuess,[0 2*tb(end)]);
    if globalVar.userInput.plotDiffCorrec == 1
        figure(fig1)
        ax = gca; ax.ColorOrderIndex = iter;
        switch orbit
            case 'lyapunov'
                plot(xb(:,1),xb(:,2))
                ax = gca; ax.ColorOrderIndex = iter;
                hold on
                plot(xb(1,1),xb(1,2),'*')
                ax = gca; ax.ColorOrderIndex = iter;
                plot(xb(end,1),xb(end,2),'o')
                ax = gca; ax.ColorOrderIndex = iter;
                axis tight
            case {'halo', 'axial', 'halo2', 'axial2'}
                %[tb,xb] = integrate(globalVar,f1,xGuess,2*tspan);

                plot3(xb(:,1),xb(:,2),xb(:,3))
                ax = gca; ax.ColorOrderIndex = iter;

                hold on; grid on
                plot3(xb(1,1),xb(1,2),xb(1,3),'*')
                ax = gca; ax.ColorOrderIndex = iter;
                plot3(xb(end,1),xb(end,2),xb(end,3),'o')
                %plot3(1-globalVar.userInput.mu,0,0,'p');
                ax = gca; ax.ColorOrderIndex = iter;
                axis tight

        end

    end
[tb,xb] = integrate(globalVar,f1,xGuess,tspan,'crossing');


    % Now get the STM(State Transistion Matrix)

    [t,PHItf,x,xf] = stm_X(globalVar,xGuess,f2,tb(end));
    
    % Get the final Vector Field from the EOM
    X_DotDotf = CR3BP([],xb(end,:),globalVar.userInput.mu);

    switch orbit
        case 'lyapunov'
            y_Dotf = X_DotDotf(2);
            x_DotDotf = X_DotDotf(4);
            DF = (PHItf(4,5)*y_Dotf - PHItf(2,5)*x_DotDotf)/y_Dotf;
            xFree = xFree - (1/DF)*Fx;
            xGuess(5) = xFree; % Corrected yDot value
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
        case 'halo2'
            y_Dotf = X_DotDotf(2);
            x_DotDotf = X_DotDotf(4);
            z_DotDotf = X_DotDotf(6);
            DF = ([PHItf(4,3) PHItf(4,5);PHItf(6,3) PHItf(6,5)] - (1/y_Dotf)*[x_DotDotf;z_DotDotf]*[PHItf(2,3) PHItf(2,5)]);
            DFn = -DF\Fx;
            xFree = xFree + DFn;
            % Guess Update
            xGuess(3) = xFree(1);
            xGuess(5) = xFree(2);
        case 'axial'
            y_Dotf = X_DotDotf(2);
            x_DotDotf = X_DotDotf(4);
            z_Dotf = X_DotDotf(3);
            DF = ([PHItf(3,1) PHItf(3,5);PHItf(4,1) PHItf(4,5)] - (1/y_Dotf)*[z_Dotf;x_DotDotf]*[PHItf(2,1) PHItf(2,5)]);
            DFn = -DF\Fx;
            xFree = xFree + DFn;
            % Guess Update
            xGuess(1) = xFree(1);
            xGuess(5) = xFree(2);
        case 'axial2'
            y_Dotf = X_DotDotf(2);
            x_DotDotf = X_DotDotf(4);
            z_Dotf = X_DotDotf(3);
            DF = ([PHItf(3,1) PHItf(3,6);PHItf(4,1) PHItf(4,6)] - (1/y_Dotf)*[z_Dotf;x_DotDotf]*[PHItf(2,1) PHItf(2,6)]);
            DFn = -DF\Fx;
            xFree = xFree + DFn;
            % Guess Update
            xGuess(1) = xFree(1);
            xGuess(6) = xFree(2);
    end
    iter = iter + 1;  % Update iter

end
xCorrec = xGuess;
tCorrec = 2*tb(end);
end

