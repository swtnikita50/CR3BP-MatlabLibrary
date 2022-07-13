%{
...
Created on  30/6/2022 17:42
taken from Karthi and then modified.

This File performs the Diffrential Correction Process.

Inputs
------
1) X_Guess - initial guess Calculated by the 'InitialGuess' function file.
2) Plot - User Supplied in MAINScript(will be 0/1)-See description in
   MAIN_LyapOrbit.m

Outputs
--------
1) tCorrec - Corrected half time value for perpendicular XZ plane crossing
2) xCorrec - Corrected Initial Value to obtain the periodic orbit.
3) DF - Final DF matrix(this is just for the user to see)
4) isMaxIterReached - 1 if No of iterations = MaxIterations is reached, 0 otheriwse

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

function [tCorrec,xCorrec,DF,isMaxIterReached] = DiffCorrec(X_Guess,Plot,G_var,orbitType)
% Extract from Global Data

if nargin<4
    orbitType = 'lyapunov';
end

f1 = G_var.IntFunc.EOM;
f2 = G_var.IntFunc.VarEqAndSTMdot;
tspan = [0 10];

isMaxIterReached = 0;
MaxIteration = 10;
iteration = 1; % set the iteration value

if X_Guess(3) == 0 && strcmp(orbitType,'halo')
    orbitType = 'lyapunov';
end

tol = 1e-12;
switch orbitType
    case 'lyapunov'
        Xfree = X_Guess(5);
    case 'halo'
        Xfree = [X_Guess(1);X_Guess(5)];
end
del_xDotf = 1;
del_zDotf = 1;
Fx = [del_xDotf;del_zDotf];
while  norm(Fx)>tol
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
    
