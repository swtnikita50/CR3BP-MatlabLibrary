%{
...
Created on  24/2/2020 18:41

This file hands over the state transistion matrix (STM)

Inputs
------
1) G_var    - Requisite Global data
2) X_guess  - Guess vector (1 x 6) or (1 x 4)
3) fun      - function to be integrated(Actually the variational equation)
4) tend     - Stop time of integration


Outputs
--------
for all families
1) t - time steps of integrations(nsteps x 1), nsteps - decided by ODE
2) PHItf - Final STM at the end of trajectory
3) x - states - (nsteps x 42) or (nsteps x 20) vector after integration 
4) xf - Final integrated vector (1 x 6) or (1 x 4) => final state after integration
Dependencies
------------
1) VarEq_Init(X_Guess)
2) Integrator(fun,x0,[0 tspan],options);
 
...
%}
function [t,PHItf,x,xf,PHI] = stm_X(globalVar,xGuess,fun,tend)

varEq_initVal = varEq_initialization(xGuess);
[tInteg,xInteg] = integrate(globalVar,fun,varEq_initVal,[0 tend]);

N = size(xInteg,2);
if N == 20
    type = 'Planar';
else 
    type = 'ThreeDim';
end

switch type
    case 'Planar'
        t = tInteg;
        PHI = xInteg;
        PHItf = reshape(xInteg(end,1:16),4,4);
        x = xInteg(:,17:20);
        xf = xInteg(end,17:20);

    case'ThreeDim'
        t = tInteg;
        PHI = xInteg;
        PHItf = reshape(xInteg(end,1:36),6,6);
        x = xInteg(:,37:42);
        xf = xInteg(end,37:42);
end

