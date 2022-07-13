%{
...
Created on 24/2/2020 16:11

This File is used for Integration uses "ODE113" - which uses "Adam
BashForth Moulton" Integration method

Inputs
1) G_var    - Requisite Global data
2) fun      - Function to be integrated.
3) x0       - Initial Condition vector
4) tspan    - Time duration of integration(Refer MATLAB's help for ODE's on
   supplying tspan)
4) type     - 'forward' or 'backward' integration.

Outputs
Same as the outputs of ODE113, 
1) t - Time values
2) x - Integrated xValues.
 
Note
-----
This file can also be used for "ODE events" function , when calling call will the options
"XYZ.IntFunc.ODEoptionsEvents" from "GlobalData". the events function is explicit for 
Circular Restrcited Three Body Problem. Change if you want.
 
Dependencies
------------
1) MATLAB inbuilt ODE113 integrator
%}
function [t,x] = Integrator(G_var,fun,x0,tspan,type)


if nargin<3
disp('provide function and/or initial condition')
elseif nargin == 3
     tspan = [0 10];
     type = 'forward';
     
elseif nargin == 4
      type = 'forward';  
end

switch type
    case 'forward'
           options  = G_var.IntFunc.ODEoptions;
           [t,x]    = ode113(fun,tspan,x0,options);
  
    case 'backward'
            options = G_var.IntFunc.ODEoptions;
            [t,x]   = ode113(fun,-tspan,x0,options);
    case 'Events'
            EventFunc = @(t,x) Events_cres(t,x,G_var);
            options = odeset('Reltol',1e-12,'Abstol',1e-12,'Events',EventFunc);
            [t,x]=ode113(fun,tspan,x0,options);
end
end