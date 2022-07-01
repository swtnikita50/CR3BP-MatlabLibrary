%{
...
Created on  24/2/2020 19:06

ODE Events function - to stop integration at required value

Inputs
------
1) t - time value
2) x - x values

Outputs
--------

for description , see about how to use ODE events

Dependencies
------------
None
...
%}
function [position,isterminal,direction] = Events_cres(t,x,G_var)


mu = G_var.Constants.mu;



isterminal = 1; 
position = x(2);

  if (x(1)<G_var.LagPts.L1(1)) || (x(1) < 1-mu &&  x(1) > G_var.LagPts.L1(1))
 direction = -1;
  elseif x(1)>G_var.LagPts.L2(1) || (x(1) > 1-mu &&  x(1) < G_var.LagPts.L2(1))
 direction = 1;
  elseif x(1)<G_var.LagPts.L3(1) || (x(1) < 0 &&  x(1) > G_var.LagPts.L3(1))
 direction = -1;
  end


   
 
 end