%{
...
Created on Feb 20, 2020 16:18 

This function calculates the jacobi value of a state in the Circular
Restricted 3-Body Problem (CR3BP)

Inputs
------
   X  - 1x6 Vector - Your current states [x,y,z,xDot,yDot,zDot]
   mu - double     - Your systems nondimenisonal mass ratio

Outputs
-------
   jacobiConst - double - The jacobi constant
   C = 2U - (xdot^2+ydot^2+zdot^2)
   U = (x^2 +y^2)/2 + (1-mu)/d + mu/r

=======================================================
Source File: Ari Rubinstein
WebSite    : https://github.com/gereshes/CR3BP-Functions/blob/master/jacobiValue3D.m
Modified to accept 2D and 3D value 
=======================================================


...
%}
function [jacobiConst] = jacobiValue3D(X,mu)


if length(X)>4
type = 'ThreeDim';
else
    type = 'Planar';
end


switch type
    case 'Planar'
        x=X(1);
        y=X(2);
       
        xDot=X(3);
        yDot=X(4);
        

r=sqrt(((x-1+mu).^2)+(y.^2));
d=sqrt(((x+mu).^2)+(y.^2));
jacobiConst = (x.^2)+(y.^2) +(2*(1-mu)./d)+(2*mu./r) - ((xDot.^2)+(yDot.^2));

%============================================================================

    case 'ThreeDim'
x=X(1);
y=X(2);
z=X(3);
xDot=X(4);
yDot=X(5);
zDot=X(6);

r=sqrt(((x-1+mu).^2)+(y.^2)+(z.^2));
d=sqrt(((x+mu).^2)+(y.^2)+(z.^2));
jacobiConst = (x.^2)+(y.^2) +(2*(1-mu)./d)+(2*mu./r) - ((xDot.^2)+(yDot.^2)+(zDot.^2));
    
     

 
        
end
