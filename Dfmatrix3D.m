%{
...
 
 Gives matrix Df(x) (i.e., the matrix of derivatives of f, where xdot=f(x) )
 for a 4-dimensional point x in the planar CR3BP's phase space 

-----------------------------------------------------------------------
 CR3BP (Circular Restricted Three-Body [Gravitational] Problem)
 with the LARGER MASS, M1 to the left of the origin at (-mu,0)
 and the smaller mass, M2, or the planet (ie. Earth), is at (1 - mu, 0)

       (rotating coords)

                 L4

 -L3------M1--+-----L1--M2--L2-

                 L5

 The following is the Jacobian matrix
 
 Df    =[  0     0    1    0 ;
           0     0    0    1 ;
       	 -Uxx  -Uxy   0    2 ;
         -Uxy  -Uyy  -2    0 ];

====================================================
 Source File       : Prof.Shane Ross (revised 2.19.04)
 WebSite(Codes)    : http://www.dept.aoe.vt.edu/~sdross/

Modification
1) Added a switch-case to calculate for planar and 3-D verrsions on 
 Feb 20,2020 by 15:30
====================================================

Inputs:
-------
1) mu - Mass parameter
2) x - X vector 

Outputs
-------
1) Df - Jacobian Matrix 
...
%}
function Df = Dfmatrix3D(x,mu)

if length(x)>4
type = 'ThreeDim';
else
    type = 'Planar';
end
% mu = mass paramater

 mu1 = 1-mu ;
 mu2 =   mu ;
 
switch type
    case 'Planar'
 r2= (x(1)+mu2)^2 + x(2)^2;      % r: distance to m1, LARGER MASS
 R2= (x(1)-mu1)^2 + x(2)^2;      % R: distance to m2, smaller mass
 
 r3= r2^1.5;
 r5= r2^2.5;
 R3= R2^1.5;
 R5= R2^2.5;
 
 %The following are three double partial derivatives of the
 % effective potential U(x,y)
 % First aprtial currently unused
Ux =  x(1) - mu1*(x(1)+mu2)/r3 - mu2*(x(1)-mu1)/R3 ;
Uy =  x(2) - mu1* x(2)     /r3 - mu2* x(2)     /R3 ;


 Uxx = 1 - mu1/r3 + (3*mu1*(x(1)+mu)^2)/r5 - mu2/R3 + (3*mu*(x(1)-1+mu)^2)/R5 ;
 Uyy = 1 - mu1/r3 + (3*mu1*x(2)^2)/r5 - mu2/R3 + (3*mu2*x(2)^2)/R5 ;
 Uxy = (3*mu1*(x(1)+mu2)*x(2))/r5 + (3*mu2*(x(1)-1+mu2)*x(2))/R5 ; 
 Uyx = Uxy;
 
 

I   = eye(2);
zer = zeros(2,2);
UXX = [Uxx Uxy; Uyx Uyy];
sig = [0 2;-2 0];

Df  = [zer I;UXX sig];
   
%==================================================================================
    case 'ThreeDim'
        
r2= (x(1)+mu2)^2 + x(2)^2 + x(3)^2;      % r: distance to m1, LARGER MASS
R2= (x(1)-mu1)^2 + x(2)^2 + x(3)^2;      % R: distance to m2, smaller mass

r3= r2^1.5;
r5= r2^2.5;
R3= R2^1.5;
R5= R2^2.5;

% The following are the two partial derivatives of the
% effective potential U(x,y)

% This is first partial currently unused
Ux =  x(1) - mu1*(x(1)+mu2)/r3 - mu2*(x(1)-mu1)/R3 ;
Uy =  x(2) - mu1* x(2)     /r3 - mu2* x(2)     /R3 ;
Uz = -(mu1*x(3))/r3 - (mu*x(3))/R3;    
    
    
    
    
Uxx = 1 - mu1/r3 + (3*mu1*(x(1)+mu)^2)/r5 - mu2/R3 + (3*mu*(x(1)-1+mu)^2)/R5 ;
Uyy = 1 - mu1/r3 + (3*mu1*x(2)^2)/r5 - mu2/R3 + (3*mu2*x(2)^2)/R5 ;
Uzz = -mu1/r3 + (3*mu1*x(3)^2)/r5 - mu2/R3 + (3*mu2*x(3)^2)/R5 ;
Uxy = (3*mu1*(x(1)+mu2)*x(2))/r5 + (3*mu2*(x(1)-1+mu2)*x(2))/R5 ; 
Uyx = Uxy;
Uxz = (3*mu1*(x(1)+mu2)*x(3))/r5 + (3*mu2*(x(1)-1+mu2)*x(3))/R5 ;
Uzx = Uxz ;
Uyz = (3*mu1*x(2)*x(3))/r5 + (3*mu2*x(2)*x(3))/R5 ;
Uzy = Uyz ;

I   = eye(3);
zer = zeros(3,3);
UXX = [Uxx Uxy Uxz; Uyx Uyy Uyz;Uzx Uzy Uzz];
sig = [0 2 0;-2 0 0;0 0 0];

Df  = [zer I;UXX sig];
end
