%{
...
Created on   26/2/2020 17:53

This file calculates the Initial Guess for Differential Correction.

Inputs
------
1) PointLoc - Lagrange Points - Provide (1/2/3) as inputs


Outputs
--------
for all families
1) XGuess(1 and 2) - Guess for the initial point, a structure file containing
        * XGuess.one - First guess
        * XGuess.two - Second Guess 
Note: 
XGuess(1) - Guess in 3 Dimension(1 x 6) vector
XGuess(2) - Guess in 2 Dimension(1 x 4) vector

Dependencies
------------
1) GlobalData - Takes "mu,Ax1,Ax2" for calculation.

Note: Change Ax1 and Ax2 values in "Globaldata"

Reference
----------
1) Wang Sang Koon, Martin W Lo,JE Marsden, Shane D Ross - "Dynamical
   Systems,the Three Body Problem and Space Mission Design", 2011
   (see chapter 4 (Sec 4.4 - Step 2) for description)
...
%}
function[XGuess] = InitialGuess(PointLoc,G_var)


Ax1 = G_var.Constants.Ax1;
Ax2 = G_var.Constants.Ax2;
mu = G_var.Constants.mu;


switch PointLoc %Location of Equilibruim (L1/L2/L3)
    case 1
        x_e = G_var.LagPts.L1;
        mu_bar = mu*abs(x_e(1) -1+mu)^-3 + (1-mu)*abs(x_e(1) + mu)^-3 ;
        nu = sqrt(-0.5*(mu_bar-2-sqrt(9*mu_bar^2-8*mu_bar)));
        Tau = -(nu^2 + 2*mu_bar+1)/(2*nu);
        nu_y0.one = -Ax1*nu*Tau;
        nu_y0.two = -Ax2*nu*Tau;
        XGuess(1,:) = [(x_e(1)-Ax1),0, 0, 0,nu_y0.one, 0];
        XGuess(2,:) = [(x_e(1)-Ax2),0, 0, 0,nu_y0.two, 0];

    case 2
        Ax1 = -Ax1;
        Ax2 = -Ax2;
        x_e = G_var.LagPts.L2;
        mu_bar = mu*abs(x_e(1) -1+mu)^-3 + (1-mu)*abs(x_e(1) + mu)^-3 ;
        nu = sqrt(-0.5*(mu_bar-2-sqrt(9*mu_bar^2-8*mu_bar)));
        Tau = -(nu^2 + 2*mu_bar+1)/(2*nu);
        nu_y0.one = -Ax1*nu*Tau;
        nu_y0.two = -Ax2*nu*Tau;
        XGuess(1,:) = [(x_e(1)-Ax1),0, 0, 0,nu_y0.one, 0];
        XGuess(2,:) = [(x_e(1)-Ax2),0, 0, 0,nu_y0.two, 0];


    case 3
        x_e = G_var.LagPts.L3;
        mu_bar = mu*abs(x_e(1) -1+mu)^-3 + (1-mu)*abs(x_e(1) + mu)^-3 ;
        nu = sqrt(-0.5*(mu_bar-2-sqrt(9*mu_bar^2-8*mu_bar)));
        Tau = -(nu^2 + 2*mu_bar+1)/(2*nu);
        nu_y0.one = -Ax1*nu*Tau;
        nu_y0.two = -Ax2*nu*Tau;
        XGuess(1,:) = [(x_e(1)-Ax1),0, 0, 0,nu_y0.one, 0];
        XGuess(2,:) = [(x_e(1)-Ax2),0, 0, 0,nu_y0.two, 0];

end