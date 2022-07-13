%{
...
Created on   26/2/2020 17:53
Modified on 2/7/22 20:21: added the initial guess part for halo orbit

This file calculates the Initial Guess for Differential Correction.

Inputs
------
1) PointLoc - Lagrange Points - Provide (1/2/3) as inputs
2) orbitType - Input the Orbit type - 'lyapunov'/'halo'
3) m - required for halo orbit, can on take values 1 or 3
1-Northern-Halo Orbits
2-Southern Halo Orbits

Outputs
--------
for all families
1) XGuess(1 and 2) - Guess for the initial point, a structure file containing
        * XGuess(1,:) - First guess
        * XGuess(2,:) - Second Guess 
Note: 
XGuess(i,:) - Guess in 3 Dimension(1 x 6) vector

Dependencies
------------
1) GlobalData - Takes "mu,Ax1,Ax2" for calculation for lyapunov orbit.
Ax1 and Ax2 for halo orbit is computed in this file

Note: Change Ax1 and Ax2 values in "Globaldata"

Reference
----------
1) Wang Sang Koon, Martin W Lo,JE Marsden, Shane D Ross - "Dynamical
   Systems,the Three Body Problem and Space Mission Design", 2011
   (see chapter 4 (Sec 4.4 - Step 2) for description)
...
%}
function[XGuess] = InitialGuess(PointLoc,G_var, orbitType, m)

if nargin <3
    orbitType = 'lyapunov';
end

Ax1 = G_var.Constants.Ax1;
Ax2 = G_var.Constants.Ax2;
mu = G_var.Constants.mu;

switch orbitType
    case 'lyapunov'
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
    case 'halo'

        delm = 2-m; %northern/southern halo
        switch PointLoc
            case 1
                dir = +1;
                gamma = G_var.LagPts.Gamma(1);
            case 2
                dir = -1;
                gamma = G_var.LagPts.Gamma(2);
        end
        XGuess = HaloLindstedtPoincare(gamma,dir,mu,delm);
end
