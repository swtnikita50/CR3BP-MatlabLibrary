function X = HaloLindstedtPoincare(gamma,dir,mu, delm,Ax)

c2 = 1/gamma^3*((dir)^2*mu+(-1)^2*(1-mu)*gamma^(2+1)/(1-dir*gamma)^(2+1));
c3 = 1/gamma^3*((dir)^3*mu+(-1)^3*(1-mu)*gamma^(3+1)/(1-dir*gamma)^(3+1));
c4 = 1/gamma^3*((dir)^4*mu+(-1)^4*(1-mu)*gamma^(4+1)/(1-dir*gamma)^(4+1));
lambda = sqrt((c2+sqrt(9*c2^2-8*c2))/2);
wp = sqrt(lambda^2+1-c2);
lambda = wp;
wv = sqrt(c2);
kappa = 2*lambda/(lambda^2+1-c2);
k = kappa;

d1 = 3*lambda^2/kappa*(kappa*(6*lambda^2-1)-2*lambda);
d2 = 8*lambda^2/kappa*(kappa*(11*lambda^2-1)-2*lambda);
a21 = 3*c3*(kappa^2-2)/(4*(1+2*c2));
a22 = 3*c3/(4*(1+2*c2));
a23 = -3*c3*lambda*(3*kappa^3*lambda-6*kappa*(kappa-lambda)+4)/(4*kappa*d1);
a24 = -3*c3*lambda*(2+3*kappa*lambda)/(4*kappa*d1);
b21 = -3*c3*lambda*(3*kappa*lambda-4)/(2*d1);
b22 = 3*c3*lambda/d1;
d21 = -c3/(2*lambda^2);
a31 = -9*lambda*(4*c3*(kappa*a23-b21)+kappa*c4*(4+kappa^2))/(4*d2)...
    + (9*lambda^2+1-c2)*(3*c3*(2*a23-kappa*b21)+c4*(2+3*kappa^2))/(2*d2);
a32 = -9*lambda*(4*c3*(kappa*a24-b22)+kappa*c4)/(4*d2)...
    - 3*(9*lambda^2+1-c2)*(c3*(kappa*b22+d21-2*a24)-c4)/(2*d2);
b31 = (3.0 * lambda * (3.0 * c3 * (k * b21 - 2.0 * a23) - c4 * (2.0 + 3.0 * k^2)) + ...
    (9.0 * lambda^2 + 1.0 + 2.0 * c2) * (12.0 * c3 * (k * a23 - b21) + ...
    3.0 * k * c4 * (4.0 + k^2)) / 8.0) / d2;
b32 = 9*lambda/d2*(c3*(kappa*b22+d21-2*a24)-c4)...
    +3*(9*lambda^2+1+2*c2)*(4*c3*(kappa*a24-b22)+kappa*c4)/(8*d2);
d31 = 3*(4*c3*a24+c4)/(64*lambda^2);
d32 = 3*(4*c3*(a23-d21)+c4*(4+kappa^2))/(64*lambda^2);
s1 = (2*lambda*(lambda*(1+kappa^2)-2*kappa))^(-1)*(3/2*c3*(2*a21*(kappa^2-2)...
    -a23*(kappa^2+2)-2*kappa*b21) - 3/8*c4*(3*kappa^4-8*kappa^2+8));
d3 = (2*lambda*(lambda*(1+kappa^2)-2*kappa));
s2 = ((3.0 / 2.0) * c3 * (2.0 * a22 * (k^2 - 2.0) + ...
    a24 * (k^2 + 2.0) + 2.0 * k * b22 + 5.0 * d21) + ...
    (3.0 / 8.0) * c4 * (12.0 - k^2)) / d3;

a2 = (3.0 / 2.0) * c3 * (a24 - 2.0 * a22) + (9.0 / 8.0) * c4;

l1 = -3/2*c3*(2*a21+a23+5*d21)-3/8*c4*(12-kappa^2)+2*lambda^2*s1;
l2 = 2.0 * s2 * lambda^2 + a2;
delta = wp^2-c2;

Ax = sqrt(abs(delta/l1));
Az = sqrt((-delta-l1*Ax^2)/l2);
nu2 = s1*Ax^2+s2*Az^2;
nu = 1+nu2;

x0 = a21*Ax^2+a22*Az^2-Ax+a23*Ax^2-a24*Az^2+a31*Ax^3-a32*Ax*Az^2;
z0 = delm*Az+delm*d21*Ax*Az*(1-3)+delm*(d32*Az*Ax^2-d31*Az^3);
ydot0 = kappa*Ax+2*(b21*Ax^2-b22*Az^2)+3*(b31*Ax^3-b32*Ax*Az^2);

x0 = gamma*x0-dir*gamma+1-mu;
z0 = gamma*z0;
ydot0 = gamma*wp*nu*ydot0;

X = [x0,0, z0, 0,ydot0, 0];
end