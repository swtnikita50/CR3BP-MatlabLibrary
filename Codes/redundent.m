clear clc close all
% System Parameters
N_x = 1000; % Number of grid points
length = 20; % Domain length
dx = length/N_x; 
dt = 0.001; 
Tf = 10; 

% Initialize the solution
x = linspace(-length/2, length/2, N_x+1); x(end) = [];
v = zeros(size(x));
v(x>0) = 1;

% MUSCL reconstruction
H = @(u) 0.5*u.^2; % Flux function
for t = 0:dt:Tf
v_prev = v;
vL = zeros(size(v));
vR = zeros(size(v));
vL(2:end-1) = v(1:end-2) + 0.25*minimod(v(2:end-1) - v(1:end-2), v(3:end) - v(2:end-1));
vR(2:end-1) = v(2:end-1) - 0.25*minimod(v(2:end-1) - v(1:end-2), v(3:end) - v(2:end-1));
fluxL = H(vL);
fluxR = H(vR);
for i = 2:N_x-1
v(i) = v_prev(i) - (dt/dx)*(fluxR(i-1) - fluxL(i));
end

% BC
v(1) = 0; 
v(N_x) = v(N_x-1); 
end

% Plotting
plot(x, v, 'k-.', 'LineWidth', 2);
xlabel('x');
ylabel('u(x,t)');
axis([-10, 10, 0, 1])
title('Solution of inviscid Burger''s equation: MUSCL reconstruction');

grid on

function k = minimod(a, b)
k = zeros(size(a));
idx = abs(a) < abs(b) & a./b > 0;
k(idx) = a(idx);
k(~idx) = b(~idx);
end

