% Plots inavriant manifolds of orbits of different families:
% 1. Lyapunov
% 2. Halo
% More to come...

function plotInvManifold(globalVar)

% Extracting Initial Conditions for Plotting the Orbit
switch globalVar.userInput.orbit
    case 'lyapunov'
        orbPar = lyapunovOrbit(globalVar);
    case 'halo'
        orbPar = haloOrbit(globalVar);
end

manifoldType = globalVar.userInput.manifoldType;
f1 = globalVar.userInput.manifoldPlot;
switch manifoldType
    case 'stable'
        propogate = 'backward';
        color = 'g';
    case 'unstable'
        propogate = 'forward';
        color = 'r';
end

% Initial Condition for invarient Manifolds
invManifoldIC = orbitInvManifoldIC(globalVar,orbPar,40, manifoldType,1);
invManifoldIC2 = orbitInvManifoldIC(globalVar,orbPar,40, manifoldType,-1);

% Plotting
figure(f1)
for i = 1:length(invManifoldIC(:,1))
    [t,x] = integrate(globalVar,globalVar.functions.systemDynamics,invManifoldIC(i,1:6),[0 2*orbPar.period],propogate);
    plot3(x(:,1),x(:,2),x(:,3),color);hold on; grid on;
end
for i = 1:length(invManifoldIC2(:,1))
    [t,x] = integrate(globalVar,globalVar.functions.systemDynamics,invManifoldIC2(i,1:6),[0 2*orbPar.period],propogate);
    plot3(x(:,1),x(:,2),x(:,3),color);hold on; grid on;
end
scatter3(1-globalVar.userInput.mu,0,0,'p','filled','b');
%scatter3(-globalVar.userInput.mu,0,0,'p','filled','r');
xlabel('x (ND)')
ylabel('y (ND)')
zlabel('z (ND)')
title(['Invariant Manifold for ', globalVar.userInput.orbit ' Orbit']);
subtitle(['L', num2str(globalVar.userInput.lagrangePt) ' \mu = ', num2str(globalVar.userInput.mu) ', Jacobi Constant = ', num2str(globalVar.userInput.jacobianConst)]);
