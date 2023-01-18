function plotOrbit(globalVar)

switch globalVar.userInput.orbit
    case 'lyapunov'
        [familyPar]              = lyapunovOrbit(globalVar);
    case 'halo'
        [familyPar]              = haloOrbit(globalVar);
end

figure()

[~,x] = integrate(globalVar,globalVar.functions.systemDynamics,familyPar.IC,[0 familyPar.period],'forward');
plot3(x(:,1),x(:,2),x(:,3));hold on; grid on;
scatter3(globalVar.lagPts.pos(globalVar.userInput.lagrangePt,1),globalVar.lagPts.pos(globalVar.userInput.lagrangePt,2),0,'p','filled');
xlabel('x (ND)')
ylabel('y (ND)')
zlabel('z (ND)')
title([globalVar.userInput.orbit ' orbit for \mu = ',num2str(globalVar.userInput.mu) ' and Jacobian = ',num2str(globalVar.userInput.jacobianConst)]);


end