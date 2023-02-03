% Plots family of orbits of different families:
% 1. Lyapunov
% 2. Halo
% More to come...

function plotFamily(globalVar)

% Extracting Initial Conditions for Plotting the Orbit
switch globalVar.userInput.orbit
    case 'lyapunov'
        [familyPar]              = lyapunovFamilyPseudoArcLengthCont(globalVar);
    case 'halo'
        [familyPar]              = haloFamilyPseudoArcLengthCont(globalVar);%haloFamily(globalVar);
end

% Plotting the Family
figure()
set(0,'DefaultAxesColorOrder',flipud(jet(length(familyPar.jacobianConst))));
for i = 1:length(familyPar.jacobianConst)
    [~,x] = integrate(globalVar,globalVar.functions.systemDynamics,familyPar.IC(i,:),[0 familyPar.period(i)],'forward');
    plot3(x(:,1),x(:,2),x(:,3));hold on; grid on;
end
scatter3(globalVar.lagPts.pos(globalVar.userInput.lagrangePt,1),globalVar.lagPts.pos(globalVar.userInput.lagrangePt,2),0,'p','filled');
caxis([familyPar.jacobianConst(end), familyPar.jacobianConst(1)]);
colorbar
xlabel('x (ND)')
ylabel('y (ND)')
zlabel('z (ND)')
title([globalVar.userInput.orbit ' orbit family for \mu = ',num2str(globalVar.userInput.mu)]);

