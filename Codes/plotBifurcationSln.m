function plotBifurcationSln(globalVar)

switch globalVar.userInput.orbit
    case 'lyapunov'
        [familyPar]              = lyapunovFamily(globalVar);
    case 'halo'
        [familyPar]              = haloFamily(globalVar);
end

%% Plotting
f1 = figure('Name','Unit Circle','NumberTitle','off');
viscircles([0,0],1), grid on; hold on;
for i = 1:globalVar.userInput.orbitCount
    scatter(real(familyPar.eigens(i).val.center(1)),imag(familyPar.eigens(i).val.center(1)),'ob','filled');
    scatter(real(familyPar.eigens(i).val.center(2)),imag(familyPar.eigens(i).val.center(2)),'ob','filled');
end
xlabel('Real')
ylabel('Imaginary')
title('Eigen Structure Associated with Center Subspace');

[xLower, xUpper] = bifurcationInterval(familyPar);
[xMid] = ContinuationBifurcation(xLower, xUpper, globalVar);

for i = 1:length(xMid.period)
    scatter(real(xMid.eigens(i).val.center(1)),imag(xMid.eigens(i).val.center(1)),'oy','filled');
    scatter(real(xMid.eigens(i).val.center(2)),imag(xMid.eigens(i).val.center(2)),'oy','filled');
end

f2 = figure('Name','Bifurcation Solutions','NumberTitle','off');
for i = 1:length(xMid.period)
[t,x] = integrate(globalVar,globalVar.functions.systemDynamics,xMid.IC(i,:),[0 xMid.period(i)],'forward');
    plot(x(:,1),x(:,2));hold on; grid on;
end
scatter3(globalVar.lagPts.pos(globalVar.userInput.lagrangePt,1),globalVar.lagPts.pos(globalVar.userInput.lagrangePt,2),0,'p','filled');

xlabel('x (ND)')
ylabel('y (ND)')
title('Bifurcation Solutions')
axis equal
end