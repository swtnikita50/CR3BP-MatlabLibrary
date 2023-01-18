function plotBifurcationSln(G_var,OrbFamPar)

[xLower, xUpper] = findBifurcationPt(OrbFamPar);
[xMid] = ContinuationBifurcation(xLower, xUpper, G_var);

figure()
for i = 1:length(xMid.time)
[t,x] = Integrator(G_var,G_var.IntFunc.EOM,xMid.IC(i,:),[0 xMid.time(i)],'forward');
    plot(x(:,1),x(:,2));hold on; grid on;
end
xlabel('x (ND)')
ylabel('y (ND)')
title('Bifurcation Solutions')
axis equal
end