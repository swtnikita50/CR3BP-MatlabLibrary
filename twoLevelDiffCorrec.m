% Incomplete
% Two Level Differential Correction
% Created on: 10-07-22 23:45


function twoLevelDiffCorrec(initPatchPts,timeStamps, G_var)
tol = 1e-3;
fun = G_var.IntFunc.VarEqAndSTMdot;
nSegments = length(initPatchPts(:,1));

% Level 1: Position Continuity
for i = 1:nSegments
    err = 1;
    X0 = initPatchPts(i,:);
    X1 = initPatchPts(i+1,:);
    while err> tol
        [t,PHIt1,x,x1,PHI] = StateTransAndX(G_var,X0,fun,timeStamps(i+1),timeStamps(i));
        L = [PHIt1(1,4) PHIt1(1,5) PHIt1(1,6) x1(4);...
            PHIt1(2,4) PHIt1(2,5) PHIt1(2,6) x1(5);...
            PHIt1(3,4) PHIt1(3,5) PHIt1(3,6) x1(6)];
        b = [X1(1)-x1(1); X1(2)-x1(2); X1(2)-x1(2)];
        u = L'*(L*L')^(-1)*b;
        err = norm(u);
        X0(4) = X0(4)+u(1);
        X0(5) = X0(5)+u(2);
        X0(6) = X0(6)+u(3);
        timeStamps(i+1) = timeStamps(i+1)+u(4);
    end
end

% Level 2: Velocity Continuity
for i = 2:nSegments
    X0 = initPatchPts(i,:);
    Xp = initPatchPts(i+1,:);
    Xm = initPatchPts(i-1,:);
    [t,PHIt0p,x,x0p,PHI] = StateTransAndX(G_var,Xp,fun,timeStamps(i+1),timeStamps(i),'backward');
    [t,PHIt0m,x,x0m,PHI] = StateTransAndX(G_var,Xm,fun,timeStamps(i),timeStamps(i-1),'forward');
    A0p = PHIt0p(1:3,1:3); B0p = PHIt0p(1:3,4:6); C0p = PHIt0p(4:6,1:3); D0p = PHIt0p(4:6,4:6);
    A0m = PHIt0m(1:3,1:3); B0m = PHIt0m(1:3,4:6); C0m = PHIt0m(4:6,1:3); D0m = PHIt0m(4:6,4:6);
    Mm  = D0m*(B0m)^(-1)*A0m-C0m;
    %get acceleration
    Mtm = a0m - D0m*(B0m)^(-1)*v0m;
    M0  = D0p*B0p^(-1)-D0m*B0m^(-1);
    Mt0 = D0m*B0m^(-1)*v0m - D0p*B0p^(-1)*v0p + a0p - a0m;
    Mp  = C0p - D0p*B0p^(-1)*A0p;
    Mtp = D0p*B0p^(-1)*v0p-a0p;
    M(:,:,i) = [Mm, Mtm, M0, Mt0, Mp, Mtp];
    Mfinal(3*(i-2+1):3*(i-1),:) = [zeros(3,4*(i-2)),M(:,:,i),zeros(3,4*(n-(i+1)))];
    %get velocity differences
    delV(i,:) = v0p-v0m;
end

delr = -Mfinal'*(Mfinal*Mfinal')^(-1)*delV;

%fix the positions and time

