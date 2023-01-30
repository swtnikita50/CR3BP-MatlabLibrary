% Uses: findEigenStructureChange
%       continuatuion 

function [xLower, xUpper] = bifurcationInterval(familyPar)


%% Computation
stabilityIdx = familyPar.stabilityIdx;
initialIdx = stabilityIdx(1).center;
if initialIdx > 1
    isIdxGreaterThan1 = 1;
else
    isIdxGreaterThan1 = 0;
end

nBifurcationPts = 0;
aftOrbitIdx = [];
for i = 2:length(stabilityIdx)
    if stabilityIdx(i).center < 1 && isIdxGreaterThan1 
        nBifurcationPts = nBifurcationPts+1;
        aftOrbitIdx(end+1) = i;
        isIdxGreaterThan1 = 0;
    elseif stabilityIdx(i).center > 1 && ~isIdxGreaterThan1
        nBifurcationPts = nBifurcationPts+1;
        aftOrbitIdx(end+1) = i;
        isIdxGreaterThan1 = 1;
    end
end

for i = 1:nBifurcationPts
    xLower.IC(i,:) = familyPar.IC(aftOrbitIdx(i)-1,:);
    xLower.period(i) = familyPar.period(aftOrbitIdx(i)-1);
    xLower.jacobianConst(i) = familyPar.jacobianConst(aftOrbitIdx(i)-1);
    xLower.monodromy(:,:,i) = familyPar.monodromy(:,:,aftOrbitIdx(i)-1);
    xLower.eigens(i) = familyPar.eigens(aftOrbitIdx(i)-1);
    xLower.stabilityIdx(i) = familyPar.stabilityIdx(aftOrbitIdx(i)-1);
    
    xUpper.IC(i,:) = familyPar.IC(aftOrbitIdx(i),:);
    xUpper.period(i) = familyPar.period(aftOrbitIdx(i));
    xUpper.jacobianConst(i) = familyPar.jacobianConst(aftOrbitIdx(i));
    xUpper.monodromy(:,:,i) = familyPar.monodromy(:,:,aftOrbitIdx(i));
    xUpper.eigens(i) = familyPar.eigens(aftOrbitIdx(i));
    xUpper.stabilityIdx(i) = familyPar.stabilityIdx(aftOrbitIdx(i));
end