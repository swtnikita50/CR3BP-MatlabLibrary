% Copyright Nikita December 2022
% Uses: findEigenStructureChange
%       continuatuion 

function [xLower, xUpper] = findBifurcationPt(OrbFamPar)

stabilityIdx = OrbFamPar.StabilityIdx;
initialIdx = stabilityIdx(1).C;
if initialIdx > 1
    isIdxGreaterThan1 = 1;
else
    isIdxGreaterThan1 = 0;
end

nBifurcationPts = 0;
aftOrbitIdx = [];
for i = 2:length(stabilityIdx)
    if stabilityIdx(i).C < 1 && isIdxGreaterThan1 
        nBifurcationPts = nBifurcationPts+1;
        aftOrbitIdx(end+1) = i;
        isIdxGreaterThan1 = 0;
    elseif stabilityIdx(i).C > 1 && ~isIdxGreaterThan1
        nBifurcationPts = nBifurcationPts+1;
        aftOrbitIdx(end+1) = i;
        isIdxGreaterThan1 = 1;
    end
end

for i = 1:nBifurcationPts
    xLower.IC(i,:) = OrbFamPar.IC(aftOrbitIdx(i)-1,:);
    xLower.time(i) = OrbFamPar.time(aftOrbitIdx(i)-1);
    xLower.Energy(i) = OrbFamPar.Energy(aftOrbitIdx(i)-1);
    xLower.Monodromy(:,:,i) = OrbFamPar.Monodromy(:,:,aftOrbitIdx(i)-1);
    xLower.Eigens(i) = OrbFamPar.Eigens(aftOrbitIdx(i)-1);
    xLower.StabilityIdx(i) = OrbFamPar.StabilityIdx(aftOrbitIdx(i)-1);
    
    xUpper.IC(i,:) = OrbFamPar.IC(aftOrbitIdx(i),:);
    xUpper.time(i) = OrbFamPar.time(aftOrbitIdx(i));
    xUpper.Energy(i) = OrbFamPar.Energy(aftOrbitIdx(i));
    xUpper.Monodromy(:,:,i) = OrbFamPar.Monodromy(:,:,aftOrbitIdx(i));
    xUpper.Eigens(i) = OrbFamPar.Eigens(aftOrbitIdx(i));
    xUpper.StabilityIdx(i) = OrbFamPar.StabilityIdx(aftOrbitIdx(i));
end