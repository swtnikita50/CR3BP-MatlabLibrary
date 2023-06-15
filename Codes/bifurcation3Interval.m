% Uses: findEigenStructureChange
%       continuatuion 

function [xLower, xUpper, nBifurcationPts] = bifurcation3Interval(familyPar)


%% Computation
stabilityIdx = familyPar.stabilityIdx;
initialIdx = stabilityIdx(1).saddle;
if isempty(initialIdx) 
    isInitSaddleEmpty = 1;
else
    isInitSaddleEmpty = 0;
end

nBifurcationPts = 0;
aftOrbitIdx = [];
for i = 2:length(stabilityIdx)
    if ~isempty(stabilityIdx(i).saddle) && isInitSaddleEmpty 
        nBifurcationPts = nBifurcationPts+1;
        aftOrbitIdx(end+1) = i;
        isInitSaddleEmpty = 0;
    elseif isempty(stabilityIdx(i).saddle) && ~isInitSaddleEmpty
        nBifurcationPts = nBifurcationPts+1;
        aftOrbitIdx(end+1) = i;
        isInitSaddleEmpty = 1;
    end
end

if nBifurcationPts~=0
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
else 
    xUpper = NaN;
    xLower = NaN;
end