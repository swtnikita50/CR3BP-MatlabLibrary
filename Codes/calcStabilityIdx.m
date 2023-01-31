% Computes stability index for orbits

function stabilityIdx = calcStabilityIdx(eigens)
    stabilityIdx.saddle = 1/2*(abs(eigens.val.stable) + abs(eigens.val.unstable));
    for i = 1:length(eigens.val.center)/2
        stabilityIdx.center(i) = 1/2*(eigens.val.center(2*i-1)+eigens.val.center(2*i));
    end
    if ~isempty(eigens.val.p)
        stabilityIdx.p = 1/2*(eigens.val.p(1)+eigens.val.p(2));
    else
        stabilityIdx.p = NaN;
    end
end
