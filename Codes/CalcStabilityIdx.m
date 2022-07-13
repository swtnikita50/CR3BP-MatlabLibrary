% Created on: 12-07-22 (16:46)
function StabilityIdx = CalcStabilityIdx(Eigens)
for i = 1:length(Eigens.S_EigVal)
    StabilityIdx.Saddle(i) = 1/2*(abs(Eigens.S_EigVal(i)) + abs(Eigens.US_EigVal(i)));
end
for i = 2:2:length(Eigens.C_Val)
    StabilityIdx.C(i/2) = 1/2*(norm(Eigens.C_Val(i))+norm(1/Eigens.C_Val(i)));
end
end