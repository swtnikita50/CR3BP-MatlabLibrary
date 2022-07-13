% Created 24/2/2020 18:33
% Variational equations that needs to be integrated
function [VarDot] = VarEqAndSTMDOT(t,VarInit,mu)
N = length(VarInit);
% G_var = GlobalData;
 %mu = G_var.Constants.mu;

if N == 20
    type = 'Planar';
elseif N == 42
    type = 'ThreeDim';
end

switch type
    case 'ThreeDim'

X =  VarInit(37:end);
Phi0 = reshape(VarInit(1:36),6,6);
Df = Dfmatrix3D(X,mu);
phidot = Df*Phi0;
PhiDot(1:36) = reshape(phidot,36,1);
PhiDot(37:42) = CRes3BP_EOM([],X,mu);

VarDot = PhiDot';

case 'Planar'

X =  VarInit(17:end);
Phi0 = reshape(VarInit(1:16),4,4);
Df = Dfmatrix3D(X,mu);
phidot = Df*Phi0;
PhiDot(1:16) = reshape(phidot,16,1);
PhiDot(17:20) = CRes3BP_EOM(t,X,mu);

VarDot = PhiDot';
end

 
end



