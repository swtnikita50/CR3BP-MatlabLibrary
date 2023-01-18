% Created 24/2/2020 18:33
% Variational equations that needs to be integrated
function [varDot] = varEq_stmDot(t,varInit,mu)
N = length(varInit);

if N == 20
    type = 'Planar';
elseif N == 42
    type = 'ThreeDim';
end

switch type
    case 'ThreeDim'

        X =  varInit(37:end);
        Phi0 = reshape(varInit(1:36),6,6);
        Df = Dfmatrix3D(X,mu);
        phidot = Df*Phi0;
        PhiDot(1:36) = reshape(phidot,36,1);
        PhiDot(37:42) = CR3BP([],X,mu);

        varDot = PhiDot';

    case 'Planar'

        X =  varInit(17:end);
        Phi0 = reshape(varInit(1:16),4,4);
        Df = Dfmatrix3D(X,mu);
        phidot = Df*Phi0;
        PhiDot(1:16) = reshape(phidot,16,1);
        PhiDot(17:20) = CR3BP(t,X,mu);

        varDot = PhiDot';
end
 
end



