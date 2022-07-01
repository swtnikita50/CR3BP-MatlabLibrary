% ceated om 24/2/2020 17:48
% Variational Equation initial values (phi_0 and x_0)
function VarEq_InitVal = VarEq_Init(x0)
N=length(x0);
Phi_t0_t0 = reshape(eye(N),N^2,1);

VarEq_InitVal(N^2+1:N^2+N) = x0;

VarEq_InitVal(1:N^2) = Phi_t0_t0;

end