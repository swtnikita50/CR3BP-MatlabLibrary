% Incomplete
%Pseudo Arc-Length Continuation Method
% Created on: 11-07-22 (13:13)

function pseudoArcLengthCont(prevIC,G_var)
delXprev = null(Dfmatrix3D(prevIC,G_var.Constants.mu));
dels = 1/2;
Xnew = prevIC + dels*delXprev;
Xnew = pseduArcDiffCorrec(Xnew,dels,Plot,G_var,orbitType);