% Incomplete
% Pseudo Arc-Length Continuation Method
% Created on: 11-07-22 (13:13)

function X_new = pseudoArcLengthCont(X_prev,dels,globalVar)
delX_prev = null(Dfmatrix3D(X_prev,globalVar.userInput.mu));
X_new = X_prev + dels*delX_prev;
X_new = diffCorrec(X_new,globalVar);