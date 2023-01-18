%{
...
Created on  02/07/22 22:10 - Nikita

This file does the computes the Halo orbit parameters.

Inputs
------
1) UserDat - Supplied UserData in main file.
2) globalVar   - Global data.
3) e       - specified jacobianContant

Outputs
--------
1) orbPar - A structure containing 
                - time      - full time period for Halo orbit.
                - IC        - Initial Conditions for Halo Orbit.
                - Energy    - Energy Value for Halo Orbit.
                - Monodromy - Monodromy matrix of Halo Orbit.
                - Eigens    - Eigen Values and Eigen Vectors for Halo Orbit.
                                    - Stable,Unstable and Center Eigen Values
                                    - Stable,Unstable and Center Eigen Vectors

Note that the size of eigen values and eigen vectors might change

Dependencies
------------
1) InitialGuess(PointLoc,globalVar)
2) diffCorrec(X_Guess,Plot,globalVar, isMaxIterReached)
                        

Reference for continuation 
--------------------------
1) Wang Sang Koon, Martin W Lo,JE Marsden, Shane D Ross - "Dynamical
   Systems,the Three Body Problem and Space Mission Design", 2011
...
%}
function [orbPar] = haloOrbit(globalVar)
% Extract the parameters
mu = globalVar.userInput.mu;
funVarEq = globalVar.functions.varEq_stmDot;
type = globalVar.userInput.type;
tol = globalVar.userInput.tolerance;
jacobianConstReq = globalVar.userInput.jacobianConst;

switch type
    case 'northern'
        m = 1;
        dir = 1;
    case 'southern'
        m = 3;
        dir = -1;
end
fprintf('\n===============================================\n')
for i = 1
    xGuess(i,:) = initGuess(globalVar);
    
    fprintf('Starting differential correction for orbit no.: %d\n',i)
    [t(i,1),x(i,:),~] = diffCorrec(xGuess(i,:),globalVar);
    jacobianConst(i,1) = jacobiValue3D(x(i,:),mu);
end

jacobianConstMax = jacobianConst(1,1);
isOrbitFound = 1;

if jacobianConst > jacobianConstMax
    fprintf(']\n===============================================\n')
    fprintf(['Error: The jacobian(energy) value enetered is more \n' ...
        'than the maximum limit for the Halo Orbit of this Libration Point\n'])
    fprintf('===============================================\n');
    isOrbitFound = 0;
else
end

i = 2;
delX0 = [0 0 1 0 0 0];
S = 0.01*dir;

while jacobianConst(i-1) > jacobianConstReq
    
    fprintf('Starting differential correction for orbit no.: %d\n',i)
    
%     if e < eMin
%         fprintf('\nError: Energy input is less than the %d.\n',eMin);
%         isOrbitFound = 0;
%         break;
%     end
    
    delta = delX0*S;
    xGuess = x(i-1,:) + delta;
    [~,~,~,isMaxIterReached] = diffCorrec(xGuess,globalVar);
    while isMaxIterReached
        delta = delta/2;
        xGuess = x(i-1,:) + delta;
        [~,~,~,isMaxIterReached] = diffCorrec(xGuess,globalVar);
    end
    [t(i,1),x(i,:),~,~] = diffCorrec(xGuess,globalVar);
    jacobianConst(i,1) = jacobiValue3D(x(i,:),mu);
    if abs(jacobianConst(i)-jacobianConst(i-1)) < tol
        disp('\n Error: The energy difference is negligible. \n');
        isOrbitFound = 0;
        break;
    end
    i = i+1;
    S = (x(i-1,3) - x(i-2,3));
end
xLowerAmp(1,:) = x(end-1,:); tLowerAmp(1,1) = t(end-1,1); jacobianConstLowerAmp(1,1) = jacobianConst(end-1,1);
xUpperAmp(1,:) = x(end,:); tUpperAmp(1,1) = t(end,1); jacobianConstUpperAmp(1,1) = jacobianConst(end,1);


switch isOrbitFound
    case 1
        fprintf('\n===============================================\n')
        fprintf('Found IC for orbit with energy lesser than and greater than %d\n',jacobianConstReq)
        fprintf('===============================================\n')


        tol = 1e-10;
        jacobioanConstNew = jacobianConstUpperAmp(1,1);
        while abs(jacobioanConstNew-jacobianConstReq)>tol
            %S = (xCorrecUpperAmp(1,3) - xCorrecLowerAmp(1,3))/2;
            %delta = delX0*S;
            delta = (xUpperAmp(1,:) - xLowerAmp(1,:))/2;
            xGuess = xLowerAmp(1,:) + delta;

            [t,x,~] = diffCorrec(xGuess,globalVar);
            jacobioanConstNew = jacobiValue3D(x,mu);
            % CHECK eNew with e and update in xcorrec and everywhere
            if jacobioanConstNew>jacobianConstReq
                xLowerAmp(1,:) = x;
                tLowerAmp(1,1) = t;
                jacobianConstLowerAmp(1,1) = jacobioanConstNew;
            else
                xUpperAmp(1,:) = x;
                tUpperAmp(1,1) = t;
                jacobianConstUpperAmp(1,1) = jacobioanConstNew;
            end
        end

        [~,monodromy,~,~] = stm_X(globalVar,x,funVarEq,t);
        [eigens.val.stable,eigens.val.unstable,eigens.val.center,eigens.val.p,eigens.vec.stable,...
            eigens.vec.unstable,eigens.vec.center,eigens.vec.p] = calcEigen(monodromy,1) ;
        stabilityIdx = calcStabilityIdx(eigens);
        
        orbPar.period      = t;
        orbPar.IC        = x;
        orbPar.jacobianConst    = jacobioanConstNew;
        orbPar.monodromy = monodromy;
        orbPar.eigens    = eigens;
        orbPar.stabilityIdx = stabilityIdx;
    case 0

        if e > jacobianConstMax
            fprintf('\n===============================================\n')
            fprintf('Error: The jacobian(energy) value enetered is greater than the maximum limit for the Libration Point\n')
            fprintf('===============================================\n');
        else
            fprintf('\nError: Energy input is less than the %d.\n',jacobianConst(end));
        end
        clear orbPar;
        orbPar = -1;
end

end