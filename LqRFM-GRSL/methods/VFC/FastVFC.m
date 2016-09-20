function VecFld = FastVFC(X, Y, conf)
% VFC  Vector Field Consensus
%   VECFLD = VFC_FAST(X, Y, CONF)
%   learns vector field from random samples with outliers. It is a fast
%   version of VFC.
%   
% Input:
%   X, Y: Original data.
%
%   conf: configuration for VFC. See VFC_init().
%
% Output:
%   VecFld: A structure type value which contains X, Y, beta, V, C, P,
%       VFCIndex. Where V = f(X), P is the posterior probability and 
%       VFCIndex is the indexes of inliers which found by VFC.
%
%   See also:: VFC_init(), VFC().

% Authors: Jiayi Ma (jyma2010@gmail.com)
% Date:    04/17/2012

fprintf('Start mismatch removal...\n');
[N, D]=size(Y); 

% Construct kernel matrix K
K=con_K(X,X,conf.beta);

% Default number of eigenvalue
NumEig = round(sqrt(N));
OPTS.issym=1;
OPTS.isreal=1;
OPTS.disp=0;
[Q,S]=eigs(K,NumEig,'lm',OPTS);
invS=spdiags(1./diag(S),0,NumEig,NumEig);

V = zeros(N,D); iter = 1;  tecr = 1; C = zeros(N,D); E = 1; 
sigma2 = sum(sum((Y-V).^2))/(N*D); gamma = conf.gamma;
%%
while (iter< conf.MaxIter) && (tecr > conf.ecr) && (sigma2 > 1e-8) 
    % E-step.
    E_old = E;
    [P, E] = get_P(Y, V, sigma2, gamma, conf.a);   

    QtC = Q'*C;
    E = E+conf.lambda/2*trace(QtC'*S*QtC);
    tecr = abs((E-E_old)/E);
    fprintf('iterate: %dth, gamma: %f, the energy change rate: %f, sigma2=%f\n', iter, gamma, tecr, sigma2);

    % M-step. Solve linear system for C.
    dP = spdiags(P,0,N,N);
    dPQ = dP*Q;
    F = dP*Y;
    C = 1/(conf.lambda*sigma2)*(F-dPQ*  ( (conf.lambda*sigma2*invS+Q'*dPQ) \ (Q'*F) ) );
    
    % update Y postions
    V = Q*(S*(Q'*C));

    Sp = sum(P);
    sigma2 = sum(P'*sum((Y-V).^2, 2))/(Sp*D);

    numcorr = length(find(P > conf.theta));
    gamma=numcorr/size(X, 1);
    if gamma > 0.95, gamma = 0.95; end
    if gamma < 0.05, gamma = 0.05; end
    iter=iter+1;
end

%%
% Fix the value of gamma, redo the EM process. 
V = zeros(N,D); iter = 1;  tecr = 1; C = zeros(N,D); E = 1; 
sigma2 = sum(sum((Y-V).^2))/(N*D);
while (iter < conf.MaxIter) && (tecr > conf.ecr) && (sigma2 > 1e-8) 
    % E-step.
    E_old = E;
    [P, E] = get_P(Y, V, sigma2, gamma, conf.a);   

    QtC = Q'*C;
    E = E+conf.lambda/2*trace(QtC'*S*QtC);
    tecr = abs((E-E_old)/E);
    fprintf('iterate: %dth, gamma: %f, the energy change rate: %f, sigma2=%f\n', iter, gamma, tecr, sigma2);

    % M-step. Solve linear system for C.
    dP = spdiags(P,0,N,N);
    dPQ = dP*Q;
    F = dP*Y;
    C = 1/(conf.lambda*sigma2)*(F-dPQ*  ( (conf.lambda*sigma2*invS+Q'*dPQ) \ (Q'*F) ) );
    
    % Update V and sigma^2
    V=Q*(S*(Q'*C));
    Sp=sum(P);
    sigma2=sum(P'*sum((Y-V).^2, 2))/(Sp*D);

    iter=iter+1;
end
%%
VecFld.X = X;
VecFld.Y = Y;
VecFld.beta = conf.beta;
VecFld.V=V;
VecFld.C=C;
VecFld.P = P;
VecFld.VFCIndex = find(P > conf.theta);

disp('Outlier removal succesfully completes.');


%%%%%%%%%%%%%%%%%%%%%%%%
function K=con_K(x,y,beta)
% CON_K constructs the kernel K, 
%   where K(i, j) = k(x, y) = exp(-beta*||x-y||^2).

[n, d]=size(x); [m, d]=size(y);

K=repmat(x,[1 1 m])-permute(repmat(y,[1 1 n]),[3 2 1]);
K=squeeze(sum(K.^2,2));
K=-beta * K;
K=exp(K);


%%%%%%%%%%%%%%%%%%%%%%%%
function [P, E]=get_P(Y,V, sigma2 ,gamma, a)
% GET_P estimates the posterior probability and part of the energy.

D = size(Y, 2);
temp1 = exp(-sum((Y-V).^2,2)/(2*sigma2));
temp2 = (2*pi*sigma2)^(D/2)*(1-gamma)/(gamma*a);
P=temp1./(temp1+temp2);
E=P'*sum((Y-V).^2,2)/(2*sigma2)+sum(P)*log(sigma2)*D/2;
