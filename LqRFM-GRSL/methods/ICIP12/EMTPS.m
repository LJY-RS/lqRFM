function Transform = EMTPS(X, Y, gamma, lambda, theta, a, MaxIter, ecr, minP)

% Authors: Jiayi Ma (jyma2010@gmail.com)
% Date:    11/17/2012

fprintf('Starting mismatch removal:\n');
[N, D]=size(X); 

% Construct kernel matrix K
K = tps_gen_K(X,X);

% Initialization
V=X; iter=1;  tecr=1; W=zeros(N,D); H=zeros(D,D); E=1; 
sigma2=sum(sum((Y-X).^2))/(N*D);
%%
newE = [];%
% newE2 = [];%
% while (iter<MaxIter) && (tecr > ecr) && (sigma2 > 1e-8) 
while (iter<MaxIter) && (tecr > ecr) && (sigma2 > 1e-8) 
    % E-step. 
    E_old=E;
    [P, E]=get_P(Y,V, sigma2 ,gamma, a);%%%%%%%%%%%%%%%%%%%%%%%%%%%% Notice: the energy function has omitted the terms containing gamma

    E=E+lambda/2*trace(W'*K*W);
    newE = [newE,E];%   
    tecr=abs((E-E_old)/E);
%     newE2 = [newE2,tecr];%
    fprintf('iterate: %dth, gamma: %f, the energy change rate: %f, sigma2=%f\n', iter, gamma, tecr, sigma2);

    % M-step. Update W, H.
    XX = X.*repmat(sqrt(P), [1, D]);
    YY = Y.*repmat(sqrt(P), [1, D]);
    [q1,q2,R] = tps_gen_qr(XX);
    P = max(P, minP);
    A = q2'.*repmat(sqrt(P)', [N-D, 1])*K*q2;
    B = q2'*K*q2;
    W = q2*((A'*A+lambda*sigma2*(B+eye(N-D)*1e-5))\(A'*q2'*YY));
    H = R\(q1'*(YY-K.*repmat(sqrt(P), [1, N])*W));
    
%     XX = X.*repmat(P, [1, D]);
%     YY = Y.*repmat(P, [1, D]);
%     [q1,q2,R] = tps_gen_qr(XX);
%     W = q2*((q2'*K*q2+lambda*sigma2*eye(N-D))\(q2'*YY));
%     H = R\(q1'*(YY-K*W));

%     LAM = diag(sqrt(P)); tmp = LAM*q2;
%     W = tmp*(inv(tmp'*K*tmp+lambda*sigma2*eye(size(tmp,2))))*tmp'*Y;
%     H = inv(R) * q1' * LAM * (Y - K*W);
    
    % Update V and sigma^2
    V=X*H + K*W;
    VD = V(:,D);    VD(abs(VD)<1e-10) = 1e-10;
    V = V./repmat(VD, [1, D]);
    Sp=sum(P);
    sigma2=sum(P'*sum((Y-V).^2, 2))/(Sp*D);

    % Update gamma
    numcorr = length(find(P > theta));
    gamma=numcorr/size(X, 1);
%     gamma = Sp/size(X, 1);
    if gamma > 0.95, gamma = 0.95; end
    if gamma < 0.05, gamma = 0.05; end
    
    iter=iter+1;
end
%%
% Fix the value of gamma, redo the EM process. 
% V=X; iter=1;  tecr=1; C=zeros(M,D); E=1; 
% sigma2=sum(sum((Y-X).^2))/(N*D);
% while (iter<MaxIter) && (tecr > ecr) && (sigma2 > 1e-8) 
%     % E-step.
%     E_old=E;
%     [P, E]=get_P(Y,V, sigma2 ,gamma, a);   
% 
%     E=E+lambda/2*trace(C'*K*C);
%     tecr=abs((E-E_old)/E);
% %     fprintf('iterate: %dth, gamma: %f, the energy change rate: %f, sigma2=%f\n', iter, gamma, tecr, sigma2);
% 
%     % M-step. Solve linear system for C.
%     P = max(P, minP);
%     C = (U'.*repmat(P', [M, 1])*U+lambda*sigma2*K)\(U'.*repmat(P', [M, 1])*Y);
% 
%     % Update V and sigma^2
%     V=U*C;
%     Sp=sum(P);
%     sigma2=sum(P'*sum((Y-V).^2, 2))/(Sp*D);
% 
%     iter=iter+1;
% end

%%
Transform.X = X(:,1:D-1);
Transform.Y = Y(:,1:D-1);
Transform.V=V(:,1:D-1);
Transform.H=H;
Transform.W=W;
Transform.P = P;
Transform.E = newE;
Transform.Index = find(P > theta);

disp('Removing outliers succesfully completed.');


%%%%%%%%%%%%%%%%%%%%%%%%
function [K] = tps_gen_K(x,z)

% Format:
[n, M] = size (x); 
[m, M] = size (z);
dim    = M  - 1;

% calc. the K matrix.
% 2D: K = r^2 * log r
% 3D: K = -r
K= zeros (n,m);

for it_dim=1:dim
  tmp = x(:,it_dim) * ones(1,m) - ones(n,1) * z(:,it_dim)';
  tmp = tmp .* tmp;
  K = K + tmp;
end;
  
if dim == 2
  mask = K < 1e-10; % to avoid singularity.
  K = 0.5 * K .* log(K + mask) .* (K>1e-10);
else
  K = - sqrt(K);
end;

%%%%%%%%%%%%%%%%%%%%%%%%
function [q1,q2,R] = tps_gen_qr(x)

[n,M] = size (x);

[q,r]   = qr(x);
q1      = q(:, 1:M);
q2      = q(:, M+1:n);
R       = r(1:M,1:M);

%%%%%%%%%%%%%%%%%%%%%%%%
function [P, E]=get_P(Y,V, sigma2 ,gamma, a)
% GET_P estimates the posterior probability and part of the energy.

D = size(Y, 2);
temp1 = exp(-sum((Y-V).^2,2)/(2*sigma2));
temp2 = (2*pi*sigma2)^(D/2)*(1-gamma)/(gamma*a);
P=temp1./(temp1+temp2);
E=P'*sum((Y-V).^2,2)/(2*sigma2)+sum(P)*log(sigma2)*D/2;
