function conf = VFC_init(conf)
% VFC  Vector Field Consensus
%   CONF = VFC_INIT(CONF) sets the default configuration for VFC.
%
%   gamma: Percentage of inliers in the samples. This is an inital value
%       for EM iteration, and it is not important. Default value is 0.9.
%
%   beta: Paramerter of Gaussian Kernel, k(x, y) = exp(-beta*||x-y||^2).
%       Default value is 0.1.
%
%   lambda: Represents the trade-off between the goodness of data fit 
%       and smoothness of the field. Default value is 3.
%
%   theta: If the posterior probability of a sample being an inlier is 
%       larger than theta, then it will be regarded as an inlier.
%       Default value is 0.75.
%
%   a: Paramerter of the uniform distribution. We assume that the outliers
%       obey a uniform distribution 1/a. Default Value is 10.
%
%   MaxIter: Maximum iterition times. Defualt value is 500.
%
%   ecr: The minimum limitation of the energy change rate in the iteration
%       process. Default value is 1e-5.
%
%   minP: The posterior probability Matrix P may be singular for matrix
%       inversion. We set the minimum value of P as minP. Default value is
%       1e-5.
%
%   method: Choose the method for outlier removal. There are three optional
%       methods: VFC, FastVFC, SparseVFC. Default value is VFC.
%
%   See also:: VFC(), FastVFC(), SparseVFC.

% Authors: Jiayi Ma (jyma2010@gmail.com)
% Date:    04/17/2012

if ~isfield(conf,'MaxIter'), conf.MaxIter = 500; end;
if ~isfield(conf,'gamma'), conf.gamma = 0.9; end;
if ~isfield(conf,'beta'), conf.beta = 0.1; end;
if ~isfield(conf,'lambda'), conf.lambda = 3; end;
if ~isfield(conf,'theta'), conf.theta = 0.75; end;
if ~isfield(conf,'a'), conf.a = 10; end;
if ~isfield(conf,'ecr'), conf.ecr = 1e-5; end;
if ~isfield(conf,'minP'), conf.minP = 1e-5; end;
if ~isfield(conf,'method'), conf.method = 'VFC'; end;