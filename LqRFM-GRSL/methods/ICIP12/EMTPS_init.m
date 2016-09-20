function conf =EMTPS_init(conf)

% Authors: Jiayi Ma (jyma2010@gmail.com)
% Date:    11/17/2012

if ~isfield(conf,'MaxIter'), conf.MaxIter = 200; end;
if ~isfield(conf,'gamma'), conf.gamma = 0.9; end;
if ~isfield(conf,'lambda'), conf.lambda = 500; end;
if ~isfield(conf,'theta'), conf.theta = 0.75; end;
if ~isfield(conf,'a'), conf.a = 5; end;
if ~isfield(conf,'ecr'), conf.ecr = 1e-5; end;
if ~isfield(conf,'minP'), conf.minP = 1e-5; end;