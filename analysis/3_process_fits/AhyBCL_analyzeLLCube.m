
function [post] = AhyBCL_analyzeLLCube(data)
%AHYBCL_ANALYZEDATACUBE Analyze log likelihood cube in DATA.

LLMat = data.LLMat;
Z = max(LLMat(:));
LLMat = LLMat - Z;

params = {'Jbar1','Jbar2','tau'};

post.params = params;

ip = 0;
idx = 1:numel(params);

for i = 1:numel(params)
    post.thetaVec{i} = log(data.([params{i} 'Vec']));
    post.dtheta(i) = diff(post.thetaVec{i}(1:2));
end

% Compute posterior
for i = 1:numel(params)    
    marg = setdiff(idx,i);
    temp = exp(LLMat);
    for j = numel(marg):-1:1
        temp = qtrapz(temp,marg(j))*post.dtheta(marg(j));
    end
    post.postpdf{i} = zeros(1,numel(temp));
    post.postpdf{i}(:) = temp/(qtrapz(temp)*post.dtheta(i));
end

% Compute marginal likelihood
post.mlike = log(sum(exp(LLMat(:)))) + Z + sum(log(post.dtheta));