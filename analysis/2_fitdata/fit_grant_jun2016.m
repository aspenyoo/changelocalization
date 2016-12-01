function fit_grant_jun2016(procid,runmax,nSamples)
%GRANT_FIT_JUN2016 Fit two subjects for the grant

% Optimization restarts
if nargin < 2 || isempty(runmax); runmax = 50; end
if nargin < 3; nSamples = []; end   % Use default

subjids = {'ALM','DR','EN','MR'};
nModels = 2;                % number of models

runlist = mod(procid-1,runmax)+1; % which of 50 starting values to use

modsubjid = floor((procid-1)/runmax) + 1; % subject and model id

subjnum = ceil(modsubjid/nModels);
subjid = subjids{subjnum};
model = mod(modsubjid-1,nModels) + 1;

fit_maximum_likelihood(subjid,model,'Contrast',runlist,runmax,nSamples);
fit_maximum_likelihood(subjid,model,'Delay',runlist,runmax,nSamples);

end