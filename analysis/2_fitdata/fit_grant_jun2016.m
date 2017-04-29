function fit_grant_jun2016(subjid,model,runlist,runmax,nSamples)
%GRANT_FIT_JUN2016 Fit two subjects for the grant

% Optimization restarts
if nargin < 2 || isempty(runmax); runmax = 50; end
if nargin < 3; nSamples = []; end   % Use default

fit_maximum_likelihood(subjid,model,'Contrast',runlist,runmax,nSamples);
fit_maximum_likelihood(subjid,model,'Delay',runlist,runmax,nSamples);

end