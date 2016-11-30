function fit_grant_jun2016(subjid,exptype,model,runmax,nSamples)
%GRANT_FIT_JUN2016 Fit two subjects for the grant

% Optimization restarts
if nargin < 4 || isempty(runmax); runmax = 50; end
if nargin < 5; nSamples = []; end   % Use default

% runlist =  mod(procid-1,runmax) + 1;
runlist = 1;


table = fit_maximum_likelihood(subjid,exptype,model,runlist,runmax,nSamples);

save(sprintf('analysis/2_fitdata/fits/changelocalization_%s_model%d_subj%s.mat',exptype,model,subjid));

end