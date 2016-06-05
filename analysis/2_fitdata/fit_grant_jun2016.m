function fit_grant_jun2016(procid,runmax,nSamples)
%GRANT_FIT_JUN2016 Fit two subjects for the grant

% Optimization restarts
if nargin < 2 || isempty(runmax); runmax = 50; end
if nargin < 3; nSamples = []; end   % Use default

runlist =  mod(procid-1,runmax) + 1;  
modsubjid = floor((procid-1)/runmax) + 1;

switch modsubjid
    case 1; nid = 1; model = 1;
    case 2; nid = 2; model = 1;
    case 3; nid = 1; model = 2;
    case 4; nid = 2; model = 2;
end

table = fit_maximum_likelihood(nid,model,runlist,runmax,nSamples);

save(['changelocalization-fit-jun2016-' num2str(procid) '.mat']);

end