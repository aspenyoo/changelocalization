function table = fit_maximum_likelihood(subjid,exptype,model,runlist,runmax,nSamples)

if nargin < 5 || isempty(runmax); runmax = 50; end
if nargin < 6 || isempty(nSamples); nSamples = 1e4; end

subjid = upper(subjid);
dataname = sprintf('processeddata_ChangeLocalization_%s_subj%s',exptype,subjid);

% # samples for high-precision estimate
if numel(nSamples) > 1
    nSamplesFinal = nSamples(2);
else
    nSamplesFinal = 4e5;
end

% Load dataset
load(dataname);

% Set parameter bounds
jbar_bounds = [0.0067,35];  % Hard bounds for JBAR1 and JBAR2
jbar_pbounds = [0.0067,1];  % Plausible bounds for JBAR1 and JBAR2
tau_bounds = [0.5,3e3];     % Hard bounds for TAU
tau_pbounds = [35,3e3];     % Plausible bounds for TAU

LB = log([jbar_bounds(1),jbar_bounds(1),tau_bounds(1)]);
UB = log([jbar_bounds(2),jbar_bounds(2),tau_bounds(2)]);
PLB = log([jbar_pbounds(1),jbar_pbounds(1),tau_pbounds(1)]);
PUB = log([jbar_pbounds(2),jbar_pbounds(2),tau_pbounds(2)]);

% BPS options
options.UncertaintyHandling = 'on';

% Generate set of starting point with a Latin hypercube design
rng(0); % Same set for all
nvars = numel(PLB);
% x0_list = lhs(runmax,nvars,PLB,PUB,[],1e3);
x0_list = bsxfun(@plus,bsxfun(@times,rand(runmax,nvars),PUB-PUB),PLB);

table.xbest = NaN(runmax,3);
table.LLbest = NaN(runmax);

for iter = 1:numel(runlist)
    runlist(iter)
    
    % Fix random seed based on iteration (for reproducibility)
    rng(iter);
    
    x0 = x0_list(runlist(iter),:);
    [xbest,fval,exitflag,output] = ...
        bps(@(x) -AhyBCL_datalikeall(exp(x),data,model,nSamples(1)),x0,LB,UB,PLB,PUB,options);    
    xbest = exp(xbest);
        
    % Evaluate function with high precision
    LLbest = AhyBCL_datalikeall(xbest,data,model,nSamplesFinal);
    
    % Store results in table
    table.xbest(runlist(iter),:) = xbest;
    table.LLbest(runlist(iter)) = LLbest;
    
end


end