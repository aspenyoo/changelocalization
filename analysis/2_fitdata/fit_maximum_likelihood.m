function table = fit_maximum_likelihood(nid,model,runlist,runmax,nSamples)

if nargin < 4 || isempty(runmax); runmax = 50; end
if nargin < 5 || isempty(nSamples); nSamples = 1e4; end

dataname = 'ChangeLocalization_grant_data.mat';

% # samples for high-precision estimate
if numel(nSamples) > 1
    nSamplesFinal = nSamples(2);
else
    nSamplesFinal = 4e5;
end

% Load dataset
temp = load(dataname);
data = temp.data_all{nid};

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
x0_list = lhs(runmax,nvars,PLB,PUB,[],1e3);

maxsubjs = 2;
maxmodels = 2;

table.xbest = NaN(maxsubjs,maxmodels,runmax,3);
table.LLbest = NaN(maxsubjs,maxmodels,runmax);

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
    table.xbest(nid,model,runlist(iter),:) = xbest;
    table.LLbest(nid,model,runlist(iter)) = LLbest;
    
end


end