function [bfp, LLVec, completedruns] = fit_parameters(data,model,runlist,runmax,nSamples)
%FIT_PARAMETERS was made by aspen oct 31 2018 to try to investigate
%different aspects of model. same as fit_maximum_likelihood but outputs
%varaibles instead of saving them to a file

if nargin < 4 || isempty(runmax); runmax = 50; end
if nargin < 5 || isempty(nSamples); nSamples = 1e4; end

% # samples for high-precision estimate
if numel(nSamples) > 1
    nSamplesFinal = nSamples(2);
else
    nSamplesFinal = 4e5;
end

% % Load dataset
% filepath = '/home/ay963/changelocalization/experiment_data/output_mat/';
% % filepath = 'experiment_data/output_mat/';
% load(sprintf('%sprocesseddata_ChangeLocalization_%s_subj%s.mat',filepath,exptype,subjid));

% Set parameter bounds
jbar_bounds = [0.0067 35];  % Hard bounds for JBAR1 and JBAR2
jbar_pbounds = [0.5 10];  % Plausible bounds for JBAR1 and JBAR2
tau_bounds = [1e-3 3e3];     % Hard bounds for TAU
tau_pbounds = [0.5 30];     % Plausible bounds for TAU
lapse_bounds = [1e-4 1];
lapse_pbounds = [1e-4 0.2];

LB = [log([jbar_bounds(1),jbar_bounds(1),tau_bounds(1)]) lapse_bounds(1)];
UB = [log([jbar_bounds(2),jbar_bounds(2),tau_bounds(2)]) lapse_bounds(2)];
PLB = [log([jbar_pbounds(1),jbar_pbounds(1),tau_pbounds(1)]) lapse_pbounds(1)];
PUB = [log([jbar_pbounds(2),jbar_pbounds(2),tau_pbounds(2)])  lapse_pbounds(2)];

% BPS options
options.UncertaintyHandling = 'on';

% Generate set of starting point with a Latin hypercube design
rng(0); % Same set for all
nvars = numel(PLB);
x0_list = lhs(runmax,nvars,PLB,PUB,[],1e3);

% % stuff for file saving
% filepath = '/home/ay963/changelocalization/analysis/2_fitdata/fits/';
% filename = sprintf('%sChangeLocalization_%s_model%d_subj%s.txt',filepath,exptype,model,subjid);
% permission = 'a+';
% formatSpec = repmat('%4.4f \t ',1,nvars+2);
% formatSpec = [formatSpec(1:end-3) '\r\n'];

[bfp, LLVec, completedruns] = deal([]);
for iter = 1:numel(runlist)
    runlist(iter)
    
    % Fix random seed based on iteration (for reproducibility)
    rng(runlist(iter));
    
    x0 = x0_list(runlist(iter),:);
    [xbest,LL,~,~] = ...
        bads(@(x) -AhyBCL_datalikeall([exp(x(1:3)) x(4)],data,model,nSamples(1)),x0,LB,UB,PLB,PUB,[],options);
%     xbest = [exp(xbest(1:3)) xbest(4)];
    
    % Evaluate function with high precision
%     LLbest = AhyBCL_datalikeall(xbest,data,model,nSamplesFinal);
    
    bfp = [bfp; xbest];
    LLVec = [LLVec; LL];
    completedruns = [completedruns; runlist(iter)];
    
end

end