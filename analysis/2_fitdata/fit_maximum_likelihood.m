function fit_maximum_likelihood(subjid,model,exptype,runlist,runmax,nSamples)

if nargin < 5 || isempty(runmax); runmax = 50; end
if nargin < 6 || isempty(nSamples); nSamples = 1e4; end

% # samples for high-precision estimate
if numel(nSamples) > 1
    nSamplesFinal = nSamples(2);
else
    nSamplesFinal = 4e5;
end

% commented out bc not good for cluster
% try
%     load('ChangeLocalization_fits_Dec2016.mat')
% end
% 
% modstr = sprintf('model%d',model); % string used for saving

% Load dataset
load(sprintf('processeddata_ChangeLocalization_%s_subj%s.mat',exptype,subjid));

% Set parameter bounds
jbar_bounds = [0.0067,35];  % Hard bounds for JBAR1 and JBAR2
jbar_pbounds = [0.0067,1];  % Plausible bounds for JBAR1 and JBAR2
tau_bounds = [0.5,3e3];     % Hard bounds for TAU
tau_pbounds = [35,3e3];     % Plausible bounds for TAU
lapse_bounds = exp([0 1]);
lapse_pbounds = exp([0 0.2]);

LB = log([jbar_bounds(1),jbar_bounds(1),tau_bounds(1) lapse_bounds(1)]);
UB = log([jbar_bounds(2),jbar_bounds(2),tau_bounds(2) lapse_bounds(2)]);
PLB = log([jbar_pbounds(1),jbar_pbounds(1),tau_pbounds(1) lapse_pbounds(1)]);
PUB = log([jbar_pbounds(2),jbar_pbounds(2),tau_pbounds(2) lapse_pbounds(2)]);

% BPS options
options.UncertaintyHandling = 'on';

% Generate set of starting point with a Latin hypercube design
rng(0); % Same set for all
nvars = numel(PLB);
x0_list = lhs(runmax,nvars,PLB,PUB,[],1e3);

% stuff for file saving
filename = sprintf('ChangeLocalization_%s_model%d_subj%s.txt',exptype,model,subjid);
permission = 'a+';
formatSpec = repmat('%4.4f \t ',1,nvars+2);
formatSpec = [formatSpec(1:end-3) '\r\n'];

for iter = 1:numel(runlist)
    runlist(iter)
    
    % Fix random seed based on iteration (for reproducibility)
    rng(runlist(iter));
    
    x0 = x0_list(runlist(iter),:);
    [xbest,fval,exitflag,output] = ...
        bps(@(x) -AhyBCL_datalikeall(exp(x),data,model,nSamples(1)),x0,LB,UB,PLB,PUB,options);    
    xbest = exp(xbest);
        
    % Evaluate function with high precision
    LLbest = AhyBCL_datalikeall(xbest,data,model,nSamplesFinal);
    
    % save file
    fileID = fopen(filename,permission);
    A1 = [xbest LLbest runlist(iter)];
    fprintf(fileID, formatSpec, A1);
    fclose(fileID);
    
%     % commented out bc not good for cluster
%     fits.(modstr).(exptype).(subjid).xbest(runlist(iter),:) = xbest;
%     fits.(modstr).(exptype).(subjid).LLbest(runlist(iter)) = LLbest;
%     
%     save('ChangeLocalization_fits_Dec2016.mat','fits')
end


end