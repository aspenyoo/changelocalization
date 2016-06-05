function LLcube(subjid, model, gridMat, nSamples, parts)
% 
% =========== INPUT VARIABLES ============
% SUBJID: subject id. string. 
% 
% MODEL: model number. (1) optimal. (2) fixed. 
%
% GRIDMAT: nParams x 3, such that each columns indicates the smallest
% value, the largest value, and the number of linear spaces between them,
% respectively, for each of the parameters. Parameters are probably Jbar1,
% Jbar2, and tau, so GRIDMAT should be 3 x 3. 
% 
% note: will be assumed that parameters are entered in log space, so will
% be exponentiating after linearly spacing gridMat(:,3) spaces between
% gridMat(:,1) and gridMat(:,2).
% 
% NSAMPLES: number of samples used to calculate log likelihood of one
% parameter combination. (default 1e4)
% 
% PARTS: vector with two elements. the first element is the number of parts
% the parameter combination matrix should be blocked into. The second
% element is which part you want to be evaulated. 
% 
% example: [3 2] would evaluate the 2nd third of the parameter combination
% meshgrid matrix.
% 
% April 28, 2016
% aspen.yoo@nyu.edu

if nargin < 4; nSamples = 1e4; end
if nargin < 5; parts = [1 1]; end

% load subject data
load(['data_ChangeLocalization_subj' subjid '_04282016.mat']);

% SET UP PARAMETER COMBINATIONS
% 1D parameters
Jbar1Vec = exp(linspace(gridMat(1,1),gridMat(1,2),gridMat(1,3)));         % Jbar for stimuli in the 1st interval
Jbar2Vec = exp(linspace(gridMat(2,1),gridMat(2,2),gridMat(2,3)));            % Jbar for stimuli in the 2nd interval
tauVec = exp(linspace(gridMat(3,1),gridMat(3,2),gridMat(3,3)));              % Scale parameter of the gamma distribution

% all possible parameter combinations
[X,Y,Z] = meshgrid(Jbar1Vec,Jbar2Vec,tauVec);
paramCombinationMat = [X(:) Y(:) Z(:)];
totalCombos = length(X(:));

% calculating which parameter combinations to calculate LL for
nParts = parts(1);
iPart = parts(2); 
paramValIdxs = (ceil((totalCombos/nParts)*(iPart-1))+1):ceil((totalCombos/nParts)*iPart);

% calculate LL
nCombo = length(paramValIdxs);
loglikeMat = nan(nCombo,2);
loglikeMat(:,1) = paramValIdxs';
for icombo = 1:nCombo;
    theta = paramCombinationMat(paramValIdxs(icombo),:); 
    
    % calculate log likelihood
    loglikeMat(icombo,2) = AhyBCL_datalikeall(theta,data,model,nSamples);    
end

% save loglikeMat
filename = sprintf('sep_LLcube_model%d_subjid%s_ipart%d_npart%d.mat',model,subjid,iPart,nParts);
save(filename,'loglikeMat');
