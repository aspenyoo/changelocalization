function [best_contrast] = calculate_bestcontrast(subjid)
% calculates the mu across conditions and different psychometric functions.
% returns "best contrast" in percents (median)

% get data from subject
load(['Training_ChangeLocalization_' subjid '.mat'])

% get mu's for each of the three psychometric functions
nPsychos = 3; % three psychometric functions
condstrVec = {'cond1','cond2','cond3','cond4'};
contrastVec = psy.(condstrVec{1}).mu;

nCond = length(condstrVec);
nMus = length(psy.(condstrVec{1}).mu);

muVec_condition = nan(nCond,nMus);
for icond = 1:nCond
    condstr = condstrVec{icond};
    
    all_mu.(condstr) = nan(nPsychos,nMus);
    for ipsycho = 1:nPsychos
        all_mu.(condstr)(ipsycho,:) = sum(sum(psy.(condstr).post{ipsycho},2),3); 
    end
    
    muVec_condition(icond,:) = bsxfunandsum(@times,psy.(condstr).psychopost',all_mu.(condstr));
    
end

mupdf = mean(muVec_condition);
mucdf = cumsum(mupdf);

% get which values are closes to 50 and interpolate between them
lowidx = find(mucdf < 0.5,1,'last'); % index of one below 0.5
highidx = find(mucdf > 0.5,1,'first'); % index of one above 0.5
slope = diff(contrastVec(lowidx:highidx))/diff(mucdf(lowidx:highidx));
best_contrast = contrastVec(lowidx) + slope*(0.5 - mucdf(lowidx)); % linear interpolation
best_contrast = exp(best_contrast)/100; % in proportion contrast

% plot(exp(contrastVec),muVec_condition)