function [loglike,prmat] = AhyBCL_datalike1(theta,Nset1,Nset2,Cset,Rmat,model,Nsamples,ds)
%AHYBCL_DATALIKE1 Log likelihood for one condition (AHY BCL experiment).
%
%   LOGLIKE = AHYBCL_DATALIKE1(THETA,NSET1,NSET2,CSET,RMAT,MODEL) computes 
%   the log likelihood of change localization response matrix RMAT given
%   parameter vector THETA. NSET1 is the set size of the 1st interval;
%   NSET2 the set size of the 2nd interval; CSET the interval that contains
%   change (either 1 or 2). RMAT is a matrix with one row per stimulus type,
%   corresponding to the amount of change (the default experiment has 18 
%   stimuli types, with changes of 5:5:90 deg, due to assumed symmetry). 
%   Each column of RMAT contains the number of trials in which the observer 
%   responded: (1st) correctly; (2nd) incorrectly, but a stimulus within the 
%   same interval as where the change occurred; (3rd) incorrectly, and a
%   stimulus in the interval where change did not occur. If a whole column
%   does not apply, fill it with zeros.
%   MODEL is a model vector (for now, MODEL(1) determines the type of model,
%   MODEL(1) is 1 for Bayesian model, 2 for max rule).
%
%   The parameter vector THETA is:
%      THETA(1)  Jbar for stimuli in the 1st interval
%      THETA(2)  Jbar for stimuli in the 2nd interval
%      THETA(3)  Scale parameter tau of the gamma distribution
%      THETA(4)  Fixed lapse probability (default 0.005 if not specified)
%
%   LOGLIKE = AHYBCL_DATALIKE1(THETA,NSET1,NSET2,CSET,RMAT,MODEL,NSAMPLES) 
%   draws NSAMPLES samples for variable precision and noisy measurements 
%   (default value NSAMPLES=1e4).
%
%   LOGLIKE = AHYBCL_DATALIKE1(THETA,NSET1,NSET2,CSET,RMAT,MODEL,NSAMPLES,DS) 
%   uses a non-standard vector of stimuli types DS (expressed in deg).
%
%   [LOGLIKE,PRMAT] = AHYBCL_DATALIKE1(...) also returns the probability of 
%   response matrix PRMAT. PRMAT has the same structure as RMAT described
%   above. If calling AHYBCL_DATALIKE1 in predictive mode, RMAT can be left
%   empty, in which case LOGLIKE will be empty.

%   Author: Luigi Acerbi <luigi.acerbi@gmail.com>
%   Date: Apr/26/2016.

if nargin < 7 || isempty(Nsamples); Nsamples = 1e4; end
if nargin < 8 || isempty(ds); ds = (5:5:90)'; end

ds = ds(:)/90*pi;               % Convert from [-90,90] deg to [-pi,pi]

Nstim = size(ds,1);             % Number of stimuli types

Jbar1 = theta(1);               % Jbar for stimuli in the 1st interval
Jbar2 = theta(2);               % Jbar for stimuli in the 2nd interval
tau = theta(3);                 % Scale parameter of the gamma distribution
if numel(theta) < 4; theta(4) = 5e-3; end
lambda = theta(4);              % Fixed lapse probability

% Check parameters and input arguments
if lambda <= 0
    error('The lapse rate LAMBDA needs to be a positive scalar.');
end
if ~any(Cset == [1,2])
    error('CSET needs to be either 1 or 2 (signifying change in resp. 1st or 2nd interval).');
end
if (Cset == 1 && Nset1 == 0) || (Cset == 2 && Nset2 == 0)
    error('The change cannot be in an interval without stimuli.');
end
if ~isempty(Rmat) && any(size(Rmat) ~= [Nstim,3])
    error(['The response matrix RMAT needs to have ' num2str(Nstim) ' rows (one per stimulus type) and 3 columns.']);
end
if isempty(Rmat) && nargout < 2
    error('RMAT is empty but the number of output arguments is less than two.')
end

if Nset1 > 0
    % Is the change in the first set?
    change = (Cset == 1);    
    % Generate noisy measurements for stimuli in the 1st interval
    [x1,kappa1] = drawnoisy(Jbar1,tau,Nstim,Nsamples,Nset1,change,ds);    
    if size(x1,1) == 1; x1 = repmat(x1,[Nstim,1,1]); end
    if size(kappa1,1) == 1; kappa1 = repmat(kappa1,[Nstim,1,1]); end
else
    x1 = []; kappa1 = [];
end

if Nset2 > 0
    % Is the change in the second set?
    change = (Cset == 2);    
    % Generate noisy measurements for stimuli in the 2nd interval
    [x2,kappa2] = drawnoisy(Jbar2,tau,Nstim,Nsamples,Nset2,change,ds);
    if size(x2,1) == 1; x2 = repmat(x2,[Nstim,1,1]); end
    if size(kappa2,1) == 1; kappa2 = repmat(kappa2,[Nstim,1,1]); end
else
    x2 = []; kappa2 = [];
end

% Decision rule
switch model(1)
    case 1  % Bayesian model
        d1 = bsxfun(@plus,-bsxfun(@times,kappa1,cos(x1)),log(besseli(0,kappa1,1))+kappa1);
        d2 = bsxfun(@plus,-bsxfun(@times,kappa2,cos(x2)),log(besseli(0,kappa2,1))+kappa2);
        
    case 2  % Non-Bayesian max rule
        d1 = abs(x1);
        d2 = abs(x2);
end

d(:,:,1:Nset1) = d1;
d(:,:,Nset1 + (1:Nset2)) = d2;

[~,idx] = max(d,[],3);  % Get location with largest decision variable
    
% Compute responses depending on interval containing change
prmat = zeros(Nstim,3);
if Cset == 1
    prmat(:,1) = sum(idx == 1, 2);
    prmat(:,2) = sum(idx > 1 & idx <= Nset1, 2);
    if Nset2 > 0
        prmat(:,3) = sum(idx > Nset1, 2);
    end
else
    prmat(:,1) = sum(idx == (Nset1+1), 2);
    prmat(:,2) = sum(idx > (Nset1+1) & idx <= (Nset1+Nset2), 2);
    if Nset1 > 0
        prmat(:,3) = sum(idx <= Nset1, 2);
    end
end
prmat = prmat/Nsamples;

% Compute lapse probabilities
if Cset == 1
    lapsepdf = lambda*[1,Nset1-1,Nset2]/(Nset1+Nset2);
else
    lapsepdf = lambda*[1,Nset2-1,Nset1]/(Nset1+Nset2);
end

% Finalize log likelihood
prmat = bsxfun(@plus, bsxfun(@times,1-lambda,prmat), lapsepdf);

% Sum log likelihood over stimuli
if ~isempty(Rmat)
    temp = Rmat(:).*log(prmat(:));
    loglike = sum(temp(isfinite(temp)));
else
    loglike = [];
end

end

%--------------------------------------------------------------------------
function [x,kappa] = drawnoisy(Jbar,tau,Nstim,Nsamples,Nset,change,ds)
%DRAWNOISY Draw noisy measurements for a stimulus interval

% Is the change in this interval?
if change
    DS(:,1,1) = ds;                         % First stimulus changes
    DS(:,1,2:Nset) = zeros(Nstim,Nset-1);   % Other stimuli do not change
    DS = repmat(DS,[1,Nsamples,1]);
else
    DS = 0;                                 % No stimulus changes
end

alpha = Jbar/tau;             % Shape parameter for the gamma
J = gamrnd(alpha,tau,[Nstim,Nsamples,Nset]);
kappa = fisher2kappa(J);
x = qrandvm(DS,kappa,[],1e-3,1e2);

end