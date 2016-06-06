function [loglike,prmat,X] = AhyBCL_datalikeall(theta,X,model,Nsamples)
%AHYBCL_DATALIKEALL Log likelihood for all conditions (AHY BCL experiment).
%
%   LOGLIKE = AHYBCL_DATALIKEALL(THETA,X,MODEL) computes the log likelihood 
%   of change localization dataset X given parameter vector THETA and model
%   MODEL. X is a cell array of data structures that contain fields Nset1,
%   Nset2, Cset, and Rmat for each experimental condition (see 
%   AHYBCL_DATALIKE1).
%
%   LOGLIKE = AHYBCL_DATALIKEALL(...,NSAMPLES) draws NSAMPLES noisy 
%   measurements per trial, stimulus and condition (if not specified
%   uses default value as per AHYBCL_DATALIKE1).
%
%   [LOGLIKE,PRMAT] = AHYBCL_DATALIKEALL(...) returns a cell array of 
%   probability of response matrices, one per condition (as per X).
%
%   [LOGLIKE,PRMAT,X] = AHYBCL_DATALIKEALL(THETA,X,MODEL) returns a cell 
%   array of data structures X filled with simulated response data. The
%   input X can be empty.
%
%   See AHYBCL_DATALIKE1.

if nargin < 4; Nsamples = []; end   % Use default

if isempty(X)   % Fill X with default experimental structure
    fakerun = 1;
    X{1}.Nset1 = 0; X{1}.Nset2 = 4; X{1}.Cset = 2;
    X{2}.Nset1 = 2; X{2}.Nset2 = 2; X{2}.Cset = 1;
    X{3}.Nset1 = 2; X{3}.Nset2 = 2; X{3}.Cset = 2;
    X{4}.Nset1 = 4; X{4}.Nset2 = 0; X{4}.Cset = 1;
    for iCnd = 1:numel(X); X{iCnd}.Rmat = []; end
else
    fakerun = 0;
end

ll = zeros(1,numel(X));

% Get log likelihood and probability response matrix for each condition
for iCnd = 1:numel(X)
    Nset1 = X{iCnd}.Nset1;
    Nset2 = X{iCnd}.Nset2;
    Cset = X{iCnd}.Cset;
    Rmat = X{iCnd}.Rmat;
    [temp,prmat{iCnd}] = AhyBCL_datalike1(theta,Nset1,Nset2,Cset,Rmat,model,Nsamples);
    if ~isempty(temp); ll(iCnd) = temp; end
end

loglike = sum(ll);

% Generate fake data if required
if fakerun && nargout > 2
    Ntrials = 24*ones(18,1);
    for iCnd = 1:numel(X)
        X{iCnd}.Rmat = mnrnd(Ntrials, prmat{iCnd});
    end
end

end