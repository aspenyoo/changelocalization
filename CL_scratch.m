%% simulate data

model = 1; % 1: optimal. 2: fixed
X = []; % data in luigiform
Nsamples = []; % default if empty

theta = [3 10 1 0.01]; % [jbar1 jbar2 tau lapse]

[~,prmat,X] = AhyBCL_datalikeall(theta,X,model,Nsamples);

%% logistic regression: pc ~ condition + delta
clear all

subjid = 'AHY';
load(sprintf('data_ChangeLocalization_subj%s.mat',subjid));

deltaVec = 2.5:5:87.5;              % possible change values
respVec = [3 0 4 0 0 0 2 0 1];      % go from response key to location index

designMat = designMat(designMat(:,2) ~= 2,:);

% set X and y matrices
y = designMat(:,7); % correct/not correct
X = [designMat(:,2) abs(designMat(:,1))]; % delay, delta
X(X(:,1) == 0,2) = 1; % delay for 0 condition is 1
X(X(:,2) == 4,2) = 3; % delay for 4 condition is 3

% logistic regression
mdl = fitglm(X,y,'Distribution','binomial')
ci = coefCI(mdl)

%% linear regression: confidence ~ condition + delta

clear all

subjid = 'AHY';
load(sprintf('data_ChangeLocalization_subj%s.mat',subjid));

deltaVec = 2.5:5:87.5;              % possible change values
respVec = [3 0 4 0 0 0 2 0 1];      % go from response key to location index

designMat = designMat(designMat(:,2) ~= 2,:);

% set X and y matrices
y = designMat(:,9); % confidence rating
X = [designMat(:,2) abs(designMat(:,1))]; % delay, delta
X(X(:,1) == 0,2) = 1; % delay for 0 condition is 1
X(X(:,2) == 4,2) = 3; % delay for 4 condition is 3

% linear regression
% X = [ones(size(designMat,1),1) X];
% [b,bint,~,~,stats] = regress(y,X) % stats: R^2, F, p, and estimate of error variance

% ordinal multinomial regression
[B,dev,stats] = mnrfit(X,y,'model','ordinal')
LL = stats.beta - 1.96*stats.se
UL = stats.beta + 1.96*stats.se

% plot regression fits
xx = linspace(5,85,50);
XX = [ones(50,1) 3*ones(50,1) xx'];
BMat = [B(1:3)'; repmat(B(4:end),1,3)];
logoddsMat = XX*BMat; % log odds for each of the 3 equations
logoddssign = sign(logoddsMat);
difflogoddssign = diff(logoddssign,1,2); % difference between the signs of log odds
[row,col] = find(difflogoddssign == 2);
predictedConf = nan(50,1);
predictedConf(row) = col+1;
predictedConf(logoddssign(:,1) == 1 ) = 1;
predictedConf(logoddssign(:,3) == -1) = 4;
hold on; plot(xx,predictedConf);


%% plot psychometric fn


