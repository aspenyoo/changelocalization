%% simulate data

model = 1; % 1: optimal. 2: fixed
X = []; % data in luigiform
Nsamples = []; % default if empty

theta = []; % [jbar1 jbar2 tau lapse]

[loglike,prmat,X] = AhyBCL_datalikeall(theta,X,model,Nsamples);