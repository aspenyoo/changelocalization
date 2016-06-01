%% simulate data

model = 1; % 1: optimal. 2: fixed
X = []; % data in luigiform
Nsamples = []; % default if empty

theta = [3 10 1 0.01]; % [jbar1 jbar2 tau lapse]

[~,prmat,X] = AhyBCL_datalikeall(theta,X,model,Nsamples);

