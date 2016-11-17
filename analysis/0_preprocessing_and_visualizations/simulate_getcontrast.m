function simulate_getcontrast

% ================= REAL PARAMETER VALUES ==============

% generate parameters from reasonable values (loosely based on carrasco paper)
lambda = .4; % lapse or 1 - max performance
mu = 10; % percent. not in log
mu = log(mu);
sigma = 1; % in log space
gamma = 0.25; % chance level

xx = linspace(0,log(100),100);
yy = psyfun_pcorrect(xx,mu,sigma,lambda,gamma,@psynormcdf);

plot(xx,yy,'k-')
ylim([0.25 0.8])


%% ============== MOC ======================
nsteps = 10;
xx = linspace(0,log(100),nsteps);
ntrialsperstep = 500/nsteps;
xx = repmat(xx,1,ntrialsperstep);
resps = getresponse(xx);        % responses

x = fminsearch(@obj_func,rand(1,4));



% =================== PSYBAYES ==================




    function resp =  getresponse(stim)
        pcorr = psyfun_pcorrect(stim,mu,sigma,lambda,gamma,@psynormcdf);
        resp = pcorr > rand(size(pcorr));
    end

    function nll = obj_func(x)
        pcorr = psyfun_pcorrect(xx,x(1),x(2),x(3),x(4),@psynormcdf);
        nll = resps.*pcorr + (1-resps).*(1-pcorr);
        
    end
end
