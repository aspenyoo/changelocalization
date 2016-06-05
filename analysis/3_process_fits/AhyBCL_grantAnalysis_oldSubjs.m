%AHYBCL_GRANTANALYSIS Analysis for the five original subjects

expflag = 0;

subjs = {'PM','AD','SM','SC','DG'};
models = 1:2;
Nsubjs = numel(subjs);
Nmodels = numel(models);
post = [];

% Load LL cube for each subject and model
for i = 1:Nsubjs
    for j = 1:Nmodels
        m = models(j);
        filename = ['LLcube_model' num2str(m) '_subjid' subjs{i} '.mat'];
        load(filename,'data');
        post{i,j} = AhyBCL_analyzeLLCube(data);
        mlike(i,j) = post{i,j}.mlike;   % Marginal likelihood
    end
end
    
fprintf('Difference in marginal likelihood per subject (higher is better for model 2)\n')
 
diff(mlike,[],2)

colmat = [  ...
    141 211 199; ... 
    251 128 114; ...   
    128 177 211; ...    
    253 180 98; ...    
    160 120 100; ...   
    70 70 233; ...    
    252 205 229; ...    
    165 211 195; ...        
    159 212 105; ...          
    188 128 189; ...        
    212 148 169; ...                    
    0 0 0 ...                   
    ]/255;

% Plot marginal posteriors
for i = 1:Nsubjs
    for j = 1:2        
        for k = 1:3
            subplot(2,3,(j-1)*3 + k);
            theta = post{i,j}.thetaVec{k};
            pdf = post{i,j}.postpdf{k};
            if expflag
                theta = exp(theta);
                pdf = pdf./theta;
            end
            
            plot(theta,pdf,'LineWidth',1,'Color',colmat(i,:));
            hold on;
            xlabel(post{1,j}.params{k});
            title(['Model ' num2str(models(j))]);
            box off;
            set(gca,'TickDir','out');
            if k == 1; ylabel('Posterior probability density'); end
            idx = 3:3:30;
            set(gca,'Xtick',theta(idx));
            
            if ~expflag
                for iidx = 1:numel(idx)
                    ticks{iidx} = num2str(exp(post{1,j}.thetaVec{k}(idx(iidx))),'%.2g');
                end
                set(gca,'Xticklabel',ticks);
            end
        end
    end
end
set(gcf,'Color','w');