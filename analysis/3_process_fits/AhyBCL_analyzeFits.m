function [LL,xbest,table] = AhyBCL_analyzeFits()

table = [];

files = dir('*.mat');
for i = 1:numel(files)
    filename = files(i).name;
    temp = load(filename);
    if isempty(table)
        table = temp.table;
    else
        % Copy not-NaN fields
        for ff = fields(table)'
            idx = ~isnan(temp.table.(ff{:}));
            table.(ff{:})(idx) = temp.table.(ff{:})(idx);
        end
    end
end

[nsubjs,nmodels,nruns] = size(table.LLbest);
nvars = size(table.xbest,4);

LL = NaN(nsubjs,nmodels);

nTrials = 1296;

% Make model and subject LL table
for i = 1:nsubjs
    for j = 1:nmodels
        [LLbest,idx] = max(table.LLbest(i,j,:));
        meanLLbest = nanmean(table.LLbest(i,j,:));
        stdLLbest = nanstd(table.LLbest(i,j,:));
        LL(i,j) = LLbest;
        aicc(i,j) = -2*LLbest + 2*nvars + 2*nvars*(nvars+1)/(nTrials-nvars-1);
        bic(i,j) = -2*LLbest + nvars*log(nTrials);
        xbest(i,j,:) = table.xbest(i,j,idx,:);
        % fprintf('Subject %d, model %d. LLbest = %.2f (mean ± SD: %.2f ± %.2f), AICc = %.2f, BIC = %.2f.\n',i,j,LLbest,meanLLbest,stdLLbest,aicc(i,j),bic(i,j));
        fprintf('Subject %d, model %d. LLbest = %.2f, AICc = %.2f, BIC = %.2f.\n',i,j,LLbest,aicc(i,j),bic(i,j));
        
    end
end

fprintf('Delta BIC: ');
for i = 1:nsubjs; fprintf('subj %d = %.2f; ', i, bic(i,2) - bic(i,1)); end
fprintf('\n');
