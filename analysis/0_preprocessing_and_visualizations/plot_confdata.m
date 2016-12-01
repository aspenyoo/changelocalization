% close all
clear all
clc

subjids = {'ALM'};

% deltaVec = 5:5:90;
deltabinVec = [0 10:10:90];
nDeltaBins = length(deltabinVec)-1;
nCond = 4;
nConf = 4;
nSubj = length(subjids);

RTthreshold = 120; % cutoff any trial with RT greater than this number

NsetVec = [0 2 2 4];            % number of stimuli in first display
CsetVec = [2 1 2 1];            % which display change target was presented

conf = nan(nSubj,nDeltaBins); conf_std = conf;
PC = nan(nSubj,nConf); PC_std = nan(nSubj,nConf);
confcond = nan(nSubj,nCond,nDeltaBins); confcond_std = confcond;
PCcond = nan(nSubj,nCond,nConf); PCcond_std = PCcond;
bincenters = nan(1,nDeltaBins);
for isubj = 1:nSubj;
    
    % load file
    subjid = subjids{isubj};
    load(sprintf('data_ChangeLocalization_subj%s.mat',subjid));
    
    stimuliMat(designMat(:,8) > RTthreshold,:) = [];
    designMat(designMat(:,8) > RTthreshold,:) = [];
    
        % get mean confidence for each delta
        for idelta = 1:nDeltaBins;
            bincenters(idelta) = mean(deltabinVec(idelta:idelta+1));
            
            idx = (abs(designMat(:,1)) <= deltabinVec(idelta+1)) & (abs(designMat(:,1)) > deltabinVec(idelta)); % how many of that bin
            conf(isubj,idelta) = mean(designMat(idx,9));
            conf_std(isubj,idelta) = std(designMat(idx,9))./sum(idx);
            
            for icond = 1:nCond;
                Nset1 = NsetVec(icond);
                Cset = CsetVec(icond);
                
                % finding which trials are in the correct conditions.
                idxx = idx & (designMat(:,2) == Nset1) & (diag(stimuliMat(:,stimuliMat(:,9))) == Cset);
                
                
                confcond(isubj,icond,idelta) = mean(designMat(idxx,9));
                confcond_std(isubj,icond,idelta) = std(designMat(idxx,9))./sqrt(sum(idxx));
            end 
        end
        
        % get PC for each confidence value
        for iconf = 1:nConf;
            
            idx = designMat(:,9) == iconf;
            PC(isubj,iconf) = mean(designMat(idx,7));
            PC_std(isubj,iconf) = std(designMat(idx,7))./sum(idx);
            
            for icond = 1:nCond;
                Nset1 = NsetVec(icond);
                Cset = CsetVec(icond);
                
                % finding which trials are in the correct conditions.
                idxx = idx & (designMat(:,2) == Nset1) & (diag(stimuliMat(:,stimuliMat(:,9))) == Cset);

                PCcond(isubj,icond,iconf) = mean(designMat(idxx,7));
                PCcond_std(isubj,icond,iconf) = std(designMat(idxx,7))./sum(idxx);
                %             PC_std(isubj,icond,iconf) = mean(designMat(idxx,7)).*(1-mean(designMat(idxx,7)))./sum(designMat(idxx,7));
                %             PC_std(isubj,icond,iconf) = (1-mean(designMat(idxx,7)))./sum(idxx);
            end
        end
end
if nSubj > 1;
    PC_std = std(PC,[],1)./sqrt(nSubj);
    PCcond_std = std(PCcond,[],1)./sqrt(nSubj);
    conf_std = std(conf,[],1)./sqrt(nSubj);
    confcond_std = std(confcond,[],1)./sqrt(nSubj);
end


% ========= PLOTS! ===========

% % mean confidence of delta for each condition
% figure;
% if nSubj == 1;
%     for icond = 1:nCond;
%         subplot(2,2,icond);
%         errorbar(bincenters,squeeze(confcond(1,icond,:)),squeeze(confcond_std(1,icond,:)),'Color','k');
%         defaultplot
%         title(sprintf('condition %d/%d',NsetVec(icond),CsetVec(icond)))
%         xlabel('delta')
%         ylabel('mean confidence')
%         ylim([1 4])
%         ax = gca;
%         ax.YTick = 1:4;
%     end
% end

% confidence for just 0/2 and 4/1 conditions
figure; hold on;
errorbar(bincenters,squeeze(confcond(1,4,:)),squeeze(confcond_std(1,4,:)),'Color',0.5*ones(1,3));
errorbar(bincenters,squeeze(confcond(1,1,:)),squeeze(confcond_std(1,1,:)),'Color','k');
% errorbar(1:4,squeeze(PCcond(1,2,:)),squeeze(PCcond_std(1,2,:)),'Color','r');
% errorbar(1:4,squeeze(PCcond(1,3,:)),squeeze(PCcond_std(1,3,:)),'Color','b');
defaultplot
xlabel('delta')
ylabel('confidence')
xlabel('delta')
ylabel('mean confidence')
ylim([1 4])
ax = gca;
ax.YTick = 1:4;
legend('first interval','second interval')



% % mean confidence of delta across condition
% figure;
% errorbar(bincenters,squeeze(conf),squeeze(conf_std),'Color','k')
% defaultplot
% title('mean confidence across condition')
% xlabel('delta')
% ylabel('mean confidence')
% ylim([1 4])
% ax = gca;
% ax.YTick = 1:4;
% 
% % PC across confidence
% figure;
% for icond = 1:nCond;
%     subplot(2,2,icond)
%     errorbar(1:4,squeeze(PCcond(1,icond,:)),squeeze(PCcond_std(1,icond,:)),'Color','k');
%     defaultplot
%     title(sprintf('condition %d/%d',NsetVec(icond),CsetVec(icond)))
%     xlabel('confidence')
%     ylabel('proportion correct')
%     ylim([.25 1])
%     ax = gca;
%     ax.YTick = [0.25 0.5 0.75 1];
%     ax.XTick = 1:4;
% end
% 
% % PC for just 0/2 and 4/1 conditions
% figure; hold on;
% errorbar(1:4,squeeze(PCcond(1,4,:)),squeeze(PCcond_std(1,4,:)),'Color',0.5*ones(1,3));
% errorbar(1:4,squeeze(PCcond(1,1,:)),squeeze(PCcond_std(1,1,:)),'Color','k');
% % errorbar(1:4,squeeze(PCcond(1,2,:)),squeeze(PCcond_std(1,2,:)),'Color','r');
% % errorbar(1:4,squeeze(PCcond(1,3,:)),squeeze(PCcond_std(1,3,:)),'Color','b');
% defaultplot
% xlabel('confidence')
% ylabel('proportion correct')
% ylim([.25 1])
% ax = gca;
% ax.YTick = [0.25 0.5 0.75 1];
% ax.XTick = 1:4;
% legend('first interval','second interval')
% 
% % PC of confidence across condition
% figure;
% errorbar(1:4,squeeze(PC),squeeze(PC_std),'Color','k')
% defaultplot
% title('PC for each confidence level')
% xlabel('confidence')
% ylabel('proportion correct')
% ylim([.25 1])
% ax = gca;
% ax.YTick = [0.25 0.5 0.75 1];
% ax.XTick = 1:4;