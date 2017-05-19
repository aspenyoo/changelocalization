function [pc, pc_std] = plot_psychometricfunction(subjids,exptype,condVec)
% plot psychometric function for single or multiple subjects

% subjids = {'AHY'};

delta = 2.5:5:87.5;
nCond = 4;
nSubj = length(subjids);

pc = nan(nSubj,nCond,length(delta));
pic = nan(nSubj,nCond,length(delta));
pii = nan(nSubj,nCond,length(delta));

pc_std = nan(nCond,length(delta)); pic_std = pc_std; pii_std = pc_std;
for isubj = 1:nSubj
    
    % load file
    subjid = subjids{isubj};
    load(sprintf('processeddata_ChangeLocalization_%s_subj%s.mat',exptype,subjid));
    
    for icond = 1:nCond
        
        v2struct(data{icond})
        
        % plot correct
        pc(isubj,icond,:) = Rmat(:,1)./sum(Rmat,2);
        
        
        % plot incorrect in correct interval
        pic(isubj,icond,:) = Rmat(:,2)./sum(Rmat,2);
        
        
        % plot incorrect in incorrect interval
        pii(isubj,icond,:) = Rmat(:,3)./sum(Rmat,2);
        
        if nSubj == 1
            pc_std(icond,:) = sqrt(squeeze(pc(isubj,icond,:).*(1-pc(isubj,icond,:))));
            pic_std(icond,:) = sqrt(squeeze(pic(isubj,icond,:)).*squeeze((1-pic(isubj,icond,:))));
            pii_std(icond,:) = sqrt(squeeze(pii(isubj,icond,:).*(1-pii(isubj,icond,:))));
        end
    end
end

combinebins = 2;
nBins_new = length(delta)/combinebins;
pc_new = nan(nSubj,nCond,nBins_new);
pic_new = nan(nSubj,nCond,nBins_new);
pii_new = nan(nSubj,nCond,nBins_new);
deltanew = nan(1,nBins_new);
for ibin = 1:nBins_new
    deltanew(ibin) = mean(delta((ibin-1)*combinebins+1 : ibin*combinebins));
    pc_new(:,:,ibin) = mean(pc(:,:,(ibin-1)*combinebins+1 : ibin*combinebins),3);
    pic_new(:,:,ibin) = mean(pic(:,:,(ibin-1)*combinebins+1 : ibin*combinebins),3);
    pii_new(:,:,ibin) = mean(pii(:,:,(ibin-1)*combinebins+1 : ibin*combinebins),3);
end

pc = pc_new;
pic = pic_new;
pii = pii_new;
delta = deltanew;

if nSubj == 1;
    pc_std = sqrt(squeeze(pc(isubj,:,:).*(1-pc(isubj,:,:))./(36.*(pc(isubj,:,:)))));
    pic_std = sqrt(squeeze(pic(isubj,:,:).*(1-pic(isubj,:,:))./(36.*(pc(isubj,:,:)))));
    pii_std= sqrt(squeeze(pii(isubj,:,:).*(1-pii(isubj,:,:))./(36.*(pc(isubj,:,:)))));
else
    pc_std = squeeze(std(pc,[],1))./sqrt(nSubj);
    pic_std = squeeze(std(pic,[],1))./sqrt(nSubj);
    pii_std = squeeze(std(pii,[],1))./sqrt(nSubj);
end

% if I don't want errorbars plotted
errorbarplot = 1;
if ~(errorbarplot)
    pc_std = zeros(nBins_new,length(delta)); pic_std = pc_std; pii_std = pc_std;
end

% subplot(2,2,1);
% colors = zeros(4,3);
colors = aspencolors(4,'passport');%['b'; 'y'; 'g'; 'r'];
% figure; hold on;
hold on;
legendd = cell(1,length(condVec));
for icond = 1:length(condVec)
    errorb(delta,squeeze(mean(pc(:,condVec(icond),:),1)),pc_std(condVec(icond),:),...
        'Color',colors(condVec(icond),:));%,'LineStyle','none');
    legendd{icond} = ['cond' num2str(condVec(icond))];
end
% errorbar(delta,squeeze(mean(pc(:,4,:),1)),pc_std(4,:),'Color',0.7*ones(1,3),'LineStyle','-');
% legend(sprintf('condition %d/%d',data{condVec(icond)}.Nset1,data{condVec(icond)}.Nset1),'1st interval')
defaultplot;
% title('pc for 0 and 4 condition')
ax = gca;
ax.XTick = [0 30 60 90];
ax.YTick = [0 0.5 1];
% legend(legendd)

% % plot pc psychometric functions for the two 2 conditions
% subplot(2,2,2);
% hold on;
% errorbar(delta,squeeze(mean(pc(:,2,:),1)),pc_std(2,:),'Color','r','LineStyle','--');
% errorbar(delta,squeeze(mean(pc(:,3,:),1)),pc_std(3,:),'Color','b','LineStyle','-');
% legend('2/1','2/2')
% defaultplot
% title('pc for 2 conditions')
% 
% % plot pic psychometric functions for both 2 condtions
% subplot(2,2,3);
% hold on;
% errorbar(delta,squeeze(mean(pic(:,2,:),1)),pic_std(2,:),'Color','r','LineStyle','--');
% errorbar(delta,squeeze(mean(pic(:,3,:),1)),pic_std(3,:),'Color','b','LineStyle','-');
% legend('2/1','2/2')
% title('incorrect item, correct interval')
% defaultplot
% 
% % plot pic psychometric functions for both 2 condtions
% subplot(2,2,4);
% hold on;
% errorbar(delta,squeeze(mean(pii(:,2,:),1))./2,pii_std(2,:),'Color','r','LineStyle','--');
% errorbar(delta,squeeze(mean(pii(:,3,:),1))./2,pii_std(3,:),'Color','b','LineStyle','-');
% legend('2/1','2/2')
% title('incorrect item, incorrect interval')
% defaultplot


% % average pc across conditions
% ave_pc = squeeze(mean(mean(pc,1),3))

