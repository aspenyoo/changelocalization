% plot psychometric function for single or multiple subjects

subjids = {'PM','SM','DG','SC'};
condidx = [4 1 2 3];

delta = 5:5:90;
nCond = 4;%length(data);
nSubj = length(subjids);

figure;
pc = nan(nSubj,nCond,length(delta));
pic = nan(nSubj,nCond,length(delta));
pii = nan(nSubj,nCond,length(delta));
for isubj = 1:nSubj;
    
    % load file
    subjid = subjids{isubj};
    load(sprintf('data_ChangeLocalization_subj%s_04282016.mat',subjid));
    
    for icond = 1:nCond;
        
        v2struct(data{icond})
        
        % plot correct
        pc(isubj,icond,:) = Rmat(:,1)./sum(Rmat,2);
        
        % plot incorrect in correct interval
        pic(isubj,icond,:) = Rmat(:,2)./sum(Rmat,2);
        
        % plot incorrect in incorrect interval
        pii(isubj,icond,:) = Rmat(:,3)./sum(Rmat,2);
        
    end
end

combinebins = 2;
nBins_new = length(delta)/combinebins;
pc_new = nan(nSubj,nCond,nBins_new);
pic_new = nan(nSubj,nCond,nBins_new);
pii_new = nan(nSubj,nCond,nBins_new);
deltanew = nan(1,nBins_new);
for ibin = 1:nBins_new;
    deltanew(ibin) = mean(delta((ibin-1)*combinebins+1 : ibin*combinebins));
    pc_new(:,:,ibin) = mean(pc(:,:,(ibin-1)*combinebins+1 : ibin*combinebins),3);
    pic_new(:,:,ibin) = mean(pic(:,:,(ibin-1)*combinebins+1 : ibin*combinebins),3);
    pii_new(:,:,ibin) = mean(pii(:,:,(ibin-1)*combinebins+1 : ibin*combinebins),3);
    
end
pc = pc_new;
pic = pic_new;
pii = pii_new;
delta = deltanew;

% plot 0 and 4 condition percent correct
subplot(2,2,1);
hold on;
errorbar(delta,squeeze(mean(pc(:,1,:),1)),squeeze(std(pc(:,1,:),[],1))/sqrt(nSubj),'Color','k','LineStyle','-');
errorbar(delta,squeeze(mean(pc(:,4,:),1)),squeeze(std(pc(:,4,:),[],1))/sqrt(nSubj),'Color',0.7*ones(1,3),'LineStyle','-');
legend('4/1','0/2')
defaultplot;
title('pc for 0 and 4 condition')

% plot pc psychometric functions for the two 2 conditions
subplot(2,2,2);
hold on;
errorbar(delta,squeeze(mean(pc(:,2,:),1)),squeeze(std(pc(:,2,:),[],1))/sqrt(nSubj),'Color','r','LineStyle','--');
errorbar(delta,squeeze(mean(pc(:,3,:),1)),squeeze(std(pc(:,3,:),[],1))/sqrt(nSubj),'Color','b','LineStyle','-');
legend('2/1','2/2')
defaultplot
title('pc for 2 conditions')

% plot pic psychometric functions for both 2 condtions
subplot(2,2,3);
hold on;
errorbar(delta,squeeze(mean(pic(:,2,:),1)),squeeze(std(pic(:,2,:),[],1))/sqrt(nSubj),'Color','r','LineStyle','--');
errorbar(delta,squeeze(mean(pic(:,3,:),1)),squeeze(std(pic(:,3,:),[],1))/sqrt(nSubj),'Color','b','LineStyle','-');
legend('2/1','2/2')
title('incorrect item, correct interval')
defaultplot

% plot pic psychometric functions for both 2 condtions
subplot(2,2,4);
hold on;
errorbar(delta,squeeze(mean(pii(:,2,:),1))./2,squeeze(std(pii(:,2,:),[],1))/sqrt(nSubj),'Color','r','LineStyle','--');
errorbar(delta,squeeze(mean(pii(:,3,:),1))./2,squeeze(std(pii(:,3,:),[],1))/sqrt(nSubj),'Color','b','LineStyle','-');
legend('2/1','2/2')
title('incorrect item, incorrect interval')
defaultplot


% average pc across conditions
ave_pc = squeeze(mean(mean(pc,1),3))

