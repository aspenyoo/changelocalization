function combineLLcube(model,subjid,nParts,marginalizeover,gridMat,filetail)
if nargin < 4; marginalizeover = 1; end
if nargin < 5; gridMat = [-0.69 3.5 30; -0.69 3.5 30; -0.69 3.5 30]; end
if nargin < 6; filetail = []; end % file index if there are multiple LLmats per subject and model. should be a string

% ========= INPUT VARIABLE ========
% MODEL: (1) optimal. (2) fixed
% SUBJID: subject ID
% NPARTS: how many parts the LL cube was chopped up into
% MARGINALIZEOVER: for plotting purposes. which variable you want to
% maginalized over. (1) Jbar1 (2) Jbar2 (3) tau
% GRIDMAT: 3 x 3 matrix. the first two columns correspond to the lowest and
% highest LOG value of the parameter, respectively. The last column
% corresponds to the number of grids between them (i.e., last term of
% linspeace). Each row corresponds to a parameter. 
% FILETAIL: should be of format '_blah', in which blah is some tail that
% the data file will have. 

% model = 1;
% nParts = 100;
% subjid = 'PM';




try
    load(sprintf('LLcube_model%d_subjid%s%s.mat',...
        model,subjid,filetail));
    v2struct(data);
catch
    
    Jbar1Vec = exp(linspace(gridMat(1,1),gridMat(1,2),gridMat(1,3)));
    Jbar2Vec = exp(linspace(gridMat(2,1),gridMat(2,2),gridMat(2,3)));
    tauVec = exp(linspace(gridMat(3,1),gridMat(3,2),gridMat(3,3)));
    totalCombos = prod(gridMat(:,3));
    
    
    LLMat = nan(totalCombos,2);
    for ipart = 1:nParts;
        load(sprintf('sep_LLcube_model%d_subjid%s_ipart%d_npart%d.mat',...
            model,subjid,ipart,nParts));
        LLMat((ceil((totalCombos/nParts)*(ipart-1))+1):ceil((totalCombos/nParts)*ipart),:) = loglikeMat;
    end
    
    LLMat = sortrows(LLMat);
    LLMat = LLMat(:,2);
    LLMat = reshape(LLMat,gridMat(:,3)');
    
    % save relevant info
    fieldNames = {'subjid','model','Jbar1Vec','Jbar2Vec','tauVec','gridMat','LLMat','fieldNames'};
    data = v2struct(fieldNames);
    save(sprintf('LLcube_model%d_subjid%s%s.mat',...
        model,subjid,filetail),'data');
end


% =============== PLOT LLmats ==================

blah = 1;
% marginalizeover = 3; % which dimension to marginalize over
keys = [KbName('left') KbName('right') KbName('esc')];
idx = 1;
clims = [min(LLMat(:)) max(LLMat(:))];
while (blah)
    switch marginalizeover
        case 1
            imagesc(Jbar2Vec,tauVec,squeeze(LLMat(idx,:,:)),clims); defaultplot;
%             colormap('Copper')
            title(sprintf('Jbar1 = %02d',Jbar1Vec(idx)))
            xlabel('Jbar2'); ylabel('tau');
            pause(0.1)
            pressedKey = waitForKeys(keys);
        case 2
            imagesc(Jbar1Vec,tauVec,squeeze(LLMat(:,idx,:)),clims); defaultplot;
%             colormap('Copper')
            title(sprintf('Jbar2 = %d',Jbar2Vec(idx)))
            xlabel('Jbar1'); ylabel('tau')
            pause(0.1)
            pressedKey = waitForKeys(keys);
        case 3
            imagesc(Jbar1Vec,Jbar2Vec,squeeze(LLMat(:,:,idx)),clims); defaultplot;
%             colormap('Copper')
            title(sprintf('tau = %d',tauVec(idx)))
            xlabel('Jbar1'); ylabel('Jbar2')
            pause(0.1)
            pressedKey = waitForKeys(keys);
            
    end
    
    % adjust which slice to be looking at
    if pressedKey == 1;
        if idx > 1;
            idx = idx - 1;
        else
            display('lowest bound hit')
        end
    elseif pressedKey == 2;
        if idx <  gridMat(marginalizeover,3);
            idx = idx + 1;
        else
            display('upper bound hit')
        end
    elseif pressedKey == length(prefs.keys);
        blah = 0;
    end
end
end

function pressedKey = waitForKeys(keys)

pressedKey=0;
while (1)
    
    [~, ~, keyCode] = KbCheck();
    if  any(keyCode(keys))
        pressedKey = find(keyCode(keys));
        break;
    end
    
end
end
