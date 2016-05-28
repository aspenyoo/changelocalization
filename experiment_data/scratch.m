%% testing stimuli
Screen('Preference', 'SkipSyncTests', 1);

grey = 128;

% screen info (visual)
screenNumber = max(Screen('Screens'));       % use external screen if exists
[w, h] = Screen('WindowSize', screenNumber)  % screen resolution of smaller display
screenResolution = [w h];                 % screen resolution

% open screen
windowPtr = Screen('OpenWindow',screenNumber,grey,[],32,2);
Screen(windowPtr,'BlendFunction',GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

nPC = 10;
pc = [60 60 60 60 60 70 80 90 80 70];

Screen('fillRect',windowPtr,128*ones(1,3));
Screen('DrawText',windowPtr,['You have earned ' num2str(pc(end)) ' points. Press any key to continue'],250,h/2 - 50,[255 255 255]);

% make graph of progress
graphwidth = round(w/3); % half of graph width
graphheight = round(h/5); % half of graph height
origin = [w/2-graphwidth/2 h/2+graphheight];
Screen('DrawLine',windowPtr,255*ones(1,3),origin(1), origin(2),origin(1)+graphwidth,origin(2)); % x-axis
Screen('DrawLine',windowPtr,255*ones(1,3),origin(1), origin(2),origin(1),origin(2)-graphheight); % y-axis
Screen('DrawText',windowPtr,'your points over time',w/2,origin(2)+100,[255 255 255]); % xlabel

% draw lines to make graph
pc = (pc-25)./100*graphheight; % rescaling to size of graph
dx = graphwidth./nPC; % equally spaced out to fill size of graph at end of exp. 
og = origin;
dc = nan(2,length(pc));
dc(:,1) = [og(1),og(2)-pc(1)];
for ipc = 2:length(pc); % draw the lines
    dc(:,ipc) = [og(1)+dx og(2)-pc(ipc)];
    Screen('DrawLine',windowPtr,255*ones(1,3),dc(1,ipc-1),dc(2,ipc-1),dc(1,ipc),dc(2,ipc)); % y-axis
    og(1) = og(1) + dx;
end
Screen('DrawDots',windowPtr,dc,5,255*ones(1,3),[],1)

% flip screen
Screen('Flip', windowPtr);
pause; 


% reliability = 1;
% StimLength = 100;
% StimSizes = [StimLength StimLength];
% im = round(makeGabor(StimLength, reliability).*255);
% StimPatch = Screen('MakeTexture',windowPtr,im);
% cuesrcrect = [0 0 StimSizes];
% destrect = CenterRectOnPoint(cuesrcrect,w/2,h/2);
% Screen('DrawTexture', windowPtr, StimPatch, cuesrcrect, destrect);
% Screen('Flip', windowPtr);
% pause; 

% linee = nan(2,1);
% for i = 1:180;
%     [linee(1), linee(2)] = lineCoord(100,i);
%     xy = [linee -linee];
%     Screen('DrawLines',windowPtr, xy, 4, [255 255 255] ,[w/2 h/2],1);
%     Screen('Flip', windowPtr);
%     pause(0.05);
% end
% 
sca;
%% 

keys = [KbName('9') KbName('7') KbName('1') KbName('3') KbName('esc')];

pressedKey=0;
while (1)
    
    [keyIsDown, secs, keyCode] = KbCheck()
    if  any(keyCode(keys))
        RT = GetSecs - tstart;
        pressedKey = find(keyCode(keys));
        break;
    end
    
end