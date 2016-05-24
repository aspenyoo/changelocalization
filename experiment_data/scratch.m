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

reliability = 1;
StimLength = 100;
StimSizes = [StimLength StimLength];
im = round(makeGabor(StimLength, reliability).*255);
StimPatch = Screen('MakeTexture',windowPtr,im);
cuesrcrect = [0 0 StimSizes];
destrect = CenterRectOnPoint(cuesrcrect,w/2,h/2);
Screen('DrawTexture', windowPtr, StimPatch, cuesrcrect, destrect);
Screen('Flip', windowPtr);
pause; 

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