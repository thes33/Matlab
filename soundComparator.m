function soundComparator(soundFiles)
% function soundComparator(soundFiles)
%
% A GUI for listening and systematic listening and comparing sound recordings.
%
% INPUT:
%  soundFiles = (optional) Cell array of a list of sound file paths to compare.  Alternatively,
%               a directory can be given to compare all .wav files in that directory.
%               [Deafult]: Use the current directory.
%
%
% @Author: Eugene Brandewie Feb 11, 2022
% @Email: eugene.brandewie@gmail.com
% @Version: 1.0



% VERSION HISTORY
%------------------------------------------
% 1.0 - Initial Release




if (nargin < 1)
    soundFiles = cd;
end



%% Process Sound Files
%===========================================

% IF directory
if (ischar(soundFiles) & isdir(soundFiles))
    sfDir = soundFiles;
    FILES = dir(soundFiles);
    soundFiles = {};
    cnt = 0;
    % FOR each file in directory
    for (ff = 1:length(FILES))
        [fPath, fName, fExt] = fileparts(FILES(ff).name);
        % IF a wav file
        if (strcmp(fExt,'.wav'))
            cnt = cnt + 1;
            soundFiles(cnt) = {[sfDir filesep FILES(ff).name]};
            fileNames(cnt) = {FILES(ff).name};
        end
    end
    
    
else % IF list of sound files
    
    for (ff = 1:length(soundFiles))
        [fPath, fName, fExt] = fileparts(soundFiles{ff});
        fileNames(ff) = {fName};
    end    
    
end

numFiles = length(soundFiles);




%% Load Sound Files
%=============================================
for (ff = 1:length(soundFiles))
    [wav, fs] = audioread(soundFiles{ff});
    WAVs(ff) = {wav};
    FS(ff) = fs;
end




%% Create GUI
%=============================================

% Create GUI
H = createSoundComparatorGUI(numFiles);
data = guidata(H);
data.soundFiles = soundFiles;
data.fileNames = fileNames;
data.WAVs = WAVs;
data.FS = FS;
data.START = 0;
data.STOP = -1;
guidata(H,data);

set(H,'CloseRequestFcn',@SoundCompareGUI_CloseFcn)

start(data.updateTimer);


end




function H = createSoundComparatorGUI(numFiles)
%% Draw UI
%=======================================================================

% Create Figure
H = figure('Name','Sound Comparator GUI','Toolbar','none');
axis off;
hold on;

% Retrieve data
data = guidata(gcf);


% Mail Panel
data.mainPanel = uipanel( ...
    'Title', '', ...
    'Units','normalized', ...
    'Position',[0.01 0.01 0.99 0.99]);


% Current Wav Text
data.currentFileLabel = uicontrol('Parent',data.mainPanel,'Style','text', 'String','...', 'Units','normalized', ...
    'Position',[0.1 .85 0.8 0.1], 'HorizontalAlignment', 'center', ...
    'BackgroundColor','white', ...
    'FontSize',22);

numPerRow = max(numFiles/3,1);
btnWidth = 0.9 / numPerRow;
lstBtnWidth = 0.9 / (numFiles - (numPerRow*2));

% Create sound file buttons
for (bb = 1:numFiles)
    btnPos(bb,1) = 0.05 + ((bb-1)*btnWidth);
    btnPos(bb,2) = 0.65;
    btnPos(bb,3) = btnWidth;
    if (bb > numPerRow)
        btnPos(bb,1) = 0.05 + ((bb-1-numPerRow)*btnWidth);
        btnPos(bb,2) = 0.55;
        btnPos(bb,3) = btnWidth;
    end
    if (bb > (numPerRow*2))  % Last Row
        btnPos(bb,1) = 0.05 + ((bb-1-numPerRow-numPerRow)*lstBtnWidth);
        btnPos(bb,2) = 0.45;
        btnPos(bb,3) = lstBtnWidth;
    end
    % Create button
    data.buttons(bb) = uicontrol('Parent',data.mainPanel,'Style','pushbutton', 'String',num2str(bb), 'Units','normalized', ...
    'Position',[btnPos(bb,1) btnPos(bb,2) btnPos(bb,3) 0.1], 'HorizontalAlignment', 'left', ...
    'FontSize',26, 'Callback',@SoundCompareGUI_SoundFileBtnCallback, ...
    'UserData',bb);
end

% Stop Button
data.stopBtn = uicontrol('Parent',data.mainPanel,'Style','pushbutton', 'String','Stop', 'Units','normalized', ...
    'Position',[0.85 0.2 0.1 0.2], 'HorizontalAlignment', 'left', ...
    'FontSize',26, 'Callback',@SoundCompareGUI_StopBtnCallback);

% Start Slider
data.startSlider = uicontrol('Parent',data.mainPanel,'Style','slider', 'Units','normalized', ...
    'Position',[0.05 0.3 0.78 0.1], 'HorizontalAlignment', 'left', ...
    'BackgroundColor','white', ...
    'Value', 0, ...
    'FontSize',26, 'Callback',@SoundCompareGUI_StartSliderCallback);

% Stop Slider
data.stopSlider = uicontrol('Parent',data.mainPanel,'Style','slider', 'Units','normalized', ...
    'Position',[0.05 0.19 0.78 0.1], 'HorizontalAlignment', 'left', ...
    'BackgroundColor','white', ...
    'Value', 1, ...
    'FontSize',26, 'Callback',@SoundCompareGUI_StopSliderCallback);

% Clock
data.clockLabel = uicontrol('Parent',data.mainPanel,'Style','text', 'String','00', 'Units','normalized', ...
    'Position',[0.05 .05 0.2 0.1], 'HorizontalAlignment', 'center', ...
    'BackgroundColor','white', ...
    'FontSize',48);


% Progress Bar
data.axis = axes('Parent',data.mainPanel,'Position',[.3 0.05 0.6 0.08],'XTick',[],'YTick',[],'YLim',[0 1],'XLim',[0 1]);
data.progressBarBG = rectangle(data.axis,'Position',[0 0 1 1],'EdgeColor','k','FaceColor',[0.8 0.8 0.8]);
data.progressBar = rectangle(data.axis,'Position',[0 0 0.5 1],'EdgeColor','k','FaceColor',[0 0.7 0]);
data.startLine = line(data.axis,[0 0],[0 1],'Color','k','LineWidth',3);
data.endLine = line(data.axis,[1 1],[0 1],'Color','k','LineWidth',3);


% Set data initial values
data.PLAYING = false;
data.MODIFICATION = false;
data.TRANSITION = false;
data.audioPlayer = [];
data.CURRENT = 1;
data.START = 0;
data.END = 1;

data.updateTimer = timer('ExecutionMode', 'fixedRate', 'Period', 0.25, 'TimerFcn', @SoundCompareGUI_TimerUpdate);

% Write data
guidata(gcf, data);
return;

end % END function createGUI






%% Sound GUI Close Function
%=========================================================================
function SoundCompareGUI_CloseFcn(object, eventdata)
data = guidata(gcf);

% Stop playing files
if (~isempty(data.audioPlayer))
    stop(data.audioPlayer);
end

stop(data.updateTimer);

delete(object);
end








%% Sound GUI Start Button Callback
%=========================================================================
function SoundCompareGUI_SoundFileBtnCallback(object, eventdata)
data = guidata(gcf);
data.CURRENT = object.UserData;

% Stop playing files
if (~isempty(data.audioPlayer))
    stop(data.audioPlayer);
end

% Set current text
data.currentFileLabel.String = data.fileNames{data.CURRENT};

% Play audioplayer
data.audioPlayer = audioplayer(data.WAVs{data.CURRENT},data.FS(data.CURRENT));
len = length(data.WAVs{data.CURRENT});
startPnt = max(1,round(data.START*len));
endPnt = min(len,round(data.END*len));
play(data.audioPlayer, [startPnt endPnt]);

guidata(gcf,data);
end







%% Sound GUI Stop Button Callback
%=========================================================================
function SoundCompareGUI_StopBtnCallback(object, eventdata)
data = guidata(gcf);

% Stop playing files
if (~isempty(data.audioPlayer))
    stop(data.audioPlayer);
end

guidata(gcf,data);
end






%% Sound GUI Start Slider Callback
%=========================================================================
function SoundCompareGUI_StartSliderCallback(object, eventdata)
data = guidata(gcf);

if (object.Value >= data.END)
    object.Value = data.START;
    return;
   % object.Value = max(data.END - 1, 0);
end
data.START = object.Value;
data.startLine.XData = [data.START data.START];


% Stop playing files
if (~isempty(data.audioPlayer))
    stop(data.audioPlayer);
end

% Set current text
data.currentFileLabel.String = data.fileNames{data.CURRENT};

% Play audioplayer
data.audioPlayer = audioplayer(data.WAVs{data.CURRENT},data.FS(data.CURRENT));
len = length(data.WAVs{data.CURRENT});
startPnt = max(1,round(data.START*len));
endPnt = min(len,round(data.END*len));
play(data.audioPlayer, [startPnt endPnt]);

guidata(gcf,data);
end




%% Sound GUI Stop Slider Callback
%=========================================================================
function SoundCompareGUI_StopSliderCallback(object, eventdata)
data = guidata(gcf);

if (object.Value <= data.START)
    object.Value = data.END;
    return;
   % object.Value = min(data.START + 1, 1);
end
data.END = object.Value;
data.endLine.XData = [data.END data.END];


% Stop playing files
if (~isempty(data.audioPlayer))
    stop(data.audioPlayer);
end

% Set current text
data.currentFileLabel.String = data.fileNames{data.CURRENT};

% Play audioplayer
data.audioPlayer = audioplayer(data.WAVs{data.CURRENT},data.FS(data.CURRENT));
len = length(data.WAVs{data.CURRENT});
startPnt = max(1,round(data.START*len));
endPnt = min(len,round(data.END*len));
play(data.audioPlayer, [startPnt endPnt]);

guidata(gcf,data);
end




%% Sound GUI Timer Callback
%=========================================================================
function SoundCompareGUI_TimerUpdate(object, eventdata)
data = guidata(gcf);

if (~isempty(data.audioPlayer))
    
    % Update Clock
    data.clockLabel.String = num2str(round(data.audioPlayer.currentSample / data.FS(data.CURRENT)));
    
    % Update Progress Bar
    wth = data.audioPlayer.currentSample / length(data.WAVs{data.CURRENT});
    data.progressBar.Position = [data.progressBar.Position(1) data.progressBar.Position(2) wth data.progressBar.Position(4)];
    drawnow;
end
guidata(gcf,data);
end










































