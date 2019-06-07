function soundTrimmer(readDirectory, writeDirectory, varargin)
% function soundTrimmer(readDirectory, writeDirectory)
%
%  Program utility for processing and trimming a batch of sound (wav) files.
%    Uses a UI interface to click on the waveform to select the beginning
%    and end of each sound file for trimming. Also allows sound file to be
%    listened to with trimmed edits in place.
%
%    *Trimming Mode* - Two points are chosen at the beginning and end of the
%       sound file. Cosine squared gain ramps are applied to remove onset/offset transients.
%
%    *Marking Mode* - Can select a number of marking points that are then saved to a '.mat'
%       file for future use. These files are saved into the 'writeDirectory' instead of
%       the processed sound files. Nearest zero-crossing is found to reduce transients.
%
%
%   INPUT:
%     readDirectory   = Directory to read sound files.
%     writeDirectory  = Directory to write processed sound files.
%
%   OPTIONAL INPUTS:
%     'nbits', [N] =   [Default: 24] - Number of bits to use for writing trimmed audio files.

%     'ramps', [N N] = [Default: [30 30]] - [StartRamp EndRamp] Duration of cos^2 ramps in ms.
%                       Applied to start and end points after trimming.
%                       Not applicable to Marking Mode.
%
%     'zero', [N]    = [Default: 0.5] - Tolerance value (in ms) for adaptive zero-crossing,
%                       which finds the nearest zero-crossing from the selected points.
%                       Only optimizes first channel. Only used in Trimming Mode if ramps are disabled.
%
%  REQUIRES: Signal Processing Toolbox
%
%
% @Author: Eugene Brandewie
% @Version: 1.0 - Dec 21, 2005 - Intial working version.
%           2.0 - Jan 22, 2019 - Refactored and added optional features.
%


% PROCESS INPUTS
%=================================================

RAMPS = [30 30];
ZERO_TOLERANCE = 0.5;
MARKING_MODE = false;
NMARKS = 2;
NBITS = 24;

if nargin >= 2
    for ii = 1:nargin-2
        S = varargin{ii};
        if (~isnumeric(S))
            % RAMPS
            if strcmpi(S,'ramps')
                RAMPS = varargin{ii+1};
                if (length(RAMPS) ~= 2)
                    if (length(RAMPS) == 1)
                        RAMPS = [RAMPS RAMPS];
                    else
                        error('''ramps'' must have two terms [StartRamp EndRamp] in ms.');
                    end
                end
            end
            % ZERO TOLERANCE
            if strcmpi(S,'zero')
                ZERO_TOLERANCE = varargin{ii+1};
                if (ZERO_TOLERANCE < 0)
                    error('''zero'' must have tolerance (in ms).');
                end
            end
            % MARKING MODE
            if strcmpi(S,'marking')
                MARKING_MODE = true;
            end
            % AUDIO BITS
            if strcmpi(S,'nbits')
                NBITS = varargin{ii+1};
                if (NBITS ~= [8 16 24 32 64])
                    error(['''nbits'' must be [8, 16, 24, 32, or 64].']);
                end
            end
        end
    end
end


% INPUT ERROR HANDLING
%============================================
if (exist(readDirectory,'dir') ~= 7)
    error(['Read directory cannot be found: ' readDirectory]);
end

if (exist(writeDirectory,'dir') ~= 7)
    error(['Write directory cannot be found: ' writeDirectory]);
end



% GET LIST OF SOUND FILES TO PROCESS
%===========================================
FILES = dir(readDirectory);
% Find only sound files
cnt = 0;
% FOR each file in directory
for (ff = 1:length(FILES))
    fName = FILES(ff).name;
    % IF is a wav file
    if (length(fName) > 4 && strcmpi(fName(end-2:end),'wav'))
        cnt = cnt + 1;
        % Add to list
        INPUT_FILES(cnt,1) = {fName};
    end
end
clear FILES cnt;



% LOOP PREPARATION
%==========================================

% Prepare UI
FIG = figure();
MODE = 'Play';

% Prepare file listing
currentNum = 1;
currentFile = INPUT_FILES{currentNum};

% Read in sound file
[WAV, FS] = audioread([readDirectory currentFile]);
TRIM_WAV = WAV;

% Prepare trims and marks
startPnt = 1;
endPnt = length(WAV);
marks = [];
currentMark = 1;
Ymin = min(min(WAV));
Ymax = max(max(WAV));
X_LIMITS = [0 length(WAV)];



% BEGIN PROCESSING LOOP
%==========================================

% WHILE endless
while(1)
    
    % Plot sound file to figure
    clf(FIG);
    hold off;
    plot(WAV);
    hold on;
    ylim([Ymin Ymax]);
    xlim(X_LIMITS);
    
    % Plot Trimmed Ends
    if (RAMPS(1) > 0)
        sOnset = round(RAMPS(1) * (FS/1000));
        rectangle('Position', [startPnt Ymin sOnset Ymax-Ymin],'EdgeColor','k','FaceColor','m');
    end
    if (RAMPS(2) > 0)
        sOffset = round(RAMPS(2) * (FS/1000));
        rectangle('Position', [endPnt-sOffset Ymin sOffset Ymax-Ymin],'EdgeColor','k','FaceColor','m');
    end
    line([startPnt startPnt],[Ymin Ymax],'Color','k','LineWidth',2);
    line([endPnt endPnt],[Ymin Ymax],'Color','k','LineWidth',2);
    
    % Plot Marks
    for (mm = 1:length(marks))
        line([marks(mm) marks(mm)],[Ymin Ymax],'Color','r','LineWidth',2);
    end
    
    
    
    % PLAYBACK MODE
    %================================
    if (strcmp(MODE,'Play'))
        
        title({['PLAYING: ' currentFile ' | File: ' num2str(currentNum) ' of ' num2str(length(INPUT_FILES))],...
            ['''P'' - Play entire sound file, ''T'' to Trim, ''M'' to Mark,'], ...
            [''' ] '' for Next,  '' [ '' for Previous''']}, 'Interpreter', 'none');
        
        
        % Get user input
        [~,~,B] = ginput(1);
        
        
        % 'P' - Play All
        %++++++++++++++++++++++++++++++++
        if (B == 112)
            soundsc(WAV,FS);
        end
        
        % 'T' - Enter Trim Mode
        %++++++++++++++++++++++++++++++++
        if (B == 116)
            MODE = 'Trim';
        end
        
        % 'M' - Enter Mark Mode
        %++++++++++++++++++++++++++++++++
        if (B == 109)
            MODE = 'Mark';
            marks = [];
            currentMark = 1;
        end
        
        % '[' - Prev File
        if (B == 91)
            % Prepare file listing
            currentNum = currentNum - 1;
            if (currentNum < 1), currentNum = length(INPUT_FILES); end
            currentFile = INPUT_FILES{currentNum};
            
            % Read in sound file
            [WAV, FS] = audioread([readDirectory currentFile]);
            
            % Prepare trims and marks
            startPnt = 1;
            endPnt = length(WAV);
            marks = [];
            currentMark = 1;
            Ymin = min(min(WAV));
            Ymax = max(max(WAV));
            X_LIMITS = [0 length(WAV)];
        end
        
        % ']' - Next File
        if (B == 93)
            % Prepare file listing
            currentNum = currentNum + 1;
            if (currentNum > length(INPUT_FILES)), currentNum = 1; end
            currentFile = INPUT_FILES{currentNum};
            
            % Read in sound file
            [WAV, FS] = audioread([readDirectory currentFile]);
            
            % Prepare trims and marks
            startPnt = 1;
            endPnt = length(WAV);
            marks = [];
            currentMark = 1;
            Ymin = min(min(WAV));
            Ymax = max(max(WAV));
            X_LIMITS = [0 length(WAV)];
        end
        
        
        % TRIMMING MODE
        %================================
    elseif (strcmp(MODE,'Trim'))
        
        title({['TRIMMING: ' currentFile ' | File: ' num2str(currentNum) ' of ' num2str(length(INPUT_FILES))],...
            ['''Space'' - Play trimmed portion,  L-Click to choose Start,  R-Click to choose End |'] ...
            ,['''X'' to return to playback menu,  ''S'' to Save'], ...
            ['''Z'' to Zoom in at point, ''R'' to Reset zoom'] ...
            }, 'Interpreter', 'none');
        
        % Get user input
        [X,~,B] = ginput(1);
        X = round(X);
        
        
        % 'X' - Cancel
        %++++++++++++++++++++++++++++++++
        if (B == 120)
            MODE = 'Play';
        end
        
        % 'R' - Reset Zoom
        %++++++++++++++++++++++++++++++++
        if (B == 114)
            X_LIMITS = [0 length(WAV)];
        end
        
        % 'Z' - Zoom In
        %++++++++++++++++++++++++++++++++
        if (B == 122)
            curX = get(gca,'XLim');
            spanX = round(0.5 * (curX(2) - curX(1))/2);
            X_LIMITS = [X-spanX X+spanX];
        end
        
        % 'S' - Save
        %++++++++++++++++++++++++++++++++
        if (B == 115)
            if (startPnt > endPnt)
                warning('Trimmed points may be reversed.');
            else
                % Create Gain Ramps
                if (RAMPS(1) > 0)
                    % Get sample lengths
                    sOnset = round(RAMPS(1) * (FS/1000));
                    onsetRamp = ((sin(pi*(3/2):pi/(sOnset-1):pi*(5/2))+1)/2).^2;
                else
                    onsetRamp = [];
                end
                if (RAMPS(2) > 0)
                    % Get sample lengths
                    sOffset = round(RAMPS(2) * (FS/1000));
                    offsetRamp = ((sin(pi*(3/2):pi/(sOffset-1):pi*(5/2))+1)/2).^2;
                else
                    offsetRamp = [];
                end
                % Apply ramps (if any)
                TRIM_WAV = WAV(startPnt:endPnt,:);
                if (~isempty(onsetRamp))
                    TRIM_WAV(1:length(onsetRamp),:) = TRIM_WAV(1:length(onsetRamp),:) .* onsetRamp';
                end
                if (~isempty(offsetRamp))
                    TRIM_WAV(end-length(offsetRamp)+1:end,:) = TRIM_WAV(end-length(offsetRamp)+1:end,:) .* offsetRamp';
                end
                
                audiowrite([writeDirectory currentFile],TRIM_WAV, FS,'BitsPerSample',NBITS);
                disp(['File written: ' currentFile]);
                MODE = 'Play';
            end
        end
        
        % 'Space' - Play (trimmed portion)
        %++++++++++++++++++++++++++++
        if (B == 32)
            
            % Create Gain Ramps
            if (RAMPS(1) > 0)
                % Get sample lengths
                sOnset = round(RAMPS(1) * (FS/1000));
                onsetRamp = ((sin(pi*(3/2):pi/(sOnset-1):pi*(5/2))+1)/2).^2;
            else
                onsetRamp = [];
            end
            if (RAMPS(2) > 0)
                % Get sample lengths
                sOffset = round(RAMPS(2) * (FS/1000));
                offsetRamp = ((sin(pi*(3/2):pi/(sOffset-1):pi*(5/2))+1)/2).^2;
            else
                offsetRamp = [];
            end
            % Apply ramps (if any)
            TRIM_WAV = WAV(startPnt:endPnt,:);
            if (~isempty(onsetRamp))
                TRIM_WAV(1:length(onsetRamp),:) = TRIM_WAV(1:length(onsetRamp),:) .* onsetRamp';
            end
            if (~isempty(offsetRamp))
                TRIM_WAV(end-length(offsetRamp)+1:end,:) = TRIM_WAV(end-length(offsetRamp)+1:end,:) .* offsetRamp';
            end
            soundsc(TRIM_WAV,FS);
        end
        
        
        % 'LClick' - Trimming - Start
        %++++++++++++++++++++++++++++++++
        if (B == 1) % Left-click
            % Zero Crossing Finder
            %---------------------------
            if (MARKING_MODE && ZERO_TOLERANCE > 0)
                % Get half-bandwidth for zero crossing
                zeroBandwidth = round(((ZERO_TOLERANCE/1000) * FS)/2);
                
                % Find nearest zero crossing (first channel only)
                subWav = (abs(WAV((X-zeroBandwidth):(X+zeroBandwidth),1)));
                indx = find(subWav == min(subWav));
                startPnt = (X-zeroBandwidth+indx-1);
            else
                % ELSE Take location as given
                startPnt = X;
            end
        end % End Left-Click
        
        
        % 'RClick' - Trimming - End
        %++++++++++++++++++++++++++++++++
        if (B == 3) % Right-click
            % Zero Crossing Finder
            %---------------------------
            if (MARKING_MODE && ZERO_TOLERANCE > 0)
                % Get half-bandwidth for zero crossing
                zeroBandwidth = round(((ZERO_TOLERANCE/1000) * FS)/2);
                
                % Find nearest zero crossing (first channel only)
                subWav = (abs(WAV((X-zeroBandwidth):(X+zeroBandwidth),1)));
                indx = find(subWav == min(subWav));
                endPnt = (X-zeroBandwidth+indx-1);
            else
                % ELSE Take location as given
                endPnt = X;
            end
        end % End Right-Click
        
        
        % MARKING MODE
        %================================
    elseif (strcmp(MODE,'Mark'))
        
        title({['MARKING: ' currentFile ' | File: ' num2str(currentNum) ' of ' num2str(length(INPUT_FILES))],...
            ['L-Click to mark point,  R-Click to return to previous mark, Mark:'] ...
            ,['''X'' to return to playback menu,  ''S'' to Save marks'], ...
            ['''Z'' to Zoom in at point, ''R'' to Reset zoom'] ...
            }, 'Interpreter', 'none');
        
        % Get user input
        [X,~,B] = ginput(1);
        X = round(X);
        
        
        
        % 'X' - Cancel
        %++++++++++++++++++++++++++++++++
        if (B == 120)
            MODE = 'Play';
        end
        
        % 'S' - Save
        %++++++++++++++++++++++++++++++++
        if (B == 115)
            save([writeDirectory currentFile(1:end-4) '.mat'],'marks');
            disp(['File written: ' [currentFile(1:end-4) '.mat'] ' | ' num2str(length(marks)) ' marks']);
            MODE = 'Play';
        end
        
        
        % 'R' - Reset Zoom
        %++++++++++++++++++++++++++++++++
        if (B == 114)
            X_LIMITS = [0 length(WAV)];
        end
        
        % 'Z' - Zoom In
        %++++++++++++++++++++++++++++++++
        if (B == 122)
            curX = get(gca,'XLim');
            spanX = round(0.66 * (curX(2) - curX(1))/2);
            X_LIMITS = [X-spanX X+spanX];
        end
        
        
        % 'Left-click' - Mark
        %++++++++++++++++++++++++++++++++
        if (B == 1) % Left-click
            
            % Zero Crossing Finder
            %---------------------------
            if (MARKING_MODE && ZERO_TOLERANCE > 0)
                % Get half-bandwidth for zero crossing
                zeroBandwidth = round(((ZERO_TOLERANCE/1000) * FS)/2);
                
                % Find nearest zero crossing (first channel only)
                subWav = (abs(WAV((X-zeroBandwidth):(X+zeroBandwidth),1)));
                indx = find(subWav == min(subWav));
                marks(currentMark) = (X-zeroBandwidth+indx-1);
            else
                % ELSE Take location as given
                marks(currentMark) = X;
            end
            currentMark = currentMark + 1;
        end % End Left-Click
        
        % 'Right-click' - Mark
        %++++++++++++++++++++++++++++++++
        if (B == 3) % Right-click
            if (currentMark > 1)
                marks = marks(1:end-1);
                currentMark = currentMark - 1;
            end
        end % End Right-Click
        
        
        
    end % END mode checking
    
    
end % END endless while loop



