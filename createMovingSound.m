function [WAV] = createMovingSound(sound, fs, IRs, instructions, varargin)
% function [WAV] = createMovingSound(sound, fs, IRs, instructions, varargin)
%
% Generate a moving sound source simulation with the given IR set.
%
% INPUTS:
%    sound = Sound waveform (samples x 1)
%    fs = Sampling frequency (Default: 44100)
%    IRs = (string path) Set of azimuth HRTFs/IRs for each position to load.
%    instructions = Cell array of locations for motion instructions.
%
%  Instruction parameters:
%       {'station', Azimuth, Distance, TimeToStay}
%       {'move', X, Y, TimeOfMovement}
%       {'rotate', R, TimeOfRotation}
%       {'hold', TimeOfHold}
%     If Time is greater than length of input waveform, set to end of waveform.
%     Must start with a 'station' instruction as starting point.
%
% Optional parameters:
%      'bufferFreq', 'N' - Times to update position (Hz) [Default: 10 Hz]
%      'overlap', 'N' - Time of overlap of two locations (msec). [Default: 10 msec]
%      'roomPath', 'S' - Path where MCRoomSim is installed.
%      'room' - {[X,Y,Z], [rX,rY,rZ]} - Set room size [X,Y,Z] and receiver position [rX,rY,rZ].
%           Room Dimensions: length (x), width (y), height (z)
%           Receiver Position within the room.  +X is forward, +Y is leftward, +Z is upward
%
% OUTPUTS:
%    WAV = the processed output waveform
%
%   REQUIRES SIGNAL PROCESSING TOOLBOX
%
% @Author: Eugene Brandewie - Sep 23th, 2021



%% INPUT HANDLING
%=========================================
timeStart = clock;

if (nargin < 3)
    error('Requires a waveform, fs, and IRs.');
end

if (size(sound,2) > 1)
    error('Too many columns. Can only process one waveform channel.');
end

% Load HARTFs
load([IRs]);

% MCRoomSim Path
mcRoomSimPath = './MCRoomSim-master/';



% OPTIONAL PARAMETERS
%=========================================
% Defaults:
bufferFreq = 10; % Times / Second (Hz)
overlap =  10; % Overlap time (mSec)



% Room Parameter Defaults
%----------------------------
ROOM_SIM = true;
% Room Dimensions: length (x), width (y), height (z)
roomDims = [5, 5.5, 3];
% Receiver Position within the room.  +X is forward, +Y is leftward, +Z is upward
receiverPosition = [2 3 1.2]; % X, Y, Z
roomParams = [];


IS_PARSED = false;
for bb = (1:length(varargin))
    if (IS_PARSED), IS_PARSED = false; continue; end
    str = varargin{bb};
    if (~isnumeric(str))
        
        % BUFFER FREQUENCY
        if strcmpi(str,'bufferFreq')
            bufferFreq = varargin{bb+1};
            IS_PARSED = true;
        end
        
        % OVERLAP TIME
        if strcmpi(str,'overlap')
            overlap = varargin{bb+1};
            IS_PARSED = true;
        end
        
        % ROOM
        if strcmpi(str,'room')
            tempPar = varargin{bb+1};
            % IF contains room sim data
            if (iscell(tempPar))
                roomDims = tempPar{1};
                if (length(roomDims) ~= 3)
                    error('Room Dimensions must be [X Y Z]');
                end
                
                receiverPosition = tempPar{2};
                if (length(receiverPosition) ~= 3)
                    error('Receiver Position must be [X Y Z]');
                end
            else %ELSE turning off room simulation
                if (~tempPar)
                    ROOM_SIM = false;
                end
            end
            IS_PARSED = true;
        end
        
        
        % ROOM SIM PATH
        if strcmpi(str,'roomPath')
            mcRoomSimPath = varargin{bb+1};
            IS_PARSED = true;
        end
        
    end
end


% Load in MCRoomSim Path
addpath([mcRoomSimPath]);



% ALLOCATE BUFFERS
%=====================================================
disp('   Processing Moving Source Simulation');
disp('--------------------------------------------');
disp(['      Update Frequency: ' num2str(bufferFreq) ' Hz']);
disp(['      Overlap: ' num2str(overlap) ' mSec']);

% Get samples per buffer
bufferSamps = round((1 / bufferFreq) * fs);
if (bufferSamps < 1), bufferSamps = 1; warning('Buffer Frequency too low, set to 1 sample.'); end
overlapSamps = overlap / 1000 * fs;

% Get number of buffers for stimuli
numberBuffers = ceil(length(sound)/bufferSamps);
starts = zeros(numberBuffers,1);
ends = zeros(numberBuffers,1);
% Process stimuli into separate buffers
bufferData = zeros(overlapSamps + bufferSamps + overlapSamps, numberBuffers);

% Modify input waveform for overlap buffers
SOUND = [zeros(overlapSamps,1); sound; zeros(overlapSamps,1)];

% FOR each buffer
for bb = 1:numberBuffers
    % Get start and end samples for this time chunk
    if (bb == 1)
        starts(bb) = 1;
    else
        starts(bb) = (bufferSamps * (bb-1)) - overlapSamps;
    end
    ends(bb) = starts(bb) + bufferSamps + overlapSamps + overlapSamps;
    % IF too long
    if (ends(bb) > length(SOUND))
        Bleft = ends(bb) - length(SOUND);
        bufferData(:,bb) = cosSquareRamps([SOUND(starts(bb):end,1); zeros(Bleft,1)],fs,overlap);
    else
        bufferData(:,bb) = cosSquareRamps(SOUND(starts(bb):ends(bb),1),fs,overlap);
    end
end



% PROCESS INSTRUCTIONS
%=====================================================

% Set parameters
currentBuffer = 1; % Current buffer being processed
WAV = zeros(ends(end) + (2*fs),2);  % Final waveform
currentAz = 0; % Current Azimuth position (deg)
currentDist = 2.0;   % Current distance position (m)


% FOR each instruction
for (bb = 1:length(instructions))
    currentInstruction = instructions{bb};
    % Switch by Type
    instructionType = currentInstruction{1};
    
    % STATION
    %----------------------------
    % Sets source as stationary at the given location.
    % Format: {'station', Azimuth, Distance, TimeToStay}
    if (strcmpi(instructionType,'station'))
        azimuth = currentInstruction{2};
        distance = currentInstruction{3};
        sampsToStay = currentInstruction{4} * fs; % convert to samples
        
        disp(['   Processsing: ' currentInstruction{1} ' : ' num2str(azimuth) ' o : ' num2str(distance) ' m : ' num2str(currentInstruction{4}) ' s']);
        
        % Calculate buffers
        sectionBuffers = ceil(sampsToStay / bufferSamps);
        if (sectionBuffers+currentBuffer > numberBuffers)
            sectionBuffers = numberBuffers-currentBuffer;
        end
        
        % Calculate required IR for this position
        [sourceIR] = calculateIR(azimuth, distance, roomDims, receiverPosition, roomParams, hrirMieLeft, hrirMieRight, fsHrir, directions);
        
        
        % FOR each section buffer
        for (bb = 0:sectionBuffers-1)
           % disp(['   .   Processing Buffer: ' num2str(bb+1) ' / ' num2str(sectionBuffers)]);
            
            % Convolve section with IR
            sectionWav = freqConvolve(bufferData(:,currentBuffer+bb),sourceIR);
            
            % Add to final output
            WAV(starts(currentBuffer+bb) : ends(currentBuffer+bb)+length(sourceIR)-2,:) = WAV(starts(currentBuffer+bb) : ends(currentBuffer+bb)+length(sourceIR)-2,:) + sectionWav;
            
        end
        
        % Update current positioning
        currentAz = azimuth;
        currentDist = distance;
        currentBuffer = currentBuffer + sectionBuffers;
        
        
        
    % HOLD
    %----------------------------
    % Sets source as stationary at the current location.
    % Format: {'hold', TimeToHold}
    elseif (strcmpi(instructionType,'hold'))
        sampsToHold = currentInstruction{2} * fs; % convert to samples
        
        disp(['   Processsing: ' currentInstruction{1} ' : ' num2str(currentInstruction{2}) ' s']);
        
        % Calculate buffers
        sectionBuffers = ceil(sampsToHold / bufferSamps);
        if (sectionBuffers+currentBuffer > numberBuffers)
            sectionBuffers = numberBuffers-currentBuffer;
        end
        
        % Calculate required IR for this position
        azimuth = currentAz;
        distance = currentDist;
        [sourceIR] = calculateIR(azimuth, distance, roomDims, receiverPosition, roomParams, hrirMieLeft, hrirMieRight, fsHrir, directions);
        
        
        % FOR each section buffer
        for (bb = 0:sectionBuffers-1)
           % disp(['   .   Processing Buffer: ' num2str(bb+1) ' / ' num2str(sectionBuffers)]);
            
            % Convolve section with IR
            sectionWav = freqConvolve(bufferData(:,currentBuffer+bb),sourceIR);
            
            % Add to final output
            WAV(starts(currentBuffer+bb) : ends(currentBuffer+bb)+length(sourceIR)-2,:) = WAV(starts(currentBuffer+bb) : ends(currentBuffer+bb)+length(sourceIR)-2,:) + sectionWav;
            
        end
        
        % Update current positioning
        currentAz = azimuth;
        currentDist = distance;
        currentBuffer = currentBuffer + sectionBuffers;
        
        
        
        
        
        % ROTATE
        %----------------------------
        % Sets source as rotating by the given amount (neg is left, pos is right).
        % Format: {'station', AzimuthChange, TimeOfRotation}
    elseif (strcmpi(instructionType,'rotate'))
        azimuthDelta = currentInstruction{2};
        sampsToRotate = currentInstruction{3} * fs; % convert to samples
        
        disp(['   Processsing: ' currentInstruction{1} ' : ' num2str(azimuthDelta) ' o : ' num2str(currentDist) ' m : ' num2str(currentInstruction{3}) ' s']);
        
        % Calculate buffers
        sectionBuffers = ceil(sampsToRotate / bufferSamps);
        if (sectionBuffers+currentBuffer > numberBuffers)
            sectionBuffers = numberBuffers-currentBuffer;
        end
        
        % Set parameters
        rotationPerBuffer = azimuthDelta / sectionBuffers;
        rotAzimuth = currentAz;
        distance = currentDist;
        
        % FOR each section buffer
        for (bb = 0:sectionBuffers-1)
            
            disp(['   .   Processing Buffer: ' num2str(bb+1) ' / ' num2str(sectionBuffers)]);
            
            % Calculate required IR for this position
            [sourceIR] = calculateIR(rotAzimuth, distance, roomDims, receiverPosition, roomParams, hrirMieLeft, hrirMieRight, fsHrir, directions);
            
            % Convolve section with IR
            sectionWav = freqConvolve(bufferData(:,currentBuffer+bb),sourceIR);
            
            % Add to final output
            WAV(starts(currentBuffer+bb) : ends(currentBuffer+bb)+length(sourceIR)-2,:) = WAV(starts(currentBuffer+bb) : ends(currentBuffer+bb)+length(sourceIR)-2,:) + sectionWav;
            
            % Update current azimuth
            rotAzimuth = rotAzimuth + rotationPerBuffer;
        end
        % Update current positioning
        currentAz = rotAzimuth;
        currentDist = distance;
        currentBuffer = currentBuffer + sectionBuffers;
        
        
        
        
        % MOVE
        %----------------------------
        % Sets source as moving linearly by the given X/Y changes.
        % Format: {'station', deltaX, deltaY, TimeOfMovement}
    elseif (strcmpi(instructionType,'move'))
        deltaX = currentInstruction{2};
        deltaY = currentInstruction{3};
        sampsToMove = currentInstruction{4} * fs; % convert to samples
        
        disp(['   Processsing: ' currentInstruction{1} ' : ' num2str(deltaX) ' m : ' num2str(deltaY) ' m : ' num2str(currentInstruction{4}) ' s']);
        
        % Calculate buffers
        sectionBuffers = ceil(sampsToMove / bufferSamps);
        if (sectionBuffers+currentBuffer > numberBuffers)
            sectionBuffers = numberBuffers-currentBuffer;
        end
        
        % Set parameters
        deltaXPerBuffer = deltaX / sectionBuffers;
        deltaYPerBuffer = deltaY / sectionBuffers;
        azimuth = currentAz;
        distance = currentDist;
        [X,Y] = pol2cart(deg2rad(-azimuth),distance);
        
        % FOR each section buffer
        for (bb = 0:sectionBuffers-1)
            
            disp(['   .   Processing Buffer: ' num2str(bb+1) ' / ' num2str(sectionBuffers)]);
            
            % Calculate required IR for this position
            [sourceIR] = calculateIR(azimuth, distance, roomDims, receiverPosition, roomParams, hrirMieLeft, hrirMieRight, fsHrir, directions);
            
            % Convolve section with IR
            sectionWav = freqConvolve(bufferData(:,currentBuffer+bb),sourceIR);
            
            % Add to final output
            WAV(starts(currentBuffer+bb) : ends(currentBuffer+bb)+length(sourceIR)-2,:) = WAV(starts(currentBuffer+bb) : ends(currentBuffer+bb)+length(sourceIR)-2,:) + sectionWav;
            
            % Update current azimuth/distance
            [X,Y] = pol2cart(deg2rad(-azimuth),distance);
            X = X + deltaXPerBuffer;
            Y = Y + deltaYPerBuffer;
            [azimuth,distance] = cart2pol(X,Y);
            azimuth = -rad2deg(azimuth);
        end
        % Update current positioning
        currentAz = azimuth;
        currentDist = distance;
        currentBuffer = currentBuffer + sectionBuffers;
        
    end
end % END for each instruction





% PROCESS FINAL POSITION UNTIL END
%=======================================
if (currentBuffer <= numberBuffers)
    
    disp(['   Processsing: stay : ' num2str(currentAz) ' o : ' num2str(currentDist) ' m']);
    
    % Calculate buffers
    sectionStart = starts(currentBuffer);
    sectionEnd = ends(end);
    
    % Calculate required IR for this position
    [sourceIR] = calculateIR(currentAz, currentDist, roomDims, receiverPosition, roomParams, hrirMieLeft, hrirMieRight, fsHrir, directions);
    
    % Convolve section with IR
    sectionWav = freqConvolve(cosSquareRamps(sound(sectionStart:end),fs,10),sourceIR);
    
    % Add to final output
    WAV(sectionStart:sectionStart + length(sectionWav) - 1,:) = WAV(sectionStart:sectionStart + length(sectionWav) - 1,:) + sectionWav;
end


disp(['   Processing Finished:  ' num2str(etime(clock,timeStart)) ' s']);
disp('--------------------------------------------');






end % END function