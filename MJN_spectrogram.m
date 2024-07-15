%---------------------------------------------------------------------------------------------%
% FUNCTION NAME AND BASIC SPECIFICATION
%
%   [specData] = ...
%           MJN_spectrogram(params,audioData)
%---------------------------------------------------------------------------------------------%
% Author:           Dr Mike Newton
% Date:             May 2024
% Location (local): [Matlab_root]/LIBRARY/MJN_Code_Library/MJN_DSP__git/
% GitHub location:  https://github.com/self-noise/Audio-Acoustics-DSP-Tools
%---------------------------------------------------------------------------------------------%
% PURPOSE OF THIS FUNCTION:
%           (1) Produces a spectrogram plot of a discrete time signal(s) that is passed to it 
%           as a vector or matrix of numbers
%---------------------------------------------------------------------------------------------%
% INPUTS:
%           params                      Structure containing all control parameters
%               .Fs                     Sample rate of audio signals
%               .windowSize_ms          The length or type of the analysis frame.
%                                       Either a number specified in ms ([10-500] ms
%                                       is typical - though note that in actual
%                                       implementation some rounding due to sample
%                                       rate will usually happen), or the text flag
%                                       'wholeSignal', in which case no windowing will
%                                       happen and a single spectrum will be returned
%                                       for each input channel
%               .params.windowType      ['Hann' or 'Rect']
%                                       Type of window applied to each
%                                       analysis frame
%               .overlapFraction        [0:1]
%                                       The fraction of the current window
%                                       that will be included in the next
%                                       window (0.5 means that the next
%                                       window will include 50% of the
%                                       current one; 0.2 means that 20%
%                                       will be included; 1 means that no
%                                       sliding windows will be used, and
%                                       hence the NSC of the whole file
%                                       will be returned [winSize_Analysis
%                                       is then ignored]; 0 means that
%                                       consecutive analysis windows do not
%                                       overlap at all)
%
%           audioData:  time domain signal (can be single or multichannel,
%                       ordered as columns in a single matrix)
%
% OUTPUTS:
%           specData (spectrogram as a matrix of frames)
%
%---------------------------------------------------------------------------------------------%
% GENERAL USAGE NOTES:
%   NOTE 1:
%---------------------------------------------------------------------------------------------%
% CHANGES TO ADD AT SOME POINT IN THE FUTURE:
%   "specData" changes size on each loop, adding a new entry to the structure for each channel in
%   sequence. This should be optimised/pre-allocated.
%---------------------------------------------------------------------------------------------%
% CHANGELOG:
%   2024-07-15: Updated filename for adding to GitHub
%   2024-05-09: Initial coding
%---------------------------------------------------------------------------------------------%
function [specData] = MJN_spectrogram(params,audioData)

N_channels      = size(audioData,2);
Ts              = 1/params.Fs;
N_tot           = length(audioData);
T_end           = N_tot/params.Fs;
tVec            = [0:Ts:T_end-Ts];


% The most common case, in which a window size is given in ms and converted to a
% number of samples
N_frame         = round((params.windowSize_ms/1000)*params.Fs);

% Hop size (samples) and number of frames needed (samples)
HA              = N_frame-round(N_frame*params.overlapFraction);
NF              = ceil(N_tot/HA);

% Post-pad with zeros so that our final frame contains at least some non-zero data
N_tot_pad       = (NF*HA)+N_frame;

% Make a time axis for the frames to be windowed: Each entry here is the
% timestamp of the 'centre' of the corresponding windowed frame
tVec_FrameCentres = 1:HA:NF*HA;

audioData = [audioData;zeros(N_tot_pad-N_tot,N_channels)];

% Define a windowing function (all ones is just a square window)
switch params.windowType
    case 'Hann'
        win                  = 0.5*(1 - cos(2*pi*(0:N_frame-1)'/N_frame)); % Manual computation of Hann window
    case 'Rect'
        win                  = ones(N_frame,1); % All ones
    otherwise
        error('Error: You have not specified a valid window type')
end

% Derived frequency domain parameters for each slice (full bandwidth as
% determined by sample rate and frame length)
df                          = params.Fs/N_frame;
fVec                        = [0:df:params.Fs-df]';

% Big empty matrices to store all the frames (the 'raw' version is used
% to compute the RMS level only, for use in estimating the signal
% level during each frame)
audioData_slices_raw       = zeros(N_frame,NF);
audioData_slices           = zeros(N_frame,NF);

% Loop over the number of channels (i.e., microphones used)
for nChan = 1:N_channels

    % Read the frames into the big matrix of slices
    for m=0:NF-1
        indRange                    = (m*HA)+1:(m*HA)+N_frame;
        audioData_slices_raw(:,m+1) = audioData(indRange,nChan);
        audioData_slices(:,m+1)     = audioData(indRange,nChan).*win;
    end

    % Compute the RMS level of each slice (RAW slice, NOT including
    % windowing)
    % Note: If this function has been given pre-calibrated acoustic
    % signals (i.e., units of Pa) then this will be the RMS sound
    % pressure in units of Pa for the given slice
    specData(nChan).sigLevel_RMS = rms(audioData_slices_raw);

    size(audioData_slices)

    % Take the DFT of all the slices (using the normed version of the DFT)
    this_DFT = fft(audioData_slices)./sqrt(N_frame);

    % Depending on whether DFT length (N_frame) is even or odd, discard upper
    % half of spectrum (assuming input time signal is real) to save memory and
    % simplify things
    even_check = rem(N_frame, 2);
    if even_check == 0
        % N_frame is even
        this_DFT = this_DFT(1:1+(N_frame/2),:);
        fVec = fVec(1:(1+(N_frame/2)));
    elseif even_check == 1
        % N_frame is odd
        this_DFT = this_DFT(1:(N_frame+1)/2,:);
        fVec = fVec(1:(N_frame+1)/2);
    end

    specData(nChan).spectrogram  = this_DFT;

    % Make a copy of the time vector (even if this creates duplicates)
    specData(nChan).tVec_FrameCentres  = [((N_frame/2)+tVec_FrameCentres)/params.Fs];

    % Make a copy of the frequency vector (even if this creates duplicates)
    specData(nChan).fVec = fVec;
end

