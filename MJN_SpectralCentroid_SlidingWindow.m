%---------------------------------------------------------------------------------------------%
% FUNCTION NAME AND SPECIFICATION
%
%   [SC_normalised_frames,SC_frames,frame_RMS_levels,tVec_FrameCentres] = ...
%           MJN_SpectralCentroid_SlidingWindow(params,audioData)
%---------------------------------------------------------------------------------------------%
% PURPOSE OF THIS FUNCTION:
%           (1) Computes the normalised spectral centroid (NSC) of a time domain signal
%               which the function first windows into a series of frames, which may (or not)
%               be overlapping. All expected parameters can be controlled.
%               Expectation is that the signal is audio, but can be anything really (see
%               input parameter that controls calibration).
%               Note that by use of appropriate input parameters (see below), you can simply
%               return the NSC for the whole file (per channel).
%---------------------------------------------------------------------------------------------%
% INPUTS:
%           params                      Structure containing all control parameters
%               .Fs                     Sample rate of audio signals
%               .params.f_expected      The reference frequency for NSC
%                                       calculations (i.e., expected
%                                       playing frequency for the given
%                                       note) in Hz
%               .windowSize_ms          The length or type of the analysis frame.
%                                       Either a number specified in ms ([10-500] ms
%                                       is typical - though note that in actual
%                                       implementation some rounding due to sample
%                                       rate will usually happen), or the text flag
%                                       'wholeSignal', in which case no windowing will
%                                       happen and a single SC and NSC will be returned
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
%               .SC_lower_bound_Hz      Lower limit of frequency range used
%                                       when evaluating the NSC (in Hz)
%               .SC_upper_bound_Hz      Upper limit of frequency range used
%                                       when evaluating the NSC (in Hz)
%
%
%           audioData                   Time domain signal (can be single or multichannel, 
%                                       ordered as columns in a single matrix)
%
% OUTPUTS:
%           NSC_frames                  Normalised spectral centroid (by f_expected)
%           SC_frames                   Spectral centroid
%           frame_RMS_levels            RMS levels of each frame
%           tVec_FrameCentres           Time indices of each frame centre
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
%   2024-07-17: Updated preamble for GitHub
%
%   2024-04-22: Edited to make genuinely multichannel compatible (i.e., more than 2 channels),
%               in a way that is dynamic with the 'audioData' sent to this function
%
%   2024-02-01: Initial coding
%---------------------------------------------------------------------------------------------%
% Author:           Dr Mike Newton
% Date:             July 2024
% Location (local): [Matlab_root]/LIBRARY/MJN_Code_Library/MJN_DSP_git/
% GitHub location:  https://github.com/self-noise/Audio-Acoustics-DSP-Tools
%---------------------------------------------------------------------------------------------%
function [SC_normalised_frames,SC_frames,frame_RMS_levels,tVec_FrameCentres] = MJN_SpectralCentroid_SlidingWindow(params,audioData)

N_channels                  = size(audioData,2);
Ts__thisFile                = 1/params.Fs;
N_tot__thisFile             = length(audioData);
T_end__thisFile             = N_tot__thisFile/params.Fs;
tVec__thisFile              = [0:Ts__thisFile:T_end__thisFile-Ts__thisFile];


%---------------------------------------------------------------------------------------------%
% WORK OUT WHETHER WE NEED SLIDING WINDOWS OR NOT
%---------------------------------------------------------------------------------------------%
if params.overlapFraction > 1
    error(['Error: You have specified an overlap fraction of ' num2str(params.overlapFraction) ' which is not permitted'])

elseif params.overlapFraction < 0
    error(['Error: You have specified an overlap fraction of ' num2str(params.overlapFraction) ' which is not permitted'])

elseif params.overlapFraction == 1
    % Overlap fraction is exactly '1', hence no sliding windows are needed
    % Instead, we compute the NSC of the whole file, per channel, and
    % including whatever window spec has been requested

    % UPDATE (09 May 2024):
    % This part needs to be reconciled with the 'if' statement later on that actually implements a
    % non-sliding window arrangement based on how the user specifies "windowSize_ms". At the moment
    % there is duplication and not as smooth/clear a way to implement this decision tree as
    % desireable.

elseif params.overlapFraction == 0
    % ******* May be able to include this in >0 case        *******
    % Overlap fraction is exactly 0, hence we need sliding windows with no
    % overlap at all
    % ******* May be able to include this in >0 case below  *******


elseif params.overlapFraction > 0
    % Overlap fraction is in the range of 0:1 (more than zero, less than 1)
    % hence we need sliding windows --> This is the most common use case

    %----------------------------------------------------------------------------%
    % FILE WINDOWING AND SLICING AND PER-FRAME DERIVED PARAMETERS
    %----------------------------------------------------------------------------%

    % Window size (samples)
    if isnumeric(params.windowSize_ms)
        % The most common case, in which a window size is given in ms and converted to a
        % number of samples
        N_frame                     = round((params.windowSize_ms/1000)*params.Fs);

        % Hop size (samples) and number of frames needed (samples)
        HA                          = N_frame-round(N_frame*params.overlapFraction);
        NF                          = ceil(N_tot__thisFile/HA);

        % Post-pad with zeros so that our final frame contains at least some non-zero data
        N_tot__thisFile_pad         = (NF*HA)+N_frame;

        % Make a time axis for the frames to be windowed: Each entry here is the
        % timestamp of the 'centre' of the corresponding windowed frame
        tVec_FrameCentres = 1:HA:NF*HA;
        tVec_FrameCentres = [((N_frame/2)+tVec_FrameCentres)/params.Fs]';

    elseif params.windowSize_ms == 'wholeSignal'
        % No windowing will be done, and instead a single overall SC (and NSC) will be
        % returned for each channel - hence we just set the frame size to equal the total number
        % of samples
        N_frame = N_tot__thisFile;

        % Hop size (samples) and number of frames needed (samples)
        HA                          = 0;
        NF                          = 1;

        % Post-pad with zeros so that our final frame contains at least some non-zero data
        N_tot__thisFile_pad         = N_frame;

        % Make a time axis for the frames to be windowed: Each entry here is the
        % timestamp of the 'centre' of the corresponding windowed frame
        %tVec_FrameCentres = 1:HA:NF*HA;
        tVec_FrameCentres = N_frame/(2*params.Fs);        

    else
        error(['Error: You have not correctly specified the window size for the spectral ...' ...
            'analysis (should be either a number [in ms], or the text flag "wholeSignal".'])
    end

    audioData = [audioData;zeros(N_tot__thisFile_pad-N_tot__thisFile,N_channels)];

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

    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Derived Spectral Centroid parameters used in examining each frame %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Frequency (Hz) to bin indices: Convert the user's lower/upper bounds
    % in Hz into a bin index (note that rounding errors will occur)
    [~,SC_lower_bin] = min(abs(fVec - params.SC_lower_bound_Hz));
    [~,SC_upper_bin] = min(abs(fVec - params.SC_upper_bound_Hz));

    % Big empty matrices to store all the frames (the 'raw' version is used
    % to compute the RMS level only, for use in estimating the signal
    % level during each frame)
    audioData_slices_raw       = zeros(N_frame,NF);
    audioData_slices           = zeros(N_frame,NF);

    % Loop over the number of channels (i.e., microphones used)
    for nChan = 1:N_channels

        % Read the frames into the big matrix of slices
        for m=0:NF-1
            indRange = (m*HA)+1:(m*HA)+N_frame;
            audioData_slices_raw(:,m+1) = audioData(indRange,nChan);
            audioData_slices(:,m+1) = audioData(indRange,nChan).*win;
        end
        
        % Compute the RMS level of each slice (RAW slice, NOT including
        % windowing)
        % Note: If this function has been given pre-calibrated acoustic
        % signals (i.e., units of Pa) then this will be the RMS sound
        % pressure in units of Pa for the given slice
        audioData_slices_Level_RMS = rms(audioData_slices_raw);

        % Take the DFT of all the slices (using the normed version of the DFT)
        audioData_slices_DFT       = fft(audioData_slices)./sqrt(N_frame);
        audioData_slices_DFT_mag   = abs(audioData_slices_DFT);

        % Depending on whether DFT length (N_frame) is even or odd, discard upper
        % half of spectrum (assuming input time signal is real) to save memory and
        % simplify things
        even_check = rem(N_frame, 2);
        if even_check == 0
            % N_frame is even
            audioData_slices_DFT_mag_halfSpec = audioData_slices_DFT_mag(1:1+(N_frame/2),:);
            fVec__halfSpec = fVec(1:(1+(N_frame/2)));
        elseif even_check == 1
            % N_frame is odd
            audioData_slices_DFT_mag_halfSpec = audioData_slices_DFT_mag(1:(N_frame+1)/2,:);
            fVec__halfSpec = fVec(1:(N_frame+1)/2);
        end

        % Apply the lower/upper bounds to the frequency vector, setting all bins
        % below the lower limit to zero, and all bins above the upper limit to zero
        fVec__halfSpec_bounds = fVec__halfSpec;
        fVec__halfSpec_bounds(1:SC_lower_bin) = 0;
        fVec__halfSpec_bounds(SC_upper_bin+1:end) = 0;
        
        % Spectral Centroid: MY METHOD (for each frame via a vector*matrix multiplication)
        these_SC_frames                 = (fVec__halfSpec_bounds'*audioData_slices_DFT_mag_halfSpec)./sum(audioData_slices_DFT_mag_halfSpec,1);
        these_SC_normalised_frames      = these_SC_frames./params.f_expected;

        % Clean up the latter parts of the file (if needed), chosing the first
        % channel (i.e., mic) as the reference point to ensure that all
        % returned data (SC_frames, in particular) have the same length
        % HOWEVER, only both to do this if we are NOT using 'wholeSignal' for the frame length
        % (regardless of whether 'params.optionCleanTail' is 0 or 1)
        if isnumeric(params.windowSize_ms)
            switch params.optionCleanTail
                case 1
                    if nChan ==1
                        [~,indexSoundLevelMax]                                = max(audioData_slices_Level_RMS);
                        these_SC_frames(indexSoundLevelMax:end)               = [];
                        these_SC_normalised_frames(indexSoundLevelMax:end)    = [];
                        audioData_slices_Level_RMS(indexSoundLevelMax:end)    = [];
                        tVec_FrameCentres(indexSoundLevelMax:end)             = [];
                    else
                        these_SC_frames(indexSoundLevelMax:end)               = [];
                        these_SC_normalised_frames(indexSoundLevelMax:end)    = [];
                        audioData_slices_Level_RMS(indexSoundLevelMax:end)    = [];
                        tVec_FrameCentres(indexSoundLevelMax:end)             = [];
                    end
            end
        end
        % Write the channel-by-channel final data
        SC_frames(:,nChan)              = these_SC_frames';
        SC_normalised_frames(:,nChan)   = these_SC_normalised_frames';
        frame_RMS_levels(:,nChan)       = audioData_slices_Level_RMS;

        clear audioData_slices_raw audioData_slices audioData_slices_Level_RMS
    end

end