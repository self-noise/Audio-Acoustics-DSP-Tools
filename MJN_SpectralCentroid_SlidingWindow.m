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
%               .params.f0_nominal      The reference frequency for NSC
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
%           NSC_frames                  Normalised spectral centroid (by f0)
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
%   2024-09-30: Lot of updates and changes. 
%               Added optional noise thresholding for SC computation
%               Normalised the 'Hann' window (if used) so its average value is 1 (when 
%               computing the magnitude spectrum)
%               Commenced adding optional power spectrum or magnitude spectrum option (not yet
%               complete as of 02 Oct 2024)
%
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
function [SC_normalised_frames,SC_frames,frame_RMS_levels,tVec_FrameCentres] = MJN_SpectralCentroid_SlidingWindow(params,audioData,n_File,Path_File,options)

N_channels                  = size(audioData,2);                                        % Number of mic channels
Ts__thisFile                = 1/params.Fs;                                              % Sample period
N_tot__thisFile             = length(audioData);                                        % Number of samples in current file
T_end__thisFile             = N_tot__thisFile/params.Fs;                                % Duration (in sec) of current file
tVec__thisFile              = [0:Ts__thisFile:T_end__thisFile-Ts__thisFile];            % Time vector for current file (including all samples)


%---------------------------------------------------------------------------------------------%
% If requested, load up a figure that is ready to receive the raw spectral data
%---------------------------------------------------------------------------------------------%
switch options.PlotSpec
    case 1
        fig(n_File) = figure;
        for nPlotChannel = 1:N_channels
            ax(nPlotChannel)  = subplot(N_channels,1,nPlotChannel);
            hold on
            grid on
            xlabel('Frequency (Hz)')

            switch options.PlotSpec_Scale
                case 'linear'
                    ylabel('DFT magnitude (Pa)')
                case 'log'
                    ylabel('DFT magnitude (dB re 1 Pa)')
            end

            title(['File: ' Path_File ' --- ' params.(['mic0' num2str(nPlotChannel)]).location],'Interpreter','none')
            xlim([params.f0_rangeLow,params.SC_upper_bound_Hz])
        end
end


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

    % Zero pad the input data so it perfectly fits an integer number of frames NF
    audioData = [audioData;zeros(N_tot__thisFile_pad-N_tot__thisFile,N_channels)];

    % Make an empty vector to record the fundamental frequency (f0) estimate for each frame
    params.f0_estimate          = zeros(NF,1);

    % Make an empty vector to record the lower frequency bound of the SC analysis window (which
    % depends on the estimated f0). Note that the upper bound is fixed as a parameter at the 
    % start of analysis.
    params.SC_lower_bound_Hz    = zeros(NF,1);

    % Define a windowing function to apply to each frame
    switch params.windowType
        case 'Hann'
            % Manual computation of Hann window
            win                  = 0.5*(1 - cos(2*pi*(0:N_frame-1)'/N_frame));
            % Normalise the window so its average value is 1, to avoid undesired amplitude
            % scaling of the frame to which it's applied
            win                  = win/(mean(win));
        case 'Rect'
            % All ones (no amplitude scaling needed) means a rectangular window
            win                  = ones(N_frame,1);
        case 'Blackman'            
            % Blackman window (manual computation
            win                  = (0.42 - 0.5*cos(2*pi*(0:N_frame-1)/(N_frame-1)) + 0.08*cos(4*pi*(0:N_frame-1)/(N_frame-1)))';
            win                  = win/(mean(win));            
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

    % % Frequency (Hz) to bin indices: Convert the user's lower/upper bounds
    % % in Hz into a bin index (note that rounding errors will occur)
    % [~,SC_lower_bin] = min(abs(fVec - params.SC_lower_bound_Hz));
    [~,SC_upper_bin] = min(abs(fVec - params.SC_upper_bound_Hz));

    % Big empty matrices to store all the frames (the 'raw' version is used
    % to compute the RMS level only, for use in estimating the signal
    % level during each frame)
    audioData_slices_raw                = zeros(N_frame,NF);
    audioData_slices                    = zeros(N_frame,NF);

    % Loop over the number of channels (i.e., microphones used)
    for nChan = 1:N_channels

        %-------------------------------------------------------------------------------------%
        % FRAME-BY-FRAME READ IN:
        % Read the frames into the big matrix of frame slices
        %-------------------------------------------------------------------------------------%        
        for m=0:NF-1
            indRange                    = (m*HA)+1:(m*HA)+N_frame;
            audioData_slices_raw(:,m+1) = audioData(indRange,nChan);
            audioData_slices(:,m+1)     = audioData(indRange,nChan).*win;

            % If requested, estimate f0 for the current frame (but only for channel 1)
            switch params.f0_Method
                case 'estimate_per_frame'
                    if nChan==1
                        this_f0                 = pitch(audioData_slices(:,m+1),params.Fs, Method='SRH', WindowLength=N_frame, OverlapLength=0, Range=[params.f0_rangeLow,params.f0_rangeHigh], MedianFilterLength=1);
                        params.f0_estimate(m+1) = this_f0;

                        % Set the lower bound for the SC analysis, based on the expected f0 of the
                        % given note. Usual expectation is that this is determined by the expected
                        % playing frequency minus a fixed offset to ensure that the DFT selection
                        % includes everything in that vicinity whole. Usual: '-100' is about right
                        params.SC_lower_bound_Hz(m+1) = params.f0_estimate(m+1) - params.f_lower_offset;
                    end
                case 'nominal'
                    params.f0_estimate(m+1) = params.f0_nominal;
                    params.SC_lower_bound_Hz(m+1) = params.f0_estimate(m+1) - params.f_lower_offset;
                otherwise
                    error('Error in "params.f0_Method": f0 can only be determined per frame ("estimate_per_frame"), or via a nominal setting ("nominal"). Set "params.f0_Method" to one of these.')
            end

            % Frequency (Hz) to bin indices: Convert the user's lower/upper bounds
            % in Hz into a bin index (note that rounding errors will occur)
            [~,SC_lower_bin] = min(abs(fVec - params.SC_lower_bound_Hz(m+1)));            
        end        

        % Compute the RMS level of each slice (RAW slice, NOT including
        % windowing)
        % Note: If this function has been given pre-calibrated acoustic
        % signals (i.e., units of Pa) then this will be the RMS sound
        % pressure in units of Pa for the given slice
        audioData_slices_Level_RMS      = rms(audioData_slices_raw);

        % Take the DFT of all the slices (using the normed version of the DFT as per:
        % https://selfnoise.co.uk/resources/signals-and-dft/dft-scaling/)
        audioData_slices_DFT            = (2/N_frame)*fft(audioData_slices);
        audioData_slices_DFT_mag        = abs(audioData_slices_DFT);

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
        fVec__halfSpec_bounds                       = fVec__halfSpec;
        fVec__halfSpec_bounds(1:SC_lower_bin)       = 0;
        fVec__halfSpec_bounds(SC_upper_bin+1:end)   = 0;


        %-------------------------------------------------------------------------------------%
        % % Testing/checking only: Save time and DFT slices
        %-------------------------------------------------------------------------------------%
        % timeSliceRAW_saved{nChan}.data  = audioData_slices_raw;
        % timeSlice_saved{nChan}.data     = audioData_slices;
        % DFTSlice_saved{nChan}.data      = audioData_slices_DFT_mag_halfSpec;
        % save('blahblah','timeSliceRAW_saved','timeSlice_saved','DFTSlice_saved')
        %-------------------------------------------------------------------------------------%

        %-------------------------------------------------------------------------------------%
        % NOISE THRESHOLDING: If requested, apply a noise floor threshold, setting all spectral 
        % magnitudes below this to zero
        %-------------------------------------------------------------------------------------%
        switch params.optionNoiseThresh
            case 0
                % Do nothing
            case 1
                % Apply thresholding
                audioData_slices_DFT_mag_halfSpec_copy = audioData_slices_DFT_mag_halfSpec;
                theseBinsBelowNoiseThreshold = find(audioData_slices_DFT_mag_halfSpec < params.NoiseThresh_Pa(nChan));
                audioData_slices_DFT_mag_halfSpec(theseBinsBelowNoiseThreshold) = 0;
            otherwise
                error('Error: the parameter "params.optionNoiseThresh" can only be set to 0 or 1')
        end

        %-------------------------------------------------------------------------------------%
        % SPECTRAL CENTROID: Actually calcualte it! And the normalised version!
        % Mike's method (each frame SC via vector*matrix multiplication)
        %-------------------------------------------------------------------------------------%
        these_SC_frames                 = (fVec__halfSpec_bounds'*audioData_slices_DFT_mag_halfSpec)./sum(audioData_slices_DFT_mag_halfSpec,1);
        these_SC_normalised_frames      = these_SC_frames'./params.f0_estimate;
        %these_SC_normalised_frames      = these_SC_frames./params.f0_nominal;


        % Clean up the latter parts of the file (if needed), chosing the first
        % channel (i.e., usually Mic01) as the reference point to ensure that all
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

        % Find the quietest frame
        [~,indexSoundLevelMin]          = min(audioData_slices_Level_RMS);

        % Make a vector to represent harmonic multiples of the expected playing frequency (use
        % a resolution of 0.5 Hz, regardless of Fs)
        fStep2                          = 0.5;
        %fVec_Harmonics                  = (0:fStep2:params.Fs - fStep2)';
        fVec_Harmonics                  = (0:fStep2:params.Fs/2)';      

        % Nearest integer bin index corresponding to expected playing frequency relative to the
        % fVec_Harmonics axis
        k_f_expected                    = (round(mean(params.f0_estimate)*2)*0.5)*2;
        %k_f_expected                    = (round(params.f_expected*2)*0.5)*2;        

        AVec_harmonics_linear           = zeros(length(fVec_Harmonics),1);        
        AVec_harmonics_linear(k_f_expected+1:k_f_expected:end)  = 1;        

        AVec_harmonics_dB               = -100*ones(length(fVec_Harmonics),1);
        AVec_harmonics_dB(k_f_expected+1:k_f_expected:end)      = 0;

        % If requested, add raw spectral data to the temporary plot figure (but only include
        % slices up until the max SPL window, and only the frequency bins included in
        % calculating the NSC)
        figLineHarmonics = 1;
        switch options.PlotSpec
            case 0
            case 1
                indRangeAnalysisWindow = find(fVec__halfSpec_bounds);
                switch options.PlotSpec_Range
                    case 'all'
                        switch options.PlotSpec_Scale
                            case 'linear'                                
                                if params.optionNoiseThresh
                                    plot([mean(params.SC_lower_bound_Hz),params.SC_upper_bound_Hz],params.NoiseThresh_Pa(nChan)*ones(2,1),'r','LineWidth',2,'Parent', ax(nChan))
                                    plot(fVec__halfSpec(indRangeAnalysisWindow),audioData_slices_DFT_mag_halfSpec_copy(indRangeAnalysisWindow,1:indexSoundLevelMax),'g','Parent', ax(nChan))
                                end
                                plot(fVec__halfSpec(indRangeAnalysisWindow),audioData_slices_DFT_mag_halfSpec(indRangeAnalysisWindow,1:indexSoundLevelMax),'Parent', ax(nChan))
                                %plot(fVec_Harmonics,AVec_harmonics_linear,'Parent', ax(nChan),'LineWidth',figLineHarmonics,'Color','k');
                            case 'log'                                
                                if params.optionNoiseThresh
                                    plot([mean(params.SC_lower_bound_Hz),params.SC_upper_bound_Hz],20*log10(params.NoiseThresh_Pa(nChan))*ones(2,1),'r','LineWidth',2,'Parent', ax(nChan))
                                    plot(fVec__halfSpec(indRangeAnalysisWindow),20*log10(audioData_slices_DFT_mag_halfSpec_copy(indRangeAnalysisWindow,1:indexSoundLevelMax)),'g','Parent', ax(nChan))
                                end
                                plot(fVec__halfSpec(indRangeAnalysisWindow),20*log10(audioData_slices_DFT_mag_halfSpec(indRangeAnalysisWindow,1:indexSoundLevelMax)),'Parent', ax(nChan))
                                %plot(fVec_Harmonics,AVec_harmonics_dB,'Parent', ax(nChan),'LineWidth',figLineHarmonics,'Color','k');
                        end

                    case 'peak'
                        switch options.PlotSpec_Scale
                            case 'linear'                                
                                if params.optionNoiseThresh
                                    plot([mean(params.SC_lower_bound_Hz),params.SC_upper_bound_Hz],params.NoiseThresh_Pa(nChan)*ones(2,1),'r','LineWidth',2,'Parent', ax(nChan))
                                    plot(fVec__halfSpec(indRangeAnalysisWindow),audioData_slices_DFT_mag_halfSpec_copy(indRangeAnalysisWindow,indexSoundLevelMax),'g','Parent', ax(nChan))
                                end
                                plot(fVec__halfSpec(indRangeAnalysisWindow),audioData_slices_DFT_mag_halfSpec(indRangeAnalysisWindow,indexSoundLevelMax),'Parent', ax(nChan))
                                %plot(fVec_Harmonics,AVec_harmonics_linear,'Parent', ax(nChan),'LineWidth',figLineHarmonics,'Color','k');
                            case 'log'                                
                                if params.optionNoiseThresh
                                    plot([mean(params.SC_lower_bound_Hz),params.SC_upper_bound_Hz],20*log10(params.NoiseThresh_Pa(nChan))*ones(2,1),'r','LineWidth',2,'Parent', ax(nChan))
                                    plot(fVec__halfSpec(indRangeAnalysisWindow),20*log10(audioData_slices_DFT_mag_halfSpec_copy(indRangeAnalysisWindow,indexSoundLevelMax)),'g','Parent', ax(nChan))
                                end
                                plot(fVec__halfSpec(indRangeAnalysisWindow),20*log10(audioData_slices_DFT_mag_halfSpec(indRangeAnalysisWindow,indexSoundLevelMax)),'Parent', ax(nChan))
                                %plot(fVec_Harmonics,AVec_harmonics_dB,'Parent', ax(nChan),'LineWidth',figLineHarmonics,'Color','k');
                        end

                    case 'min'
                        switch options.PlotSpec_Scale
                            case 'linear'                                
                                if params.optionNoiseThresh
                                    plot([mean(arams.SC_lower_bound_Hz),params.SC_upper_bound_Hz],params.NoiseThresh_Pa(nChan)*ones(2,1),'r','LineWidth',2,'Parent', ax(nChan))
                                    plot(fVec__halfSpec(indRangeAnalysisWindow),audioData_slices_DFT_mag_halfSpec_copy(indRangeAnalysisWindow,indexSoundLevelMin),'g','Parent', ax(nChan))
                                end
                                plot(fVec__halfSpec(indRangeAnalysisWindow),audioData_slices_DFT_mag_halfSpec(indRangeAnalysisWindow,indexSoundLevelMin),'Parent', ax(nChan))
                                %plot(fVec_Harmonics,AVec_harmonics_linear,'Parent', ax(nChan),'LineWidth',figLineHarmonics,'Color','k');
                            case 'log'
                                if params.optionNoiseThresh
                                    plot([mean(params.SC_lower_bound_Hz),params.SC_upper_bound_Hz],20*log10(params.NoiseThresh_Pa(nChan))*ones(2,1),'r','LineWidth',2,'Parent', ax(nChan))
                                    plot(fVec__halfSpec(indRangeAnalysisWindow),20*log10(audioData_slices_DFT_mag_halfSpec_copy(indRangeAnalysisWindow,indexSoundLevelMin)),'g','Parent', ax(nChan))
                                end
                                plot(fVec__halfSpec(indRangeAnalysisWindow),20*log10(audioData_slices_DFT_mag_halfSpec(indRangeAnalysisWindow,indexSoundLevelMin)),'Parent', ax(nChan))
                                %plot(fVec_Harmonics,AVec_harmonics_dB,'Parent', ax(nChan),'LineWidth',figLineHarmonics,'Color','k');
                        end

                    otherwise
                        error('Error: You can only set "options.PlotSpec_Range" to "all", "peak", or "min".')
                end
                %ax.YLim = [10,-140];
                set(ax,'ylim',[-140,10])
            otherwise
                error('Error: You can only set "options.PlotSpec" to "0" or "1".')
        end

        clear audioData_slices_raw audioData_slices audioData_slices_Level_RMS
    end

end