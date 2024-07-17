%---------------------------------------------------------------------------------------------%
% FUNCTION NAME AND SPECIFICATION
%
%   [params] = MJN_micCal__norm_to_Pa(params)
%---------------------------------------------------------------------------------------------%
% PURPOSE OF THIS FUNCTION:
%   Microphone calibration for acoustics research
%   Calculate a conversion factor(s) to transform digital numbers (scaled betweel -1:+1 by ADC)
%   into Pascals
%       ---->   i.e., take digital signal(s) and multiply by this number(s) to obtain signals
%               expressed in real world acoustic units
%---------------------------------------------------------------------------------------------%
% INPUTS:
%   params.
%       micCalType: 'singleRef' (only 1 mic calibration is provided, any others done using
%                           manufacturer specs relative to this)
%                   'individual' (each mic has its own calibration data as measured on the day
%                               -- BEST OPTION)
%                   Note that when only a single mic is specified, 'singleRef' and 'individual' are
%                   the same (either flag can be used)
%       N_Mics:     The number of microphones to be calibrated (this should match the number of
%                   unique 'params.micXX' entries as per item below
%       micXX:      At least one entry 'params.micXX' (starting with 'params.mic01') that provides
%                   metatdata for a microphone to be calibrated
%
% OUTPUTS:
%   params.
%       micXX.mic_calibration_factor :  This is a single number and provides scaling for the given
%                                       mic from digital/normalised units to real world acoustic
%                                       pressure (in Pa)
%       mic_cal_vector               :  A compiled vector of all the individual mic calibration
%                                       factors in a single place
%
%       Other entries within 'params' are not touched by this function
%---------------------------------------------------------------------------------------------%
% GENERAL USAGE NOTES:
%   NOTE 1: The input structure 'params' usually contains a bunch of stuff, but this function
%           only looks for and deals with the structure arguments noted below
%   NOTE 2: 'params.mic01' MUST be specified or calibration will not proceed
%   NOTE 3: You MUST manually specify the number of mics to be included in the calibration
%           ('params.N_Mics'), as a whole positive number. If you specify '3' mics
%           (params.N_Mics=3) but include fewer or more sets of mic metadata, an error will be
%           thrown
%---------------------------------------------------------------------------------------------%
% CHANGES TO ADD AT SOME POINT IN THE FUTURE:
%   "specData" changes size on each loop, adding a new entry to the structure for each channel in
%   sequence. This should be optimised/pre-allocated.
%---------------------------------------------------------------------------------------------%
% CHANGELOG:
%   2024-07-17:     Updated preamble for GitHub
%
%   2024-04-22:     Created this function as a standalone place to calibrate measurement
%                   microphones used in acoustics research
%---------------------------------------------------------------------------------------------%
% Author:           Dr Mike Newton
% Date:             July 2024
% Location (local): [Matlab_root]/LIBRARY/MJN_Code_Library/MJN_DSP/
% GitHub location:  https://github.com/self-noise/Audio-Acoustics-DSP-Tools
%---------------------------------------------------------------------------------------------%
function [params] = MJN_micCal__norm_to_Pa(params)


%----------------------------------------------------------------------------%
% Basic error checking on the number of mics specified
%----------------------------------------------------------------------------%
if ~isnumeric(params.N_Mics)
    % If the 'N_Mics' field isn't a number, throw an error
    error('Error: You have not specified a numeric field for the "params.N_Mics" field. Edit your metadata to make it a number.')
end
if params.N_Mics<=0
    error('Error: You have specified a nonsensical number of mics via the "params.N_Mics" field. Edit your metadata to make this an whole number of value 1 or larger.')
end
if floor(params.N_Mics)~=params.N_Mics
    error('Error: You have specified a nonsensical number of mics via the "params.N_Mics" field. Edit your metadata to make this an whole number of value 1 or larger.')
end


%----------------------------------------------------------------------------%
% SIGNAL CONVERSION FROM DIGITAL INTO ACOUSTIC UNITS
%----------------------------------------------------------------------------%

% Whatever workflow is employed via the main loop below, we *know* that at least 1 mic is being
% calibrated, and its specs and reference level are provided - therefore, do it right away
% Calibration for mic01 is straightforward, since we used a calibrator with known level
params.mic01.mic_calibration_factor   = params.mic01.ac_pressure_ref/params.mic01.sigLvl_dgtl_RMS;

% Initialise calibration vector (usually used to turn multichannel audio into scaled acoustic
% pressure signals using a single vector multiplication, applied elsewhere)
params.mic_cal_vector = [params.mic01.mic_calibration_factor];

if params.N_Mics>1
    % In cases where more than a single mic is to be calibrated, we need to select the correct
    % workflow to follow
    switch params.micCalType
        case 'individual'
            % In this workflow, we just repeat the process for mic01 with the remaining mics,
            % since we have exactly the data required to do so
            for n=2:params.N_Mics
                params.(['mic' num2str(n,'%02.f')]).mic_calibration_factor ...
                    = params.(['mic' num2str(n,'%02.f')]).ac_pressure_ref/params.(['mic' num2str(n,'%02.f')]).sigLvl_dgtl_RMS;

                params.mic_cal_vector = [params.mic_cal_vector,params.(['mic' num2str(n,'%02.f')]).mic_calibration_factor];
            end            

        case 'singleRef'
            % 'mic01' provides a single master calibration reference, and calibration of
            % remaining microphones is achieved by applying the manufacturer's sensitivity
            % specifications (which must be provided by metadata for each mic), plus the
            % recorded preamp gain used on the day the audio was recorded
            for n=2:params.N_Mics
                delta_mVPa      = params.(['mic' num2str(n,'%02.f')]).mVPA/params.mic01.mVPA;
                delta_preamp    = 10^((params.(['mic' num2str(n,'%02.f')]).preampGain-params.mic01.preampGain)/20);

                % Calculate the calibration factor for the other microphones based
                % on that identified for the reference (taking into account
                % sensitivities and applied preamp gains)
                params.(['mic' num2str(n,'%02.f')]).mic_calibration_factor = params.mic01.mic_calibration_factor*delta_mVPa*delta_preamp;

                % Fill out the calibration vector
                params.mic_cal_vector = [params.mic_cal_vector,params.(['mic' num2str(n,'%02.f')]).mic_calibration_factor];

            end
    end
end
