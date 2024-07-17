%---------------------------------------------------------------------------------------------%
% FUNCTION NAME AND SPECIFICATION
%
%   out = MJN_conv_fft(y,ir)
%---------------------------------------------------------------------------------------------%
% PURPOSE OF THIS FUNCTION:
%   Performs linear (i.e., 'acyclic') convolution between a pair of vectors in the frequency 
%   domain, via the FFT algorithm (with appropriate zero-padding). In typical use this would be 
%   convolution of a 'dry' audio signal with an impulse response of some kind 'ir'.
%
%   The function accepts either single-channel or two-channel files (for both 'dry' and 'ir')
%---------------------------------------------------------------------------------------------%
% INPUTS:
%    y          : Mono/stereo input sound data
%    ir         : Mono/stereo input IR data
%
% OUTPUTS:
%   out         : Convolved mono/stereo sound data
%---------------------------------------------------------------------------------------------%
% GENERAL USAGE NOTES:
%   Length of the output is exactly the same as when using Matlab's inbuilt "conv" function
%   --> i.e. N_Y (number of input samples) + N_IR (number of ir samples) - 1
%   --> More info, e.g., https://ccrma.stanford.edu/~jos/sasp/Acyclic_Convolution_Matlab.html
%---------------------------------------------------------------------------------------------%
% CHANGES TO ADD AT SOME POINT IN THE FUTURE:
%   TBC
%---------------------------------------------------------------------------------------------%
% CHANGELOG:
%   2024-07-17:     Changed function name and updated preamble comments
%   2018-04-02:     Created this function as a quick way to do linear convolution
%---------------------------------------------------------------------------------------------%
% Author:           Dr Mike Newton
% Date:             July 2024
% Location (local): [Matlab_root]/LIBRARY/MJN_Code_Library/MJN_DSP/
% GitHub location:  https://github.com/self-noise/Audio-Acoustics-DSP-Tools
%---------------------------------------------------------------------------------------------%
function out = MJN_conv_fft(y,ir)

nY       = max(size(y));
nChans   = min(size(y));

if nChans > 2
    error('Error: Your input sound file has more than 2 channels.')
end
if nChans>=nY
    error('Error: Your stereo input sound data has rows and columns ordered incorrectly')
end

nIR      = max(size(ir));
nIRChans = min(size(ir));

if nIRChans > 2
    error('Error: Your impulse response file has more than 2 channels.')
end
if nIRChans>=nIR
    error('Error: Your stereo impulse response data has rows and columns ordered incorrectly')
end

if nChans==nIRChans && nChans == 2
    disp('Note: Your input sound file is stereo, and your impulse response is stereo; full stereo convolution carrried out.')
elseif nChans==nIRChans && nChans == 1
    disp('Note: Your input sound file is mono, and your impulse response is mono; mono convolution carrried out.')
end

if nChans<nIRChans
    %disp('Note: Your impulse response file is stereo, but your input is mono; proceeding with Left channel only.')
    disp('Note: Your impulse response file is stereo, but your input is mono; convolving mono input with each IR channel to produce stereo output file.')
elseif nChans>nIRChans
    disp('Note: Your input sound file is stereo, but your impulse response is mono; convolving both channels with the mono IR.')
end


switch nChans
    case 1
        y   = [y;zeros(nIR-1,1)];
        switch nIRChans
            case 1
                ir  = [ir(:,1);zeros(nY-1,1)];
                out = ifft(fft(y).*fft(ir));
            case 2
                irL            = [ir(:,1);zeros(nY-1,1)];
                irR            = [ir(:,2);zeros(nY-1,1)];
                out(:,1)       = ifft(fft(y).*fft(irL));
                out(:,2)       = ifft(fft(y).*fft(irR));
        end
        
    case 2
        switch nIRChans
            case 1
                yL             = [y(:,1);zeros(nIR-1,1)];
                yR             = [y(:,2);zeros(nIR-1,1)];
                out            = zeros(nY+nIR-1,2);
                ir             = [ir(:,1);zeros(nY-1,1)];
                out(:,1)       = ifft(fft(yL).*fft(ir));
                out(:,2)       = ifft(fft(yR).*fft(ir));
            case 2
                yL             = [y(:,1);zeros(nIR-1,1)];
                yR             = [y(:,2);zeros(nIR-1,1)];
                out            = zeros(nY+nIR-1,2);
                irL            = [ir(:,1);zeros(nY-1,1)];
                irR            = [ir(:,2);zeros(nY-1,1)];
                out(:,1)       = ifft(fft(yL).*fft(irL));
                out(:,2)       = ifft(fft(yR).*fft(irR));
        end
end


% ir       = [ir;zeros(nY-1,1)];
% out = zeros(nY+nIR-1,nChans);
%
% for jChans = 1:nChans
%     for jIRChans = 1:nIRChans
%         thisY         = [y(:,jChans);zeros(nIR-1,1)];
%         out(:,jChans) = ifft(fft(thisY).*fft(ir));
%     end
% end
