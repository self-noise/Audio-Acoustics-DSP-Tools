%---------------------------------------------------------------------------------------------%
% FUNCTION NAME AND SPECIFICATION
%
%   out = MJN_deconv_fft(in,out)
%---------------------------------------------------------------------------------------------%
% PURPOSE OF THIS FUNCTION:
%   Performs linear (i.e., 'acyclic') deconvolution between a pair of vectors in the frequency 
%   domain, via the FFT algorithm (with appropriate zero-padding). In typical use this would be 
%   deconvolution of the 'ouput' of a system from an 'input' to it, in order to recover the 
%   system's impulse response
%
%   The function accepts either single-channel or two-channel files (for both 'dry' and 'ir')
%---------------------------------------------------------------------------------------------%
% INPUTS:
%    in          : Mono/stereo input data
%    out         : Mono/stereo output data
%
% OUTPUTS:
%   ir           : Deconvolved mono/stereo sound data (usually an impulse response)
%---------------------------------------------------------------------------------------------%
% GENERAL USAGE NOTES:
%   Length of the output is exactly the same as when using Matlab's inbuilt "conv" function
%   --> i.e. N_Y (number of input samples) + N_IR (number of ir samples) - 1
%   --> More info, e.g., https://ccrma.stanford.edu/~jos/sasp/Acyclic_Convolution_Matlab.html
%---------------------------------------------------------------------------------------------%
% CHANGES TO ADD AT SOME POINT IN THE FUTURE:
%   TBA
%---------------------------------------------------------------------------------------------%
% CHANGELOG:
%   2024-07-17:     Changed function name and updated preamble comments
%
%   2018-04-02:     Created this function as a quick way to do linear deconvolution
%---------------------------------------------------------------------------------------------%
% Author:           Dr Mike Newton
% Date:             July 2024
% Location (local): [Matlab_root]/LIBRARY/MJN_Code_Library/MJN_DSP/
% GitHub location:  https://github.com/self-noise/Audio-Acoustics-DSP-Tools
%---------------------------------------------------------------------------------------------%
function ir = MJN_deconv_fft(in,out)

nIN         = max(size(in));
nIN_Chans   = min(size(in));

if nIN_Chans > 2
    error('Error: Your input sound file has more than 2 channels.')
end
if nIN_Chans>=nIN
    error('Error: Your stereo input sound data has rows and columns ordered incorrectly')
end

nOUT       = max(size(out));
nOUT_Chans = min(size(out));

if nOUT_Chans > 2
    error('Error: Your output sound file has more than 2 channels.')
end
if nOUT_Chans>=nOUT
    error('Error: Your stereo output sound file has rows and columns ordered incorrectly')
end

% if nIN_Chans==nOUT_Chans && nIN_Chans == 2
%     disp('Note: Your input sound file is stereo, and your impulse response is stereo; full stereo convolution carrried out.')
% elseif nIN_Chans==nOUT_Chans && nIN_Chans == 1
%     disp('Note: Your input sound file is mono, and your impulse response is mono; mono convolution carrried out.')
% end

% if nIN_Chans<nOUT_Chans
%     disp('Note: Your input sound file is stereo, but your output sound file is mono; proceeding to deconvolve with Left channel only.')
% elseif nChans>nIRChans
%     disp('Note: Your input sound file is stereo, but your output sound file is mono; convolving both channels with the mono IR.')
% end


switch nIN_Chans
    case 1
        in      = [in;zeros(nOUT-1,1)];        
        out     = [out(:,1);zeros(nIN-1,1)];
        ir      = ifft(fft(out)./fft(in));
    case 2
        error('More than 1 input channels supplied - cannot proceed.')
%         switch nIRChans
%             case 1
%                 inL             = [in(:,1);zeros(nIR-1,1)];
%                 inR             = [in(:,2);zeros(nIR-1,1)];
%                 out            = zeros(nY+nIR-1,2);
%                 out             = [ir(:,1);zeros(nY-1,1)];
%                 ir(:,1)       = ifft(fft(yL).*fft(ir));
%                 ir(:,2)       = ifft(fft(yR).*fft(ir));
%             case 2
%                 yL             = [y(:,1);zeros(nIR-1,1)];
%                 yR             = [y(:,2);zeros(nIR-1,1)];
%                 out            = zeros(nY+nIR-1,2);
%                 irL            = [ir(:,1);zeros(nY-1,1)];
%                 irR            = [ir(:,2);zeros(nY-1,1)];
%                 out(:,1)       = ifft(fft(yL).*fft(irL));
%                 out(:,2)       = ifft(fft(yR).*fft(irR));
%         end
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
