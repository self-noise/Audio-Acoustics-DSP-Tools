%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$
% Convolution operation - Carried out using FFTs
%
%           out = conv_FFT__MJN(y,ir)
%
% Inputs:
%         y - mono/stereo input sound data
%         ir - mono/stereo input IR data
%
% Outputs:
%         out - convolved mono/stereo sound data
%
% Length of output is exactly the same as using Matlab's "conv" function
%   --> i.e. nY + nIR - 1
%   --> As per Understanding DSP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$
function out = conv_FFT__MJN(y,ir)

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
