%---------------------------------------------------------------------------------------------%
% FUNCTION NAME AND SPECIFICATION
%
%   [~] = fftplot__MJN(x,Fs,plotStyle,varargin)
%---------------------------------------------------------------------------------------------%
% PURPOSE OF THIS FUNCTION:
%   This function takes a time input signal and calculates the DFT using Matlab's FFT function. 
%   It also calculates a frequency vector (f) and plots f versus FFT(x) over over the range 
%   0:Fs/2.
%   The program works for any number of input channels, as long as they are stored in a single 
%   vector x(:,channelIndex). Normally for a stereo signal channelIndex runs from 1:2. But if 
%   this program is being called to compare multiple different signals it may run from 1:N 
%   where N is the number of different signals. 
%   Each signal is plotted in a new subplot(N,1,channelIndex).
%---------------------------------------------------------------------------------------------%
% INPUTS:
%           x   :   Input signal to be analysed (mono or stereo or
%                   multichannel)
%           Fs  :   Sampling rate (Hz)
%           plot_style : 'stem' or 'plot'
%                       --> Choose to use a line plot or a signal
%                             processing style stem plot
%           varargin(1) : Cell array of strings for plot titles
%
% OUTPUTS:
%           A new figure is ploted showing the DFT of the input signal
%---------------------------------------------------------------------------------------------%
% GENERAL USAGE NOTES:
%   NOTE 1: 
%---------------------------------------------------------------------------------------------%
% CHANGES TO ADD AT SOME POINT IN THE FUTURE:
%   TBC
%---------------------------------------------------------------------------------------------%
% CHANGELOG:
%   2024-07-15:     Updated preamble for GitHub publication
%
%   2012-03-03:     Created this function as a quick way to do reasonably nice spectrum plots
%---------------------------------------------------------------------------------------------%
% Author:           Dr Mike Newton
% Date:             July 2024
% Location (local): [Matlab_root]/LIBRARY/MJN_Code_Library/MJN_DSP/
% GitHub location:  https://github.com/self-noise/Audio-Acoustics-DSP-Tools
%---------------------------------------------------------------------------------------------%
function MJN_fftplot(x,Fs,plot_style,varargin)

%------------------------------------------------------------------------%
% Error checking
%------------------------------------------------------------------------%

% Make sure there are enough titles if specified
switch isempty(varargin)
    case 0
        if length(varargin{1})~=size(x,2)
            error('Error: You have not specified enough plot titles as *varargin{1}* arguments. Use more title entries or remove *varargin{1}* from the function call.')
        end
end

%------------------------------------------------------------------------%
% Calculate frequency vector
%------------------------------------------------------------------------%
N   = size(x,1);
df  = Fs/N;
f   = 0:df:Fs-df;

%------------------------------------------------------------------------%
% Calculate and plot FFTs
%------------------------------------------------------------------------%
figFont = 20;
y_fft       = zeros(size(x));
nChannels   = size(x,2);
fig_ffts    = figure;
ax          = zeros(nChannels,1);
for jChannels = 1:nChannels
    y_fft(:,jChannels) = fft(x(:,jChannels));
    ax(jChannels)  = subplot(nChannels,1,jChannels);
    switch plot_style
        case 'stem'
            stem(f,abs(y_fft(:,jChannels)),'k','linewidth',2)
        case 'plot'
            plot(f,abs(y_fft(:,jChannels)),'k','linewidth',2)
        otherwise
            error('Plot style not specified correctly. Please use -stem- or -plot-')
    end
    xlim([0 Fs/2])
    %xlim([0 4000])
    xlabel('Frequency (Hz)','FontSize',figFont)
    ylabel('Amplitude','FontSize',figFont)
    switch isempty(varargin)
        case 0
            title(varargin{1}(jChannels),'FontSize',figFont)
        otherwise
            title(['DFT of input signal channel ' num2str(jChannels)],'FontSize',figFont)
    end
    grid on
end
linkaxes(ax,'x');

% Maximise figure
set(fig_ffts,'Position',get(0,'Screensize'));
set(ax,'FontSize',figFont-2)

