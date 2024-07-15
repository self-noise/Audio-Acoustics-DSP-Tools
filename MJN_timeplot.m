%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Program          : timeplot__MJN.m
% Program category : AMT MSc DSP Helper program
% Author           : Michael J. Newton
% Date             : 03/2012
%
% Overview: This function takes any number of time domain input signals and
%           simply plots them, either against sample number or against time
%           in seconds.
%           The program works for any number of input channels, as long as
%           they are stored in a single vector x(:,channelIndex). Normally
%           for a stereo signal channelIndex runs from 1:2. But if this
%           program is being called to compare multiple different signals
%           it may run from 1:N where N is the number of different signals.
%           Each signal is plotted in a new subplot(N,1,channelIndex).
%
% Inputs
%           x       :   Input signal to be plotted (mono or stereo or
%                       multichannel)
%           Fs      :   Sampling rate (Hz)
%           units   :   'seconds' or 'samples'
%           normOp  :   0 for no signal normalisation, 1 for normalisation
%           varargin:   y limits for each plot
%                       --> varargin(1) is a vector [ylow yhigh] for
%                           signal x(:,1)
%                       --> If setting any varargin entries, make sure to
%                       pad with [] for all other signals not controlled,
%                       e.g. for 3 signals, only first 2 to be controlled,
%                       use:
%                           [ylow1 yhigh1],[ylow1 yhigh1],[]
%
% Outputs
%           A new figure is ploted showing the time domain plots of the
%           input signals
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function timeplot__MJN(x,Fs,units,normOp,varargin)


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Length of input signals
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
N = size(x,1);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Perform normalisation if required
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
switch normOp
    case 0
    case 1
        % Find maximum abs value of each column, then tile vertically and
        % perform element by element division to normalise
        w = repmat(max(abs(x)),[N 1]);
        x = x./w;
    otherwise
        error('Error: Supplied normOp is not 0 (no normalisation) or 1 (normalisation)')
end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plot the time signals, customising the plot limits if required
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figFont = 20;
nChannels   = size(x,2);
fig_time    = figure;
ax          = zeros(nChannels,1);
for jChannels = 1:nChannels
    ax(jChannels)  = subplot(nChannels,1,jChannels);
    switch units
        case 'samples'
            plot(x(:,jChannels),'k','linewidth',2)
            xlabel('Time (samples)','FontSize',figFont)
        case 'seconds'
            plot(0:1/Fs:(N/Fs)-(1/Fs),x(:,jChannels),'k','linewidth',2)
            xlabel('Time (s)','FontSize',figFont)
    end
    ylabel('Amplitude','FontSize',figFont)
    title(['Time domain plot of input signal channel ' num2str(jChannels)],'FontSize',figFont)
    grid on
    set(gca,'FontSize',figFont-2)
    
    % Apply ylimits for current signal, if requested        
    if ~isempty(varargin)
        ylim(varargin{jChannels})        
    end
end
linkaxes(ax,'x');

% Maximise figure
set(fig_time,'Position',get(0,'Screensize'));



