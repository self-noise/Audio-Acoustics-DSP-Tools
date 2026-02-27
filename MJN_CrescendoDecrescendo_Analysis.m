
function MJN_CrescendoDecrescendo_Analysis(options,params,paths)


%----------------------------------------------------------------------------%
% MICROPHONE CALIBRATION
%----------------------------------------------------------------------------%
params = MJN_micCal__norm_to_Pa(params);


%----------------------------------------------------------------------------%
% FILE READING AND ANALYSIS
%----------------------------------------------------------------------------%

% How many files to analyse?
N_Files = numel(paths.Path_File);

% Empty structure to store all the analysed data field (SC, sig level, etc)
data_Analysed = struct('NSC_frames', {},'SC_frames', {}, 'sigLevel_Pa', {}, 'tVec_FrameCentres', {});

% A record of the largest single NSC and SPL values in this particular run
% of the script, to enable quick plotting of the data if required with
% appropriate axes limits
recordMaxNSC    = 0;
counterMaxNSC   = 0;
recordMaxSPL    = 0;
counterMaxSPL   = 0;

% Main loading and analsing loop
for n_File = 1:N_Files

    % Read in the current WAV file
    paths.Path_Full = [paths.Path_Folder_Main '/' paths.Path_Folder_Main_Data '/' paths.Path_Folder_thisExp '/' paths.Path_Folder_thisNote '/' ...
        paths.Path_Folder_thisConfig '/' paths.Path_Folder_thisPlayer '/' paths.Path_File{n_File}];
    [thisSignal,params.Fs]   = audioread(paths.Path_Full);

    % Apply the microphone calibration vector channel-by-channel to produce
    % a calibrated set of signals in units of Pa
    thisSignal_Pa = thisSignal.*params.mic_cal_vector;

    % If the recording is a decrescendo, then time-flip the data before sending for analysis (a
    % record of whether or not this is done is kept in 'params')
    switch paths.Path_Folder_thisConfig
        case 'Decrescendo'
            thisSignal_Pa = flipud(thisSignal_Pa);
    end

    %-----------------------------------------------------------------%
    % SC analysis function call (file located elsewhere on Matlab path)
    %-----------------------------------------------------------------%
    [data_Analysed(n_File).NSC_frames,data_Analysed(n_File).SC_frames,data_Analysed(n_File).sigLevel_Pa,data_Analysed(n_File).tVec_FrameCentres,...
        params.f0_estimate(n_File,1).value,params.SC_lower_bound_Hz(n_File,1).value] = ...
        MJN_SpectralCentroid_SlidingWindow(params,thisSignal_Pa,options,n_File,paths.Path_File{n_File});
    %-----------------------------------------------------------------%

    % This keeps track of the max value of NSC to aid in plotting
    switch counterMaxNSC
        case 1
            recordMaxNSC = max(data_Analysed(counterMaxNSC).NSC_frames(:));
        otherwise
            if  max(data_Analysed(n_File).NSC_frames(:))>recordMaxNSC
                recordMaxNSC = max(data_Analysed(n_File).NSC_frames(:));
            end
    end

    % This keeps track of the max value of SPL to aid in plotting
    switch counterMaxSPL
        case 1
            recordMaxSPL = max(data_Analysed(counterMaxSPL).NSC_frames(:));
        otherwise
            if  max(data_Analysed(n_File).sigLevel_Pa(:))>recordMaxSPL
                recordMaxSPL = max(data_Analysed(n_File).sigLevel_Pa(:));
            end
    end
end

% Save .MAT file if required
params.metadataExperiment.Experiment    = paths.Path_Folder_thisExp;
params.metadataExperiment.Note          = paths.Path_Folder_thisNote;
params.metadataExperiment.Config        = paths.Path_Folder_thisConfig;
params.metadataExperiment.Player        = paths.Path_Folder_thisPlayer;
params.metadataExperiment.Takes         = paths.Path_File;

% Prepare output filenames, whether or not they are used
filenameStem_SavedData = [...
    paths.Path_Folder_Main '/' ...
    paths.Path_Folder_Main_Data '/' ...
    paths.Path_Folder_thisExp '/' ...
    paths.Path_Folder_thisNote '/' ...
    paths.Path_Folder_thisConfig '/' ...
    paths.Path_Folder_thisPlayer];
filenameFile_SavedData = [...
    paths.Path_Folder_thisExp '_' ...
    paths.Path_Folder_thisNote '_' ...
    paths.Path_Folder_thisConfig '_' ...
    paths.Path_Folder_thisPlayer '__' ...
    [params.spectrumType] '_' ...
    ['win' num2str(params.windowSize_ms) 'ms'] '_' ...
    ['ovlp' num2str(round(params.windowSize_ms*params.overlapFraction)) 'ms'] '_' ...
    ['SC' num2str(params.SC_upper_bound_Hz) 'Hz'] '__' ...
    paths.Path_File{1,1}(1:end-4) '-' paths.Path_File{N_Files,1}(1:end-4)...
    ];
params.metadataExperiment.Filename_Output = filenameFile_SavedData;

% If requested write a Matlab .MAT data file with output data
if options.SaveFile==1
    %params.metadataExperiment.N_Mics        = N_Mics;
    save([filenameStem_SavedData '/' filenameFile_SavedData '.mat'],'data_Analysed','params')
    disp(['Data has been saved to: ' filenameFile_SavedData '.mat'])
end


% Do a quick plot of everything if required
switch options.Plot
    case 1

        %SP_RMS_LowerRangePlot   = 0;
        %SP_RMS_UpperRangePlot   = ceil(recordMaxSPL*4)/4;
        % SP_RMS_LowerRangePlot   = 0;
        % SP_RMS_UpperRangePlot   = 4;

        %NSC_LowerRangePlot      = 1;        
        %NSC_UpperRangePlot      = ceil(recordMaxNSC*5)/5;
        %NSC_UpperRangePlot      = 30;

        % NSC_LowerRangePlot      = options.PlotRangeNSC(1);
        % NSC_UpperRangePlot      = options.PlotRangeNSC(2);

        figFont = 16;
        figLine = 2;

        fig1 = figure;

        for nMic = 1:params.N_Mics
            ax(nMic) = subplot(params.N_Mics,1,nMic);
            for nPlot = 1:N_Files
                hold on
                grid on
                plotHandle1(nPlot) = plot(data_Analysed(nPlot).sigLevel_Pa(:,nMic),data_Analysed(nPlot).NSC_frames(:,nMic),'Linewidth',figLine,'MarkerSize',5);

                %xlim([SP_RMS_LowerRangePlot,SP_RMS_UpperRangePlot])
                %ylim([NSC_LowerRangePlot,NSC_UpperRangePlot])
                %ylim([NSC_LowerRangePlot,8])
                %xlim([SP_RMS_LowerRangePlot,3])
                % xlim([SP_RMS_LowerRangePlot,SP_RMS_UpperRangePlot])
                %ylim([NSC_LowerRangePlot,NSC_UpperRangePlot])
                %ylim([NSC_LowerRangePlot,5])
                xlim([options.PlotRangeLev])
                ylim([options.PlotRangeNSC])

                ylabel(['Normalised spectral centroid'],'FontSize',figFont,'Interpreter','Latex')
                set(gca,'FontSize',figFont,'TickLabelInterpreter','latex')
            end

            titleString = {[params.metadataExperiment.Experiment ' (' paths.Path_Folder_thisNote ', ' paths.Path_Folder_thisPlayer ') - ' paths.Path_Folder_thisConfig]; [params.(['mic' num2str(nMic,'%02.f')]).location]};
            %title(titleString,'FontSize',figFont+1,'Interpreter','Latex','FontName', 'Impact')
            title(titleString,'FontSize',figFont+1,'Interpreter','none')%'FontName', 'Impact')
            %$\textrm{t}$
        end
        %legend(ax(1),Path_File,'Interpreter','none');
        xlabel(['Sound level [Pa RMS]'],'FontSize',figFont,'Interpreter','Latex')

        set(fig1,'Position', [236 65 950 790],'Color','w');
        tightfig(fig1);

        % If requested, also save figures
        switch options.PlotSave
            case 1
                saveas(fig1, [filenameStem_SavedData '/' filenameFile_SavedData], 'fig')
                export_fig(fig1, [filenameStem_SavedData '/' filenameFile_SavedData '.pdf']);
                %export_fig(fig1, [filenameStem_SavedData '/' filenameFile_SavedData '.png'],'-m3');
                %export_fig(fig1, [filenameStem_SavedData '/' filenameFile_SavedData '.eps']);
        end
end

% Speak to the user via terminal
switch options.PlotSpec
    case 1
        switch options.PlotSpec_Range
            case 'all'
                disp('Note: Multiple figures have been plotted which show the raw spectra for all analysis frames used to obtain NSC values')
            case 'peak'
                disp('Note: Multiple figures have been plotted which show the raw spectra used to obtain NSC values for the frame corresponding to peak playing level')
            case 'min'
                disp('Note: Multiple figures have been plotted which show the raw spectra used to obtain NSC values for the frame corresponding to minimum playing level')
        end
end

end