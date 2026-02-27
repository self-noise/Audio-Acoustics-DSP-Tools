function MJN_CrescendoDecrescendo_CompiledPlotting(options,paths)
%---------------------------------------------------------------------------------------------%
% FUNCTION NAME AND SPECIFICATION
%
%           MJN_CrescendoDecrescendo_CompiledPlotting(options,paths)
%---------------------------------------------------------------------------------------------%
% PURPOSE OF THIS FUNCTION:
%           (1) Loads crescendo/decrescendo SC vs sigLev data files and plots the requested set
%           of data from within this/these file(s)
%---------------------------------------------------------------------------------------------%
% INPUTS:
%           options                     Structure containing all options to control the plot
%               .optionPlotAvg       
%               .optionErrorPlot     
%               .sigLev_delta        
%               .optionSaveFigure    
%               .N_Mics              
%               .listInstruments    
%               .listPlayers         
%               .SC_PlotLim          
%               .pressure_PlotLim    
%               .figFont             
%               .plotLineWidth       
%               .plotMarkerSize      
%               .plotMarkerStyleVec  
%               .plotMarkerSizeVec   
%               .listSEP             
%
%           paths                       Structure containing all paths to data files to load
%               .currentFolder
%               .dataFolder
%               .files
%               .outputFilename
%
% OUTPUTS:
%           A nicely formatted plot containing the requested SC vs sigLev data
%---------------------------------------------------------------------------------------------%
% GENERAL USAGE NOTES:
%   NOTE 1:
%---------------------------------------------------------------------------------------------%
% CHANGES TO ADD AT SOME POINT IN THE FUTURE:
%   TBC
%---------------------------------------------------------------------------------------------%
% CHANGELOG:
%
%   2024-02-01: Initial coding (via transfer from what was previously a script coped into each
%               experiment's folder
%---------------------------------------------------------------------------------------------%
% Author:           Dr Mike Newton
% Date:             December 2024
% Location (local): [Matlab_root]/LIBRARY/MJN_Code_Library/MJN_DSP_git/
% GitHub location:  https://github.com/self-noise/Audio-Acoustics-DSP-Tools
%---------------------------------------------------------------------------------------------%


%---------------------------------------------------------------------------------------------%
% Rename all the option entries (so we can keep syntax in rest of this function clean and clear
%---------------------------------------------------------------------------------------------%
optionPlotAvg       = options.optionPlotAvg;
optionErrorPlot     = options.optionErrorPlot;
sigLev_delta        = options.sigLev_delta;
optionSaveFigure    = options.optionSaveFigure;
N_Mics              = options.N_Mics;
listInstruments     = options.listInstruments;
plotDisplayNames    = options.plotDisplayNames;
listPlayers         = options.listPlayers;
SC_PlotLim          = options.SC_PlotLim;
pressure_PlotLim    = options.pressure_PlotLim;
figFont             = options.figFont;
plotLineWidth       = options.plotLineWidth;
plotMarkerSize      = options.plotMarkerSize;
plotMarkerStyleVec  = options.plotMarkerStyleVec;
plotMarkerSizeVec   = options.plotMarkerSizeVec;
listSEP             = options.listSEP;

%---------------------------------------------------------------------------------------------%
% Load datasets so they are easily accessible (datasets associated with 'listPlayers')
%---------------------------------------------------------------------------------------------%
N_DataSets      = size(paths.files,1);            % Total number of datasets we will load
for nDataSet = 1:N_DataSets
    dataSets{nDataSet,1} = load([paths.currentFolder '/' paths.dataFolder '/' paths.files{nDataSet}]);
end

%---------------------------------------------------------------------------------------------%
% Main plotting code
%---------------------------------------------------------------------------------------------%
N_Instruments   = size(options.listInstruments,1);  % Number of instruments we will plot
N_Players       = size(options.listPlayers,1);      % Number of individual players we will plot

% Create a list of nice colours, and rotate throuh the list if required (to permit unique 
% colour combinations)
colVec = MJN_Colours;
colVec.col = circshift(colVec.col,[0,-options.startColourIncrement]);
fig1 = figure(1);

% Each column is a different microphone. Each row is a different player. Optionally add a row
% for showing averages across all players at a given microphone.
switch optionPlotAvg
    case 0
        % There are as many rows in the plotted figure as players
        N_Players_TOT = N_Players;
        T_Layout = tiledlayout(N_Players_TOT,N_Mics,'padding','compact','TileSpacing','compact');
    case 1
        % There are as many rows in the plotted figure as players + 1 (for averages)
        N_Players_TOT = N_Players+1;
        T_Layout = tiledlayout(N_Players_TOT,N_Mics,'padding','compact','TileSpacing','compact');
    otherwise
        error('Error: "optionPlotAvg" can only be set to 0 or 1')
end

% Calculate the total number of 'subplot' tiles in the main plotting figure
countTotal      = N_Players_TOT*N_Mics;

% Main plotting loop
% Loop indices cycle through:
%           PLAYER-->MICROPHONE-->INSTRUMENT-->DATASET TYPE (crescendo and/or decrescendo)
count = 0;
for nPlayer = 1:N_Players_TOT
    for nMic = 1:N_Mics
        count = count+1;

        ax(count)           = nexttile(count);
        ax(count).FontSize  = figFont-2;
        hold on
        grid on

        countLeg = 0;
        for nInstrument=1:N_Instruments

            % Pull out the current player+instrument+dataset and plot it
            for nDataSet = 1:N_DataSets
                %-----------------------------------------------------------------------------%
                % This switch loop is key in determining which data gets pulled out prepared
                % for plotting
                %-----------------------------------------------------------------------------%
                countLeg = countLeg + 1;
                switch optionPlotAvg
                    case 0
                        % We aren't planning to do any per-instrument averaging across players,
                        % so just pull out data for the current player
                        % countLeg = countLeg + 1;
                        theseRows       = dataSets{nDataSet}.TableData(strcmp(dataSets{nDataSet}.TableData.Instrument, string(listInstruments{nInstrument,:})) & strcmp(dataSets{nDataSet}.TableData.Player, string(listPlayers{nPlayer,:})), :);
                    case 1
                        if nPlayer < N_Players_TOT
                            % We will be plotting the data averaged by instrument (i.e., across
                            % all players), but for the moment just load the relevant player's
                            % data as we haven't plotted all of it yet
                            % countLeg = countLeg + 1;
                            theseRows   = dataSets{nDataSet}.TableData(strcmp(dataSets{nDataSet}.TableData.Instrument, string(listInstruments{nInstrument,:})) & strcmp(dataSets{nDataSet}.TableData.Player, string(listPlayers{nPlayer,:})), :);

                        elseif nPlayer == N_Players_TOT
                            % All individual player data had been plotted, so now pull out the
                            % data for the current instrument and ALL players in a single dump,
                            % ready to be averaged
                            theseRows   = dataSets{nDataSet}.TableData(strcmp(dataSets{nDataSet}.TableData.Instrument, string(listInstruments{nInstrument,:})), :);
                        end
                end

                % Stack together the data into single matrices for sound pressure and SC (which
                % can then be averaged, analysed, etc)
                thisConcat_sigLev_Pa  = vertcat(theseRows.sigLev_Pa{:});
                thisConcat_SC_norm    = vertcat(theseRows.SC_norm{:});

                % thisLegendLabel{countLeg,:} = [char(listInstruments{nInstrument,:}) ' (' dataSets{nDataSet}.TableData.Cresc_or_Decresc{1} ')'];
                plotHandle(countLeg)        = plot(thisConcat_sigLev_Pa(:,1),thisConcat_SC_norm(:,nMic),plotMarkerStyleVec{nDataSet},'Parent', ax(count),'MarkerSize',plotMarkerSizeVec(nDataSet),'LineWidth', plotLineWidth,'MarkerEdgeColor',colVec.col{nInstrument},'MarkerFaceColor',colVec.col{nInstrument});
                plotHandle(countLeg).Marker = 'none';

                % On the first run through, make empty arrays to store concatenated datasets
                % together in the same place. This is relevant when multiple datasets are being
                % loaded together, such as a set of data for crescendos, and a separate set for
                % decresencos. This part of the code allows all crescendos+decrescendos for a
                % given combinations of player+microphone+instrument to be combined and then
                % analysed and plotted together
                switch nDataSet
                    case 1
                        thisConcat_sigLev_Pa_RETAIN = [];
                        thisConcat_SC_norm_RETAIN   = [];
                end

                switch optionErrorPlot
                    case 1
                        % Append current dataset to master 'retained' set
                        thisConcat_sigLev_Pa_RETAIN = [thisConcat_sigLev_Pa_RETAIN;thisConcat_sigLev_Pa];
                        thisConcat_SC_norm_RETAIN   = [thisConcat_SC_norm_RETAIN;thisConcat_SC_norm];

                        % Only upon loading the final dataset do we do the stats work
                        switch nDataSet
                            case N_DataSets

                                % Compute stats via nested function
                                % [sigLev_average,SC_norm_average,SC_norm_std] = SC_Stats(thisConcat_sigLev_Pa,thisConcat_SC_norm,N_Mics,sigLev_delta);
                                [sigLev_average,SC_norm_average,SC_norm_std] = SC_Stats(thisConcat_sigLev_Pa_RETAIN,thisConcat_SC_norm_RETAIN,N_Mics,sigLev_delta);

                                % Use 'shadedErrorBar' to add them to the current axis
                                sss(countLeg) = shadedErrorBar(sigLev_average, SC_norm_average(:,nMic), SC_norm_std(:,nMic),'plotAxes',ax(count),'lineprops', {'color',colVec.col{nInstrument}});

                                % Set face and edge properties of the stats on the plot
                                %set(sss(1).edge,'LineWidth',2,'LineStyle',':')
                                sss(countLeg).mainLine.LineWidth = 3;
                                sss(countLeg).patch.FaceColor = colVec.col{nInstrument};

                                %thisLegendLabel_STATS{nInstrument,:} = [char(listInstruments{nInstrument,:}) ' (crescendo+decresendo combined)'];
                                
                                %thisLegendLabel_STATS{nInstrument,:} = [char(listInstruments{nInstrument,:})];
                                thisLegendLabel_STATS{nInstrument,:} = [char(plotDisplayNames{nInstrument,:})];
                        end
                end
            end

            xlim([0,pressure_PlotLim])
            ylim([1,SC_PlotLim])
        end

        switch optionPlotAvg
            case 0
                switch nPlayer
                    case N_Players
                        %xlabel({'RMS acoustic pressure (at'; 'microphone 1, 3m from player)'},'FontSize',figFont)                        
                        xlabel(['$\bar{p}_{\mathrm{(mic \,1)}}$'],'FontSize',figFont+6,'Interpreter','latex')                        
                end
            case 1
                switch nPlayer
                    case N_Players_TOT
                        % xlabel({'RMS acoustic pressure (at'; 'microphone 1, 3m from player)'},'FontSize',figFont)                        
                        xlabel(['$\bar{p}_{\mathrm{(mic \,1)}}$'],'FontSize',figFont+6,'Interpreter','latex')
                end

        end

        switch nMic
            case 1
                %ylabel({'Normalised';'spectral centroid'},'FontSize',figFont)
                ylabel(['$\hat{S}_c$'],'FontSize',figFont+6,'Interpreter','latex')
        end

        switch optionPlotAvg
            case 0
                % title(ax(count),['Mic ' num2str(nMic) ' (' char(listPlayers{nPlayer,:}) ')'],'FontSize',figFont)
                title(ax(count),['Mic ' num2str(nMic) ' (player: ' char(listPlayers{nPlayer,:}) ')'],'FontSize',figFont)
            case 1
                if nPlayer < N_Players_TOT
                    % title(ax(count),['Mic ' num2str(nMic) ' (' char(listPlayers{nPlayer,:}) ')'],'FontSize',figFont)
                    title(ax(count),['Mic ' num2str(nMic) ' (player: ' char(listPlayers{nPlayer,:}) ')'],'FontSize',figFont)
                else
                    %title(ax(count),['Mic ' num2str(nMic) ' (averaged over ' num2str() ' player)' ],'FontSize',figFont)
                    title(ax(count),['Mic ' num2str(nMic) ' (average of all players)' ],'FontSize',figFont)
                end


        end
    end
end

% Associate the legend with first subplot tile (though all tiles will be the same labelling) and
% place it in a dedicated area. Also reorder the lines/labels so the legend appears in descending
% order of expected 'brassiness' (based on SEP)
% lh = legend(ax(1),plotHandle,thisLegendLabel,'Location','NorthEastOutside','Orientation','Vertical','FontSize',figFont);

% [listSEP_ordered,idx_SEP] = sortrows(listSEP);
statsLines = [sss.mainLine];
% statsLines = statsLines(flipud(idx_SEP));
% thisLegendLabel_STATS = thisLegendLabel_STATS(flipud(idx_SEP))
lh = legend(ax(1),statsLines,thisLegendLabel_STATS,'Location','NorthEastOutside','Orientation','Horizontal','FontSize',figFont,'Interpreter','None');
lh.Layout.Tile = 'North';

fig1.Color = 'w';



% Figure size and positioning
% [Lower-left-x, Lower-left-y, width (total), height (total)]
switch N_Players_TOT
    case 1
        fig1.Position = [433 1 915 1000*(1/4)*1.02];
    case 2
        fig1.Position = [433 1 915 1000*(2/4)*1.02];
    case 3
        fig1.Position = [433 1 915 1000*(3/4)*1.02];
    case 4
        fig1.Position = [433 1 915 1000];
end



% %---------------------------------------------------------------------------------------------%
% % SECONDARY PLOT: SPECTRAL ENRICHMENT PARAMETER FOR EACH INSTRUMENT
% %---------------------------------------------------------------------------------------------%
%
% SEP_dataSort1   = [listInstruments,num2cell(listSEP),colVec.col(1:N_Instruments)'];
% SEP_dataSort1   = sortrows(SEP_dataSort1,2);
% SEP_dataSort2   = [listInstruments2,num2cell(listSEP2),colVec.col(N_Instruments+1+1:N_Instruments_TOT+1)'];
% SEP_dataSort2   = sortrows(SEP_dataSort2,2);
% %SEP_dataSort    = [SEP_dataSort1;SEP_dataSort2];
% %SEP_dataSort    = sortrows(SEP_dataSort,2);
%
% fig2 = figure(2);
%
% %ax2(1) = subplot('Position', [0.1, 0.28, 0.485, 0.665]); % [left, bottom, width, height]
% ax2(1) = subplot('Position', [0.1, 0.28, 0.55, 0.665]); % [left, bottom, width, height]
% hold on
% grid on
% for nInstrument=1:N_Instruments
%     bar(nInstrument,SEP_dataSort1{nInstrument,2},'FaceColor',SEP_dataSort1{nInstrument,3},'EdgeColor',SEP_dataSort1{nInstrument,3})
% end
%
% ylabel('Spectral enrichment parameter (E)','FontSize',figFont+4)
% title(["8' instruments"],'FontSize',figFont+4)
% xlim([0.5,N_Instruments+0.5])
% ylim([0,8])
% xticks((1:N_Instruments))
% xticklabels(SEP_dataSort1(:,1))
% ax2(1).FontSize = figFont+2;
%
%
%
% %ax2(2) = subplot('Position', [0.66, 0.28, 0.3, 0.665]); % Adjust left and width to fit
% ax2(2) = subplot('Position', [0.725, 0.28, 0.22, 0.665]); % Adjust left and width to fit
% hold on
% grid on
% for nInstrument=1:N_Instruments2
%     bar(nInstrument,SEP_dataSort2{nInstrument,2},'FaceColor',SEP_dataSort2{nInstrument,3},'EdgeColor',SEP_dataSort2{nInstrument,3})
% end
%
% title(["Horns"],'FontSize',figFont+4)
% xlim([0.5,N_Instruments2+0.5])
% ylim([0,8])
% xticks((1:N_Instruments2))
% xticklabels(SEP_dataSort2(:,1))
% ax2(2).FontSize = figFont+2;
%
%
% fig2.Color = 'w';
% fig2.Position = [655 816 1009 582];




%---------------------------------------------------------------------------------------------%
% OPTIONAL: Figure export/saving
%---------------------------------------------------------------------------------------------%
switch optionSaveFigure
    case 1
        export_fig(fig1, paths.outputFilename,'-pdf','-fig');
        % export_fig(fig2, 'All players -- All instruments -- Crescendo + Decrescendo -- SEP (shaded STD V3)' ,'-pdf','-fig');
end



%---------------------------------------------------------------------------------------------%
% Nested functions below this line
%---------------------------------------------------------------------------------------------%
    function [sigLev_average,SC_norm_average,SC_norm_std] = SC_Stats(thisConcat_sigLev_Pa,thisConcat_SC_norm,N_Mics,sigLev_delta)

        % Sort rows by ascending sigLev (mic01)
        [thisConcat_sigLev_Pa,idx]  = sortrows(thisConcat_sigLev_Pa);
        thisConcat_SC_norm          = thisConcat_SC_norm(idx,:);

        % Now window the signal level and SC data by discrete steps (set via sigLev_delta) and perform
        % basic statistics on each window of data

        % What is the total range of sigLev (mic01)?
        sigLev_range = thisConcat_sigLev_Pa(end,1) - thisConcat_sigLev_Pa(1,1);

        % How many discrete windows do we need in order to cover the total sigLev (mic01) range?

        %thisConcat_sigLev_Pa(1,1)
        sigLev_lowerRounded = floor(thisConcat_sigLev_Pa(1,1)/sigLev_delta)*sigLev_delta;


        %thisConcat_sigLev_Pa(end,1)
        sigLev_upperRounded = ceil(thisConcat_sigLev_Pa(end,1)/sigLev_delta)*sigLev_delta;


        sigLev_windowVec = sigLev_lowerRounded:sigLev_delta:sigLev_upperRounded;
        N_sigLev_windows = length(sigLev_windowVec)-1;

        sigLev_average  = zeros(N_sigLev_windows,1);
        SC_norm_average = zeros(N_sigLev_windows,N_Mics);
        SC_norm_std     = zeros(N_sigLev_windows,N_Mics);
        for nWin=1:N_sigLev_windows
            % Find all measurements in the current sigLev window range
            indices = find(thisConcat_sigLev_Pa(:,1) >= sigLev_windowVec(nWin) & thisConcat_sigLev_Pa(:,1) < sigLev_windowVec(nWin+1));

            % Using the relevant indices, average together the corresponding rows of SC data,
            % preserving it by microphone (i.e., column 1 indices are averaged to produce a single
            % number, column 2 indices are averaged to produce a single number, etc)

            sigLev_average(nWin)    = mean(thisConcat_sigLev_Pa(indices,1),1);

            SC_norm_average(nWin,:) = mean(thisConcat_SC_norm(indices,:),1);
            SC_norm_std(nWin,:)     = std(thisConcat_SC_norm(indices,:),0,1);
            % SC_norm_std(nWin,:)     = iqr(thisConcat_SC_norm(indices,:),1);

        end
    end

end

