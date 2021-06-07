% By Stephen Zhang 2019/05/20
%% General variables
clear

% Force overwrite?
FORCE_OVERWRITE = true;

% Edges (place holder)
edges = [30,60,10,10];
%% Get file info
% Get directory
default_add = '\\nasquatch\data\2p\stephen';
fpath = uigetdir(default_add);

flist = dir(fullfile(fpath,'**','*.sbx'));
nfiles = length(flist);

% Only process channel 1
flist2keep = ones(nfiles,1);
for i = 1 : nfiles
    % Only keep odd runs
    flist2keep(i) = mod(flist(i).name(end-4),2);    
end

flist =  flist(flist2keep == 1);
nfiles = sum(flist2keep);

%% Loop through
for file_ind = 1 : nfiles
    %% Extraction
    disp(['-------------- Processing File ', num2str(file_ind), ' out of ', num2str(nfiles), '. --------------'])
    
    % Paths
    movpath = fullfile(flist(file_ind).folder,flist(file_ind).name);
    savepath_ica = sprintf('%s_ica.mat', movpath(1:end-4));
    if ~exist(savepath_ica, 'file') || FORCE_OVERWRITE
        icaguidata = sbxSlideExtractROIs(movpath, savepath_ica, edges, 'type', 'pcaica', 'force', true, ...
            'axons', false, 'downsample_t', 2, 'downsample_xy', 1, ...
            'chunksize', 1000, 'npcs',50, 'temporal_weight', 0.99, ...
            'smoothing_width', 1.2, 'spatial_threshold_sd', 2,'firstpctokeep',1);

    
        save(savepath_ica, 'icaguidata', '-v7.3');
    else
        disp('ICA data already exist and will be used.')
        icaguidata = load(savepath_ica, 'icaguidata');
        icaguidata = icaguidata.icaguidata;
    end
    %% Ica morphological filter
    savepath_ROI = sprintf('%s_ROI.mat', movpath(1:end-4));
    
    if ~exist(savepath_ROI, 'file') || FORCE_OVERWRITE
        % Use ICA/std based morphological filter
        ROIdata = sbxSlideICAMorphFilter(icaguidata, 'PlotOrNot', false, 'SizeThresh', [14 200],...
            'MorphThresh', 0.25, 'RatioThresh', -1, 'TwoD_RedundencyRemoval', false, 'UseGroupedSTDImage',...
            true, 'IterativeSubtraction', false, 'FlattenSTDImage', true);
        
        % Apply filters
        ROIdata.signals = sbxSlideApplyFilters(ROIdata, movpath, '', 'ProgressBar', true,...
            'ProcessRegions', false, 'makefilterRGB', true, 'DoNeuropil', true, 'NeuropilRingSize', 1,...
            'MinNeuropilPixels', 3);
        
        % Filter based on signal and neuropil signal
        ROIdata = sbxSlideROIFilterByNp(ROIdata, 'Method', 'LocalExtremaToMedian', 'Threshold', 2000);
%         test = ROIdata.SigPass.rgbstack;
%         test(:,:,3) = 0;
%         imtool(test)
        
        % Area overlap and correlation redundancy removal
        ROIdata = sbxSlideROIFilterByAreaCorrPeaks(ROIdata, 'AreaOverlap', 0.5,...
            'PeakProminenceThreshold', 1000);
%         test = ROIdata.AreaCorrPass.rgb_stack;
%         test(:,:,3) = 0;
%         imtool(test)

        % Write ROIs
        ROIdata = sbxSlideROIWrite(ROIdata, icaguidata, movpath, 'UseGroupedSTDImage', true,...
            'NtoCombine', 18, 'UseROIEdge', true);
        
        % Saving stuff
        if ~exist(savepath_ROI, 'file') || FORCE_OVERWRITE
            save(savepath_ROI, 'ROIdata', '-v7.3');
        end

    end

    
end
