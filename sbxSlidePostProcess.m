% By Stephen Zhang 2019/05/20
%% General variables
clear

% Force overwrite?
FORCE_OVERWRITE = false;

%% Get file info
% Get directory
default_add = '\\anastasia\data\2p\stephen';
fpath = uigetdir(default_add);

flist = dir(fullfile(fpath,'**','*-1.tif'));
nfiles = length(flist);

% Only process channel 1
flist2keep = ones(nfiles,1);
for file_ind = 1 : nfiles
    % Only keep odd runs
    flist2keep(file_ind) = mod(flist(file_ind).name(end-4),2);    
end

flist =  flist(flist2keep == 1);
nfiles = sum(flist2keep);

%% Loop through
for file_ind = 1 : nfiles
    % Display
    disp(['-------------- Processing File ', num2str(file_ind), ' out of ', num2str(nfiles), '. --------------'])
    
    % Parse the file paths
    % Edited Tif address
    savepath_editedROIstack = fullfile(flist(file_ind).folder,flist(file_ind).name);
    
    % ROI address, original and corrected
    savepath_ROI = sprintf('%s.mat',savepath_editedROIstack(1:end-10));
    savepath_ROI_corr = [savepath_ROI(1:end-4), '_corr.mat'];
    
    % Sbx addresses
    savepath_sbx_Ch1 = sprintf('%s.sbx',savepath_editedROIstack(1:end-14));
    savepath_sbx_Ch2 = [savepath_sbx_Ch1(1:end-5),num2str(str2double(savepath_sbx_Ch1(end-4))+1),'.sbx'];
    
    % Signal address
    savepath_sig = sprintf('%s_sig.mat',savepath_editedROIstack(1:end-14));
    
    % General info address
    savepath_geninfo = sprintf('%s_GenInfo.mat',savepath_editedROIstack(1:end-14));
    
    % Results address
    savepath_results = sprintf('%s_results.mat',savepath_editedROIstack(1:end-14));
    
    %% Use user file to fix ROIs
    if ~exist(savepath_ROI_corr, 'file') || FORCE_OVERWRITE
        % If no file, fix and save
        % Fix the ROIs using user edits
        ROIdata_corr = sbxSlideFixROIs(savepath_ROI, savepath_editedROIstack, 'Region_size_thresh', 500, 'SizeThresh', 14);
        
        % Save the fixed ROIs
        save(savepath_ROI_corr, 'ROIdata_corr', '-v7.3');
    else
        % If the fixed ROIs mat file exst, use it
        ROIdata_corr = load(savepath_ROI_corr);
        ROIdata_corr = ROIdata_corr.ROIdata_corr;
    end
    
    %% Apply ROIs to both channels
    if ~exist(savepath_sig, 'file') || ~exist(savepath_geninfo, 'file') || FORCE_OVERWRITE
        % If no file, apply filters and save
        [ROIsignal, GenInfo] = sbxSlideApplyFilters(ROIdata_corr, savepath_sbx_Ch1, savepath_sbx_Ch2, 'gapsize', 10);
        
        % Save the signals
        save(savepath_sig, 'ROIsignal', '-v7.3');
        save(savepath_geninfo, 'GenInfo', '-v7.3');
    else
        % If the fixed ROIs mat file exst, use it
        ROIsignal = load(savepath_sig);
        ROIsignal = ROIsignal.ROIsignal;
        
        GenInfo = load(savepath_geninfo);
        GenInfo = GenInfo.GenInfo;
    end
    
    %% Determine colocalization
    if ~exist(savepath_results, 'file') || FORCE_OVERWRITE
        % If no file, apply colocalization and save
        ROIResults = sbxSlideColocal(ROIsignal, GenInfo, 'HMMode', 'LastPointAboveHM', 'Tolerance', 1);
        
        % Save the results
        save(savepath_results, 'ROIResults', '-v7.3');
        
        % Generate RGB stacks of the results
        [rgb_colocal, rgb_ortho] = sbxSlideMakeRGBfromROI(ROIdata_corr, ROIsignal, GenInfo, ROIResults);
        
        % Write the images
        imwrite(rgb_colocal, sprintf('%s_colocal.tif',savepath_results(1:end-4)));
        imwrite(rgb_ortho, sprintf('%s_orthogonal.tif',savepath_results(1:end-4)));
        
        % Generate csv output
        CSV_table = sbxSlideResultTable(ROIResults, ROIsignal);
        
        % WRite outputtable
        writetable(CSV_table, sprintf('%s.csv',savepath_results(1:end-4)));
    else
        % If the fixed ROIs mat file exst, use it
        ROIResults = load(savepath_results);
        ROIResults = ROIResults.ROIResults;
    end
    
end

% Say done
disp(['Done.'])