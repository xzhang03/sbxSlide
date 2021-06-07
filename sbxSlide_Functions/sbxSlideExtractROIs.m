function icaguidata = sbxSlideExtractROIs(movpaths, savepath, edges, varargin)
%SBXEXTRACTROIS Extract ROIs from a series of SBX files (easy to add TIFF
%   capabilities) and save into the icaguidata format of a '.ica' file at
%   savepath. movpaths should be a cell array, savepath is a path to which
%   '.ica' will be appended if it does not exist

    p = inputParser;
    addOptional(p, 'pmt', 0, @isnumeric);  % Which PMT to use for analysis, 0-green, 1-red
    addOptional(p, 'type', 'pcaica');  % Can be 'nmf' or 'pcaica'. Will eventually add the ability for 'both'
    addOptional(p, 'force', false);  % Overwrite a save file if it already exists
    addOptional(p, 'axons', false);  % Set to true if running on axon data. Will be forced into PCA/ICA
    addOptional(p, 'downsample_t', 5);  % The number of pixels to downsample in time/frames
    addOptional(p, 'downsample_xy', 2);  % The number of pixels to downsample in space
    addOptional(p, 'chunksize', 1000);  % The size of a parallel chunk to be read in

    % PCA/ICA only
    addOptional(p, 'npcs', 1000, @isnumeric);  % The number of principal components to keep
                                               % WARNING: Divided by 4 for axons
    addOptional(p, 'firstpctokeep', 4, @isnumeric);
    addOptional(p, 'badframes', []);
    addOptional(p, 'temporal_weight', 0.1, @isnumeric);  % The temporal (versus spatial) weight for PCA. Default from Schnitzer was 0.5
    addOptional(p, 'smoothing_width', 2, @isnumeric);  % Standard deviation of Gaussian smoothing kernel (pixels)
                                                       % WARNING: Divided by 2 for axons
    addOptional(p, 'spatial_threshold_sd', 2, @isnumeric);  % Threshold for the spatial filters in standard deviations. Default from Schnitzer was 5
    
    % NMF only
    addOptional(p, 'cellhalfwidth', 2.5, @isnumeric);  % The half-width of cells in pixels
    addOptional(p, 'mergethreshold', 0.8, @isnumeric);  % The threshold for merging two ROIs that are neighboring
    addOptional(p, 'patchsize', -1, @isnumeric);  % Size of a patch to examine in parallel in pixels, calculated from cellhalfwidth
    addOptional(p, 'ncomponents', -1, @isnumeric);  % The number of rois to find in a patch, for advanced users only
    addOptional(p, 'minarea', 7, @isnumeric);  % The min area of a cell in pixels, will be pushed to 25 if left at 7 for PCA/ICA
    addOptional(p, 'maxarea', 500, @isnumeric);  % The max area of a cell in pixels, ignored for PCAICA unless different from 500
    addOptional(p, 'seeds', []);  % Seed locations for cells, replaces greedy algorithm
    
    % Parameters for both
    addOptional(p, 'overlap', 0.9, @isnumeric);  % The fraction of overlap allowed, combined with crosscorr
    addOptional(p, 'crosscorr', 0.9, @isnumeric);  % The fraction of correlation allowed, combined with overlap. The lower SNR of overlapping ROIs is removed
    addOptional(p, 'suffix', '.ica');  % Suffix to save into
    
    % Parameters for grouped std
    addOptional(p, 'grouopedstd', true); % Also calculate grouped std (grouped along the z axis)
    addOptional(p, 'groupedstd_stepsz', 10); % step size for grouped std
    
    parse(p, varargin{:});
    p = p.Results;
    
    if ~strcmp(savepath(end-3:end), p.suffix), savepath = [savepath p.suffix]; end
    if ~p.force && exist(savepath, 'file')
        fprintf('icaquidata already exists: %s\n', savepath);
        return;
    end
    if p.axons, p.type = 'pcaica'; end
    if strcmpi(p.type(1:3), 'pca') && p.minarea == 7, p.minarea = 25; end
    
    %% Read in all files and combine them into a single movie
    
    if isnumeric(movpaths)
        mov = movpaths;
    else

        info = sbxInfo(movpaths);
        nframes = info.max_idx + 1;
        optosteps = length(info.otwave);

        % Just read the movie
        fprintf('Reading movie...');
        mov = sbxReadPMT(movpaths, 0, nframes, p.pmt);
        fprintf(' Done.\n');

        % Rearrange frames
        % Throw out first frame and copy the top of the first volume to
        % the end
        mov = mov(:,:,[2:nframes, optosteps]);

        % Binning
        fprintf('Processing movie stack...')
        if p.downsample_xy > 1
            mov = binxy(mov, p.downsample_xy);
        end
        if p.downsample_t
            mov = bint(mov, p.downsample_t);
        end
        fprintf(' Done.\n')
    end
    
    %% Grouped std
    % Show
    fprintf('Grouped std...')
    
    % optotune steps and wave
    noptowaves = nframes / optosteps;
    optowave = repmat(1 : p.downsample_t : optosteps, [1 noptowaves]);
    
    % Grouped std
    ngroups = ceil(optosteps / p.groupedstd_stepsz);
    groupedstd_mat = (1 : p.groupedstd_stepsz : optosteps)';
    groupedstd_mat(:,2) = (p.groupedstd_stepsz : p.groupedstd_stepsz : optosteps)';
    
    % Grouped waves
    Grouped_optowave = repmat(optowave, [ngroups, 1]);
    
   for i = 1 : ngroups
       Grouped_optowave(i,:) = (Grouped_optowave(i,:) >= groupedstd_mat(i,1)) .*...
           (Grouped_optowave(i,:) <= groupedstd_mat(i,2));
   end
   
   % Grouped std frames
   Grouped_std_stack = zeros(info.sz(1), info.sz(2), ngroups);
   for i = 1 : ngroups
       Grouped_std_stack(:,:,i) = std(mov(:,:,Grouped_optowave(i,:) > 0), [], 3);
   end
   fprintf(' Done.\n')
   
    
    %% Extract ROIs
    fprintf('Extracting ROIs...\n')
    if strcmp(p.type, 'pcaica')
       
        maxarea = [];
        if p.maxarea ~= 500, maxarea = p.maxarea; end
        icaguidata = sbxSlidePreprocessPCAICA(mov, 'axons', p.axons, ...
            'npcs', p.npcs, 'temporal_weight', p.temporal_weight, 'smoothing_width', p.smoothing_width, ...
            'spatial_threshold_sd', p.spatial_threshold_sd, 'minarea', p.minarea, 'maxarea', maxarea, ...
            'overlap', p.overlap, 'crosscorr', p.crosscorr, 'firstpctokeep', p.firstpctokeep,...
            'badframes', p.badframes);
    end
    
    p.edges = edges;
    icaguidata.pars = p;
    fprintf('Done.\n')
    
    % Load up mov std
    icaguidata.movstd_grouped = Grouped_std_stack;
end

