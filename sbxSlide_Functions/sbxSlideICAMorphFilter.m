function ROIdata = sbxSlideICAMorphFilter(icaguidata, varargin)
% SbxICAMorphFilter apply a morphological filter to extra ROIs from ICAs
% ROIdata = sbxICAMorphFilter(icaguidata, varargin)
% By Stephen Zhang 2019/05/20
p = inputParser;
addOptional(p, 'PlotOrNot', false);  % Plot the figure or not
addOptional(p, 'SizeThresh', [14 200]); % Minimal and maximal areas to be considered an ROI
addOptional(p, 'MorphThresh', 0.03); % Threshold for morphological filter
addOptional(p, 'RatioThresh', 1.4); % Ratio threshold
addOptional(p, 'TwoD_RedundencyRemoval', true); % Applying a method to remove filters based on 2D data
                                              % Good for thin slices
addOptional(p, 'IterativeSubtraction', true); % Iteratively subtract existing filters
                                              % Good for thin slices

% Std image
addOptional(p, 'UseSTDImage', true); % Use STD image of the volume to enhance contrast
addOptional(p, 'UseGroupedSTDImage', false); % Iteratively apply grouped STD images to enhance contrast (will increase the number of filter-frames). 
addOptional(p, 'FlattenSTDImage', false); % Locally normalize the STD image or not

parse(p, varargin{:});
p = p.Results;

%% Initialize all filters 
Allfilters_flat = zeros(size(icaguidata.movstd));

% Determine the number of filter frames
if p.UseGroupedSTDImage
    ngroupedfilters = size(icaguidata.movstd_grouped, 3);
    Allfilters_stack = ...
        repmat(Allfilters_flat, [1 1 size(icaguidata.ica_filters,3) * 2 * ngroupedfilters]);
    ncell_vec = zeros(size(icaguidata.ica_filters,3) * 2 * ngroupedfilters,1);
else
    Allfilters_stack = repmat(Allfilters_flat, [1 1 size(icaguidata.ica_filters,3) * 2]);
    ncell_vec = zeros(size(icaguidata.ica_filters , 3) * 2,1);
end

% Flatten grouped std images
if p.UseGroupedSTDImage && p.FlattenSTDImage
    n=8; %assumes neurons are around 8 pixels wide
    m=150; % Ring size
    
    for j = 1 : ngroupedfilters
        stdim = icaguidata.movstd_grouped(:,:,j);
        stdim_prime = double(stdim)-double(imgaussfilt(double(stdim),n));
        icaguidata.movstd_grouped(:,:,j) = ...
            stdim_prime ./ (imgaussfilt(stdim_prime.^2,m) .^ (1/2));
    end
end

%% Processing positive ICA filters
hwait = waitbar(0, 'Processing');
for i = 1 : size(icaguidata.ica_filters,3) * 2 % Loop through and alternate between positive and negative filters
    % Use Morphological filter
    % Update waitbar
    waitbar((i-1)/size(icaguidata.ica_filters,3)/2, hwait,....
        sprintf('Processing %i/%i', i, size(icaguidata.ica_filters,3) * 2))

    Morphvaragin = {'downsample_xy', 1, 'PlotOrNot', p.PlotOrNot, 'subtractMask', Allfilters_flat > 0,...
        'Threshold', p.MorphThresh, 'SizeThresh', p.SizeThresh - 1};

    % Alternate between positive and negative filters
    signfunc = (mod(i,2) - 0.5) / 0.5; % Sign function based on odd/even
    ica_ind = ceil(i/2); % Figure out which ica to segment
    
    if ~p.UseGroupedSTDImage % Using overall std image if any
        % Apply filter
        if p.UseSTDImage && p.FlattenSTDImage
            % Use dot product with std image to enhance contrast
            % Flatten
            stdim = icaguidata.movstd;
            n=8; %assumes neurons are around 8 pixels wide
            m=150; % Ring size
            stdim_prime = double(stdim)-double(imgaussfilt(double(stdim),n));
            stdim = stdim_prime ./ (imgaussfilt(stdim_prime.^2,m) .^ (1/2));
            CurrIm = signfunc *  stdim .* icaguidata.ica_filters(:,:,ica_ind);
        elseif p.UseSTDImage
            % Use dot product with std image to enhance contrast
            CurrIm = signfunc * icaguidata.ica_filters(:,:,ica_ind);
        end
        
        % Apply filter
        MorphFilters = sbxMorphologicalFilterExtractROIs_SZ(CurrIm, Morphvaragin);

        % Apply ratio filter (Center/Surround) and size filter
        % Check if any ROI is detected
        if isfield(MorphFilters, 'ica')
            if p.IterativeSubtraction
                RatioFiltervaragin = {'PlotOrNot', p.PlotOrNot, 'RatioThresh', p.RatioThresh, 'npix', 11,...
                    'PlotGrid', [10 10], 'SizeThresh', p.SizeThresh, 'subtractMask', Allfilters_flat > 0};
            else
                RatioFiltervaragin = {'PlotOrNot', p.PlotOrNot, 'RatioThresh', p.RatioThresh, 'npix', 11,...
                    'PlotGrid', [10 10], 'SizeThresh', p.SizeThresh, 'subtractMask', []};
            end
            MorphFilters = CenterSurroundFilter(icaguidata.movstd, MorphFilters, RatioFiltervaragin);

            % Update all filters
            % Update both flat and stack filters
            Allfilters_flat = Allfilters_flat + (MorphFilters.AllFilters > 0) * i;
            Allfilters_stack(:,:,i) = MorphFilters.AllFilters;

            % Update the cell numbers
            ncell_vec(i) = length(MorphFilters.ica);
        else
            Allfilters_stack(:,:,i) = 0;
            ncell_vec(i) = 0;

        end
    else % Use grouped std image
        for j = 1 : ngroupedfilters
            % Calculate an overall index
            ind = (i-1)*ngroupedfilters + j;
            
            % Use dot product with std image to enhance contrast
            % Current image
            CurrIm = signfunc * icaguidata.movstd_grouped(:,:,j)...
                .* icaguidata.ica_filters(:,:,ica_ind);
            MorphFilters = sbxMorphologicalFilterExtractROIs_SZ(CurrIm, Morphvaragin);
            
            % Apply ratio filter (Center/Surround) and size filter
            % Check if any ROI is detected
            if isfield(MorphFilters, 'ica')
                if p.IterativeSubtraction
                    RatioFiltervaragin = {'PlotOrNot', p.PlotOrNot, 'RatioThresh', p.RatioThresh, 'npix', 11,...
                        'PlotGrid', [10 10], 'SizeThresh', p.SizeThresh, 'subtractMask', Allfilters_flat > 0};
                else
                    RatioFiltervaragin = {'PlotOrNot', p.PlotOrNot, 'RatioThresh', p.RatioThresh, 'npix', 11,...
                        'PlotGrid', [10 10], 'SizeThresh', p.SizeThresh, 'subtractMask', []};
                end
                MorphFilters = CenterSurroundFilter(icaguidata.movstd_grouped(:,:,j), MorphFilters, RatioFiltervaragin);

                % Update all filters
                % Update both flat and stack filters
                Allfilters_flat = Allfilters_flat + (MorphFilters.AllFilters > 0) * ind;
                Allfilters_stack(:,:,ind) = MorphFilters.AllFilters;

                % Update the cell numbers
                ncell_vec(ind) = length(MorphFilters.ica);
            else
                Allfilters_stack(:,:,ind) = 0;
                ncell_vec(ind) = 0;

            end
        end
    end



end
close(hwait)
%% Filter out redundant ROIs
if p.TwoD_RedundencyRemoval
    Allfilters_stack2 = Allfilters_stack > 0;
    Allfilters_stack3 = cumsum(Allfilters_stack2,3) > 0;
    Allfilters_stack4 = Allfilters_stack;
    Allfilters_stack4(:,:,2:end) = (Allfilters_stack2(:,:,2:end) - Allfilters_stack3(:,:,1:end-1)) .* Allfilters_stack(:,:,2:end);

    % Initial the corrected filters
    Allfilters_stack_corr = zeros(size(Allfilters_stack2));
    ncell_vec_corr = zeros(size(Allfilters_stack2,3),1);

    for ii = 1 : size(Allfilters_stack2,3)
        currframe = Allfilters_stack4(:,:,ii);
        nROIs = max(currframe(:));

        for i = 1 : nROIs
            currROI = currframe == i;
            if sum(currROI(:)) >= p.SizeThresh(1)
                Allfilters_stack_corr(:,:,ii) = (currframe == i)*i + Allfilters_stack_corr(:,:,ii);
                ncell_vec_corr(ii) = ncell_vec_corr(ii) + 1;
            end
        end
    end
else
    Allfilters_stack_corr = Allfilters_stack;
    ncell_vec_corr = ncell_vec;
end

%% Make plot
% Plot all filters
rgbstack = repmat(mat2gray(icaguidata.movstd),[1 1 3]);
rgbstack(:,:,3) = mat2gray(icaguidata.movmax);
rgbstack(:,:,1) = edge(sum(Allfilters_stack_corr,3) > 0);

if p.PlotOrNot
    figure, imshow(rgbstack);
end

% Add up the cell numbers
ncells = sum(ncell_vec_corr);


%% Output
ROIdata.rgbstack = rgbstack;
ROIdata.ncells = ncells;
ROIdata.ncells_vec = ncell_vec_corr;
ROIdata.ROIs = Allfilters_stack_corr;
    
end