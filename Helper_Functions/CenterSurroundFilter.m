function icaoutput = CenterSurroundFilter(inputim, icainput, varargin)
% CenterSurroundFilter apply an object/surround threshold to remove ROIs
% that are too dark
% icaoutput = CenterSurroundFilter(inputim, icainput, varargin)

p = inputParser;
addOptional(p, 'RatioThresh', 1.3); % Default 1.3 as the center/surround threshold
addOptional(p, 'PlotOrNot', false); % Default don't plot anything
addOptional(p, 'PlotGrid', [10 10]); % Default to plot 10x10 grid of segmentations, if PlotOrNot is true.
addOptional(p, 'npix', 11);  % Default 11 pixel wide squares
addOptional(p, 'SizeThresh', [14 200]); % Size threshold
addOptional(p, 'subtractMask', []); % Default don't subtract anything from the mask

if length(varargin) == 1 && iscell(varargin{1}), varargin = varargin{1}; end
parse(p, varargin{:});
p = p.Results;

% Radius of the square
npixr = floor(p.npix/2);

% Initiate figure if needed
if p.PlotOrNot
    figure
end

% Initalize book-keeping variables
ncells = 0;
FlagPass = zeros(length(icainput.ica),1) == 1;
totalFilter = zeros(size(icainput.AllFilters));
Sizes = zeros(length(icainput.ica),1);

% Loop through to apply the center/surround threshold as well as the size
% threshold again
for i = 1 : length(icainput.ica)
    
    % Load in basic variables
    CurrFilt = icainput.ica(i).filter > 0;
    centroid = icainput.Centroids(i,:);
    
    % Remove previous filters
    if ~isempty(p.subtractMask)
        CurrFilt = CurrFilt & (p.subtractMask == 0);
    end
    
    % Draw a predefined square around the centroid (but not go over the bounts)
    im_rlim = max(centroid(1)-npixr : centroid(1)+npixr, 1);
    im_rlim = min(im_rlim, size(inputim,1));
    im_clim = max(centroid(2)-npixr : centroid(2)+npixr, 1);
    im_clim = min(im_clim, size(inputim,2));
    
    % Crop the image
    inputim_crop = inputim(im_rlim, im_clim);
    CurrFilt_crop = CurrFilt(im_rlim, im_clim);
    
    % Calculate center, surround and size
    v_cen = median(inputim_crop(CurrFilt_crop));
    v_surr = median(inputim_crop(~CurrFilt_crop));
    Sizes(i) = sum(CurrFilt(:));
    
    if v_cen/v_surr > p.RatioThresh && Sizes(i) >= p.SizeThresh(1)
        % Propagage the number of cells
        ncells = ncells + 1;
        
        % Flag
        FlagPass(i) = true;
        
        % Update total filter
        totalFilter = totalFilter + CurrFilt * ncells;
                
        % Make plot if called for
        if ncells <= p.PlotGrid(1) * p.PlotGrid(2) && p.PlotOrNot
            % Calculate a croped, filtered image
            inputim_crop_filt = inputim_crop .* CurrFilt_crop;
            subplot(p.PlotGrid(1),  p.PlotGrid(2), ncells);
            imshow(inputim_crop_filt,[]);
            text(0,-npixr,num2str(i));
        end
        
    end
    
end

% Update ica data
icaoutput = icainput;
icaoutput.AllFilters = totalFilter;
icaoutput.CellAreas = Sizes(FlagPass);
icaoutput.Centroids = icaoutput.Centroids(FlagPass,:);
icaoutput.ica = icaoutput.ica(FlagPass);
end