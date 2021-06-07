function ROIdata_corr = sbxSlideFixROIs(savepath_ROI, savepath_editedROIstack, varargin)
% sbxSlideFixROIs use user-edited ROI stacks to fix the ROIs
% ROIdata_corr = sbxSlideFixROIs(savepath_ROI, savepath_editedROIstack, varargin)
% By Stephen Zhang 2019/05/20

% Default area threshold to be considered a region
p = inputParser;
addOptional(p, 'Region_size_thresh', 500); 
addOptional(p, 'SizeThresh', 14); 
parse(p, varargin{:});
p = p.Results;

% Load the ROI file
ROIdata = load(savepath_ROI, 'ROIdata');
ROIdata = ROIdata.ROIdata;

% Load the RGB stack
stack_corr = imread(savepath_editedROIstack);

% Process the original stack
stack_og = uint8(ROIdata.rgbstack * 255);

% Areas to subtract
areas2subtract = stack_og - stack_corr;
areas2subtract = areas2subtract(:,:,2) > 0;

% Areas to add
areas2add = stack_corr - stack_og;
areas2add_thresh = areas2add(:,:,2) >= 1;
[areas2add_labeled, nareas] = bwlabel(areas2add_thresh);

% Parse through the areas
areas2add_small = zeros(size(stack_og(:,:,1)));
regions = zeros(size(stack_og(:,:,1)));

% Number of regions and new areas
n_new_area = 0;
n_regions = 0;

for i = 1 : nareas
    % Current area
    Current_area = imfill(areas2add_labeled == i, 'holes');
    
    % Determine if a new area is a region or not
    if sum(Current_area(:)) >= p.Region_size_thresh
        % Region
        n_regions = n_regions + 1;
        regions = regions + Current_area * n_regions;   
    else
        % Area
        n_new_area = n_new_area + 1;
        areas2add_small = areas2add_small + Current_area * n_new_area;
    end
end

% Watershed new areas
D = bwdist(areas2add_small < 1);
D = -D;
D(areas2add_small < 1) = Inf;
L = watershed(D);
L(areas2add_small < 1) = 0;
L2 = areathresh(L, 1, p.SizeThresh, 3, 8);

% Relabel new areas
areas2add_small2 = bwlabel(L2);

% Copy the ROI data
ROIdata_corr = ROIdata;

% Add the areas
ROIdata_corr.ROIs(:,:,end+1) = areas2add_small2;

% Record the number of areas per slide
narea_vec = nan(size(ROIdata.ROIs,3), 1);

% Loop through and remove the areas
for i = 1 : size(ROIdata_corr.ROIs,3)
    % Remove areas
    currentframe = ROIdata_corr.ROIs(:,:,i);
    currentframe(areas2subtract) = 0;
    ROIdata_corr.ROIs(:,:,i) = currentframe; 
    
    % Count the number of areas per slide
    narea_vec(i) = max(currentframe(:));
end

% Update the cell counts
ROIdata_corr.ncells = sum(narea_vec);
ROIdata_corr.ncells_vec = narea_vec;

% Update regions
ROIdata_corr.regions = regions;
ROIdata_corr.n_regions = n_regions;

% Fix the RGB stack
ROIdata.rgbstack(:,:,1) = edge(sum(ROIdata_corr.ROIs,3) > 0);
end