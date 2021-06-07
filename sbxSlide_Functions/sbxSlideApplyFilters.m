function [ROIsignal, GenInfo] = sbxSlideApplyFilters(ROIdata_corr, savepath_sbx_Ch1, savepath_sbx_Ch2, varargin)
% sbxSlideApplyFilter applies ROI filters to the original sbx files to
% calculate signal
% [ROIsignal, im_stds] = sbxSlideApplyFilters(ROIdata_corr, savepath_sbx_Ch1, savepath_sbx_Ch2, varargin)
% By Stephen Zhang 2019/05/20

% Gap size in generating the rgb
p = inputParser;

% Individual filter RGB stack
addOptional(p, 'makefilterRGB', true); % Make filter RGB
addOptional(p, 'gapsize', 10); % Padding for cell images

% General settings
addOptional(p, 'ProgressBar', false); % Default no progress bar
addOptional(p, 'ProcessRegions', true); % Categorize cells by their brain regions (user supplied)

% Pre-bin before calculations (not used for rgb stacks)
addOptional(p, 'binxy', 2);

% Neuropil
addOptional(p, 'DoNeuropil', false); % Calculate neuropil signal
addOptional(p, 'NeuropilRingSize', 2); % The first ring size ty (Will increase until the number of pixels suffice)
addOptional(p, 'MinNeuropilPixels', 10); % Minimal number of pixels to calculate neuropil signal

parse(p, varargin{:});
p = p.Results;

% Determine which channel to be done
Ch1_do = ~isempty(savepath_sbx_Ch1);
Ch2_do = ~isempty(savepath_sbx_Ch2);

% Get sbx file metainfo
if Ch1_do
    % Channel 1
    sbxstruct_1 = sbxInfo(savepath_sbx_Ch1);
    nframes_1 = sbxstruct_1.max_idx + 1;
end
if Ch2_do
    % Channel 2
    sbxstruct_2 = sbxInfo(savepath_sbx_Ch2);
    nframes_2 = sbxstruct_2.max_idx + 1;
end

% Read the movie
if Ch1_do
    fprintf('Reading Channel 1...')
    imdata_Ch1 = sbxReadPMT(savepath_sbx_Ch1, 0, nframes_1, 0);
    fprintf(' Done.\n')
end
if Ch2_do
    fprintf('Reading Channel 2...')
    imdata_Ch2 = sbxReadPMT(savepath_sbx_Ch2, 0, nframes_2, 0);
    fprintf(' Done.\n')
end

% Take optotune levels
if Ch1_do
    sweeps_1 = nframes_1/sbxstruct_1.otparam(3);
end
if Ch2_do
    sweeps_2 = nframes_2/sbxstruct_2.otparam(3);
end

% Take median between sweeps
% Rearrange (x y optotune sweep)
if Ch1_do
    imdata_Ch1_re = reshape(imdata_Ch1, [sbxstruct_1.sz(1), sbxstruct_1.sz(2),...
        sbxstruct_1.otparam(3), sweeps_1]);
end
if Ch2_do
    imdata_Ch2_re = reshape(imdata_Ch2, [sbxstruct_2.sz(1), sbxstruct_2.sz(2),...
        sbxstruct_2.otparam(3), sweeps_2]);
end

% Take median values between sweeps
if Ch1_do
    fprintf('Calculating median Channel 1...')
    imdata_Ch1_med = median(imdata_Ch1_re,4);
    fprintf(' Done.\n')
end
if Ch2_do
    fprintf('Calculating median Channel 2...')
    imdata_Ch2_med = median(imdata_Ch2_re,4);
    fprintf(' Done.\n')
end

% Rearrange the optotune order
% Channel 1
if Ch1_do
    imdata_Ch1_med2 = imdata_Ch1_med;
    imdata_Ch1_med2(:,:,1:end-1) = imdata_Ch1_med(:,:,2:end);
    imdata_Ch1_med2(:,:,end) = imdata_Ch1_med(:,:,1);
end

% Channel 2
if Ch2_do
    imdata_Ch2_med2 = imdata_Ch2_med;
    imdata_Ch2_med2(:,:,1:end-1) = imdata_Ch2_med(:,:,2:end);
    imdata_Ch2_med2(:,:,end) = imdata_Ch2_med(:,:,1);
end

% Generate the general curvature
if Ch1_do
    gen_curve_Ch1 = squeeze(median(double(imdata_Ch1_med2), 1));
end
if Ch2_do
    gen_curve_Ch2 = squeeze(median(double(imdata_Ch2_med2), 1));
end

% Generate medians of the two channels to flatten them
if Ch1_do
    std_Ch1 = std(double(imdata_Ch1_med2), [], 3);
end
if Ch2_do
    std_Ch2 = std(double(imdata_Ch2_med2), [], 3);
end

% Bin the two stacks for calcultions
if Ch1_do
    bin_Ch1 = uint16(binxy(double(imdata_Ch1_med2), p.binxy));
end
if Ch2_do
    bin_Ch2 = uint16(binxy(double(imdata_Ch2_med2), p.binxy));
end

% Get an image of all ROIs flattened
if p.DoNeuropil
    ROIs_all = sum(ROIdata_corr.ROIs,3);
    ROIs_all = binxy(ROIs_all, p.binxy) > 0;
end

% Initialize signal structure
ROIsignal = struct('Section', [], 'ROI_index', [], 'Region', [],...
    'Size',[], 'bbox',[], 'Signal_Ch1', [], 'Signal_Ch2', [],...
    'image_Ch1', [], 'image_Ch2', [], 'Np_Signal_Ch1', [], 'Np_Signal_Ch2', []);
ROIsignal = repmat(ROIsignal,[ROIdata_corr.ncells, 1]);

% Keep track which ROIs are empty
FilledROI_vec = zeros(ROIdata_corr.ncells, 1);

% Loop through to apply filters
overall_ind = 0;

if p.ProgressBar
    hwait = waitbar(0, 'Processing');
end

for j = 1 : size(ROIdata_corr.ROIs,3)
    % Grab bounding box
    curr_bbox = regionprops(ROIdata_corr.ROIs(:,:,j),'BoundingBox');
    
    if p.ProgressBar
        hwait = waitbar(j/size(ROIdata_corr.ROIs,3), hwait,...
            sprintf('Processing %i/%i', j, size(ROIdata_corr.ROIs,3)));
    end
    
    for i = 1 : ROIdata_corr.ncells_vec(j)
       % Progress the index
       overall_ind = overall_ind + 1;
       
       % Fill in the basic ROI info
       ROIsignal(overall_ind).Section = j;
       ROIsignal(overall_ind).ROI_index = i;
       
       % Get the ROI
       Current_ROI = ROIdata_corr.ROIs(:,:,j) == i;
       
       % Bin ROI
       Current_ROI_bin = binxy(Current_ROI, p.binxy);
       
       % Get the ROI size
       ROIsignal(overall_ind).Size = sum(Current_ROI(:));
       
       % Get the ROI size (binned)
       size_binned = sum(Current_ROI_bin(:));
       
       % Determine if ROI is filled
       if ROIsignal(overall_ind).Size > 0
           % Flag the ROI filled
           FilledROI_vec(overall_ind) = 1;
           
           % Find out which region it belongs to
           if p.ProcessRegions
               ROIsignal(overall_ind).Region = round(median(ROIdata_corr.regions(Current_ROI)));
           end
           
           % Neuropil
           if p.DoNeuropil
               % Connect dots of the current ROI
               ROI_con = imclose(Current_ROI_bin, [1 1; 1 1]) > 0;
               NP_size = p.NeuropilRingSize;
               
               % Make the first neuropil
               Current_NP = imdilate(ROI_con, strel('disk', NP_size));
               Current_NP(ROIs_all) = 0;
               
               % Repeat if pixel count is not enough
               while sum(Current_NP(:)) < p.MinNeuropilPixels
                   % Make a bigger ring
                   NP_size = NP_size + 1;
                   Current_NP = imdilate(ROI_con, strel('disk', NP_size));
                   Current_NP(ROIs_all) = 0;
               end
           end
           
           % Individual filter RGB stack
           if p.makefilterRGB
               % Grab bounding box
               ROIsignal(overall_ind).bbox = curr_bbox(i).BoundingBox;

               % Add the gaps
               rgb_bbox_x = round((curr_bbox(i).BoundingBox(2) - p.gapsize) : ...
                   (curr_bbox(i).BoundingBox(2)+curr_bbox(i).BoundingBox(4) + p.gapsize));
               rgb_bbox_y = round((curr_bbox(i).BoundingBox(1) - p.gapsize) : ...
                   (curr_bbox(i).BoundingBox(1)+curr_bbox(i).BoundingBox(3) + p.gapsize));

               % Make sure the bounding box does not get out of range
               rgb_bbox_x = rgb_bbox_x(rgb_bbox_x > 0);
               rgb_bbox_y = rgb_bbox_y(rgb_bbox_y > 0);
               rgb_bbox_x = rgb_bbox_x(rgb_bbox_x <= size(Current_ROI, 1));
               rgb_bbox_y = rgb_bbox_y(rgb_bbox_y <= size(Current_ROI, 2));

               % Generathe rgb stacks
               % Channel 1
               if Ch1_do
                   rgb_ch1 = zeros(length(rgb_bbox_x),length(rgb_bbox_y), 3);
                   rgb_ch1(:,:,2) = mat2gray(std_Ch1(rgb_bbox_x, rgb_bbox_y));
                   rgb_ch1(:,:,1) = edge(Current_ROI(rgb_bbox_x, rgb_bbox_y));
               end

               % Channel 2
               if Ch2_do && Ch1_do
                   % Just take the stack from ch1 to go faster
                   rgb_ch2 = rgb_ch1;
               elseif Ch2_do
                   % Make Ch2 stack from scratch
                   rgb_ch2 = zeros(length(rgb_bbox_x),length(rgb_bbox_y), 3);
                   rgb_ch1(:,:,1) = edge(Current_ROI(rgb_bbox_x, rgb_bbox_y));
                   rgb_ch2(:,:,2) = mat2gray(std_Ch2(rgb_bbox_x, rgb_bbox_y));
               end

               % Load the rgb stacks
               if Ch1_do
                   ROIsignal(overall_ind).image_Ch1 = rgb_ch1;
               end
               if Ch2_do
                   ROIsignal(overall_ind).image_Ch2 = rgb_ch2;
               end
           end
           
           % Compute a signal tensor
           if Ch1_do
               Sig_tensor_Ch1 = ...
                   bin_Ch1 .* uint16(repmat(Current_ROI_bin, [1 1 sbxstruct_1.otparam(3)]));
               if p.DoNeuropil
                   Np_tensor_Ch1 = ...
                       bin_Ch1 .* uint16(repmat(Current_NP, [1 1 sbxstruct_1.otparam(3)]));
               end
           end
           if Ch2_do
               Sig_tensor_Ch2 = ...
                   bin_Ch2 .* uint16(repmat(Current_ROI_bin, [1 1 sbxstruct_2.otparam(3)]));
               if p.DoNeuropil
                   Np_tensor_Ch2 = ...
                       bin_Ch2 .* uint16(repmat(Current_NP, [1 1 sbxstruct_1.otparam(3)]));
               end
           end
           
           % Compute the signals
           if Ch1_do
               ROIsignal(overall_ind).Signal_Ch1 = ...
                   squeeze(sum(sum(Sig_tensor_Ch1,1),2)) / size_binned;
               if p.DoNeuropil
                   ROIsignal(overall_ind).Np_Signal_Ch1 = ...
                       squeeze(sum(sum(Np_tensor_Ch1,1),2)) / sum(Current_NP(:));
               end
           end
           if Ch2_do
               ROIsignal(overall_ind).Signal_Ch2 = ...
                   squeeze(sum(sum(Sig_tensor_Ch2,1),2)) / size_binned;
               if p.DoNeuropil
                   ROIsignal(overall_ind).Np_Signal_Ch2 = ...
                       squeeze(sum(sum(Np_tensor_Ch2,1),2)) / sum(Current_NP(:));
               end
           end
       end
    end
end

% Close
close(hwait)

% Remove the empty ROIs
ROIsignal = ROIsignal(FilledROI_vec == 1);

% Generate outputs for stds
if Ch1_do
    GenInfo.Ch1_std = std_Ch1;
    GenInfo.Ch1_curve = gen_curve_Ch1;
end
if Ch2_do
    GenInfo.Ch2_std = std_Ch2;
    GenInfo.Ch2_curve = gen_curve_Ch2;
end

end