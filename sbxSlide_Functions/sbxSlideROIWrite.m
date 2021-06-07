function ROIdata = sbxSlideROIWrite(ROIdata, icaguidata, movpath, varargin)
% sbxSlideROIWrite writes images of the ROIs
% ROIdata = sbxSlideROIWrite(ROIdata, icaguidata, movpath, varargin)

%% Parse inputs
p = inputParser;
addOptional(p, 'inputType', 'ACP'); % 'SigPass', 'Raw', 'ACP'

addOptional(p, 'UseGroupedSTDImage', true);
addOptional(p, 'NtoCombine', 1); % Combine ROIs from different frames or from the same STD Image

addOptional(p, 'UseROIEdge', true); % Show the edge of the ROI or not.

parse(p, varargin{:});
p = p.Results;

%% IO
% Number of ROI frames in total
switch p.inputType
    case 'SigPass'
        ROIs = ROIdata.SigPass.ROIs;
        RGB_stacks = ROIdata.rgbstack;
    case 'Raw'
        ROIs = ROIdata.SigPass.ROIs;
        RGB_stacks = ROIdata.SigPass.rgbstack;
    case 'ACP'
        ROIs = ROIdata.AreaCorrPass.ROIs;
        RGB_stacks = ROIdata.AreaCorrPass.rgb_stack;
end
nROI_frames = size(ROIs,3);

% Number of images to write (not counting the overall image)
if isfield(icaguidata, 'movstd_grouped')
    % Number of std frames
    Nstdim = size(icaguidata.movstd_grouped, 3);
else
    % Can't use grouped STD images since there are none
    p.UseGroupedSTDImage = false;
end

%% Group ROI frames
% Convert ROI data to 4D
if isfield(icaguidata, 'movstd_grouped')
    ROIs = reshape(ROIs, ...
        [size(ROIs,1), size(ROIs,2), Nstdim, nROI_frames/Nstdim]);
end

% Group as instructed first (4D)
if p.NtoCombine > 1 && isfield(icaguidata, 'movstd_grouped')
    % Get a grouping matrix
    gmat1 = 1 : p.NtoCombine : size(ROIs,4);
    gmat2 = p.NtoCombine: p.NtoCombine : size(ROIs,4);
    if length(gmat2) == length(gmat1) - 1
        gmat2 = [gmat2, size(ROIs,4)];
    end
    gmat = [gmat1;gmat2];
    gmat = gmat';
    
    % Initialize
    ROIs_g = zeros(size(ROIs,1), size(ROIs,2), size(ROIs,3), size(gmat,1));
    
    for i = 1 : size(gmat,1)
        % Indices
        i1 = gmat(i,1);
        i2 = gmat(i,2);
        
        % Combine frames
        ROIs_g(:,:,:,i) = sum(ROIs(:,:,:,i1:i2),4) > 0;
    end
    
    % Get the number of images (before splitting to different grouped stds)
    nIm = size(ROIs_g,4);
    
elseif isfield(icaguidata, 'movstd_grouped')
    % Just copy
    ROIs_g = ROIs;
    
    % Get the number of images (before splitting to different grouped stds)
    nIm = size(ROIs_g,4);
end

% Combine ROI frames if user electing to not using grouped std frames (4D)
if ~p.UseGroupedSTDImage && isfield(icaguidata, 'movstd_grouped')
    % Sum all the std frames together and squeeze
    ROIs_g = squeeze(sum(ROIs_g,3) > 0);
    
    % Get the number of images
    nIm = size(ROIs_g,3);
end

% Group as instructed first (3D)
if p.NtoCombine > 1 && ~isfield(icaguidata, 'movstd_grouped')
    % Get a grouping matrix
    gmat1 = 1 : p.NtoCombine : size(ROIs,3);
    gmat2 = p.NtoCombine: p.NtoCombine : size(ROIs,3);
    if length(gmat2) == length(gmat1) - 1
        gmat2 = [gmat2, size(ROIs,3)];
    end
    gmat = [gmat1;gmat2];
    gmat = gmat';
    
    % Initialize
    ROIs_g = zeros(size(ROIs,1), size(ROIs,2), size(gmat,1));
    
    for i = 1 : size(gmat,1)
        % Indices
        i1 = gmat(i,1);
        i2 = gmat(i,2);
        
        % Combine frames
        ROIs_g(:,:,i) = sum(ROIs(:,:,i1:i2),3) > 0;
    end
    
    % Get the number of images
    nIm = size(ROIs_g,3);
    
elseif ~isfield(icaguidata, 'movstd_grouped')
    % Just copy
    ROIs_g = ROIs;
    
    % Get the number of images
    nIm = size(ROIs_g,3);
end


%% Make RGB stacks
if p.UseGroupedSTDImage
    % Make a stack for each STD image
    RGB_stacks = repmat(RGB_stacks, [1 1 1 Nstdim]);
    
    for i = 1 : Nstdim
        % Replace the green channel
        RGB_stacks(:,:,2,i) = mat2gray(icaguidata.movstd_grouped(:,:,i));
        RGB_stacks(:,:,3,i) = 0;
    end
    
else
    % Clear the blue channel if not already
    RGB_stacks(:,:,3) = 0;
end

%% Paths
% Work out path
[fp, fn, ~] = fileparts(movpath);

% Output parth
fp_out = fullfile(fp, 'ROIs');
if ~exist(fp_out, 'dir')
    mkdir(fp_out);
end

% Output ROI filename
fn_ROI_all = sprintf('%s_Tot.tif',fn);
if p.UseGroupedSTDImage
    fn_ROI_gen =...
        sprintf('%s_%%0%id_%%0%id_ROI.tif',...
        fn, ceil(log(Nstdim)/log(10)), ceil(log(nIm)/log(10)));
    
    fn_ROI_gen_stdoverall = sprintf('%s_%%0%id_Sec.tif',...
        fn, ceil(log(Nstdim)/log(10)));
else
    fn_ROI_gen =...
        sprintf('%s_%%0%id_ROI.tif', fn, ceil(log(nIm)/log(10)));
end

%% Write
% Individual ROI frames (or combined)
% 4D or 3D
if p.UseGroupedSTDImage
    % 4D
    for j = 1 : Nstdim
        % Volume
        rgb = RGB_stacks(:,:,:,j);

        for i = 1 : nIm
            % Fill the red channel
            if p.UseROIEdge
                rgb(:,:,1) = edge(ROIs_g(:,:,j,i)) > 0;
            else
                rgb(:,:,1) = ROIs_g(:,:,j,i) > 0;
            end

            % Filename
            fn_ROI = sprintf(fn_ROI_gen, j, i);

            % Write
            imwrite(rgb, fullfile(fp_out, fn_ROI));
        end
    end
else
    % 3D
    % Copy over
    rgb = RGB_stacks;
    
    for i = 1 : nIm
        % Fill the red channel
        if p.UseROIEdge
            rgb(:,:,1) = edge(ROIs_g(:,:,i)) > 0;
        else
            rgb(:,:,1) = ROIs_g(:,:,i) > 0;
        end
        
        % Filename
        fn_ROI = sprintf(fn_ROI_gen, i);
        
        % Write
        imwrite(rgb, fullfile(fp_out, fn_ROI));
    end
end

% Overall per grouped std frame
% 4D or 3D
if p.UseGroupedSTDImage
    for j = 1 : Nstdim
        % Volume
        rgb = RGB_stacks(:,:,:,j);
        
        % Fill the red channel
        if p.UseROIEdge
            rgb(:,:,1) = edge(sum(ROIs_g(:,:,j,:),4)) > 0;
        else
            rgb(:,:,1) = sum(ROIs_g(:,:,j,:),4) > 0;
        end

        % Filename
        fn_ROI = sprintf(fn_ROI_gen_stdoverall, j);

        % Write
        imwrite(rgb, fullfile(fp_out, fn_ROI));
    end
end

% Overall
switch p.inputType
    case 'SigPass'
        imwrite(ROIdata.rgbstack, fullfile(fp_out, fn_ROI_all));
    case 'Raw'
        imwrite(ROIdata.SigPass.rgbstack, fullfile(fp_out, fn_ROI_all));
    case 'ACP'
        imwrite(ROIdata.AreaCorrPass.rgb_stack, fullfile(fp_out, fn_ROI_all));
end

% Output ROI
ROIdata.ROIwrite.gmat = gmat;
ROIdata.ROIwrite.ROIs_g = ROIs_g;
ROIdata.ROIwrite.par = p;
end