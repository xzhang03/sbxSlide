function [rgb_colocal, rgb_ortho] = sbxSlideMakeRGBfromROI(ROIdata_corr, ROIsignal, GenInfo, ROIResults)
% sbxSlideMakeRGBfromROI generates two RGB stacks, one for colocalized
% cells and the other for orthogonal cells.
%[rgb_colocal, rgb_ortho] = sbxSlideMakeRGBfromROI(ROIdata_corr, ROIsignal, im_stds, ROIResults)
% Stephen Zhang 2019/05/21

% Initialize outputs
rgb_colocal = repmat(GenInfo.Ch2_std, [1 1 3]);
rgb_colocal(:,:,2) = 0;
rgb_colocal(:,:,3) = 0;

rgb_ortho = rgb_colocal;

% Loop through to parse the ROIs
for i = 1 : length(ROIResults.Colocal)
    % Grab ROI info
    Section_ind = ROIsignal(i).Section;
    ROI_ind = ROIsignal(i).ROI_index;
    
    % Grab current ROI
    ROI = ROIdata_corr.ROIs(:,:,Section_ind) == ROI_ind;
    
    % Parse based on if the ROI is colocalized or not
    if ROIResults.Colocal(i)
        rgb_colocal(:,:,3) = rgb_colocal(:,:,3) + ROI;
    else
        rgb_ortho(:,:,3) = rgb_ortho(:,:,3) + ROI;
    end
end

% Add the regions
rgb_colocal(:,:,3) = edge(rgb_colocal(:,:,3)) + edge(ROIdata_corr.regions);
rgb_ortho(:,:,3) = (edge(rgb_ortho(:,:,3)) + edge(ROIdata_corr.regions)) > 0;

% Scale each channel
for i = 1 : 3
    rgb_colocal(:,:,i) = mat2gray(rgb_colocal(:,:,i)) + 0.2;
    rgb_ortho(:,:,i) = mat2gray(rgb_ortho(:,:,i)) + 0.2;
end



end