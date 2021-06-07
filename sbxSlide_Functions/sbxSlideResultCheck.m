function sbxSlideResultCheck(GenInfo, ROIsignal, ROIResults)
% sbxResultCheck shows the colocalization results visually.
% sbxSlideResultCheck(GenInfo, ROIsignal, ROIResults)

% Making RGB
rgb_stack = repmat(mat2gray(GenInfo.Ch2_std) + 0.1, [1 1 3]);
rgb_stack(:,:,2) = mat2gray(GenInfo.Ch1_std) + 0.1;
rgb_stack(:,:,3) = 0;
%%
% ROI to look
ROI2look = 1;

% Make figure
figure;

% Main figure
subplot(2, 7, [1:6, 8:13]);
imshow(rgb_stack);
hr = rectangle('Position', ROIsignal(ROI2look).bbox, 'EdgeColor','w');

% Make title
titlestring = ['ROI ', num2str(ROI2look),'/',num2str(size(ROIsignal,1))];
if ROIResults.Colocal(ROI2look)
    titlestring = [titlestring, ': POSITIVE'];
else
    titlestring = [titlestring, ': NEGATIVE'];
end
ht = title(titlestring);

% Green channel
subplot(2,7,7)
h1 = imshow(ROIsignal(ROI2look).image_Ch1, []);

% Red channel
subplot(2,7,14)
h2 = imshow(ROIsignal(ROI2look).image_Ch2, []);

ROI2look = input('Next ROI to look: ');
while ROI2look > 0
    % Update images
    set(hr, 'Position', ROIsignal(ROI2look).bbox);
    set(h1, 'CDATA', ROIsignal(ROI2look).image_Ch1);
    set(h2, 'CDATA', ROIsignal(ROI2look).image_Ch2);
    
    % Update title
    titlestring = ['ROI ', num2str(ROI2look),'/',num2str(size(ROIsignal,1))];
    if ROIResults.Colocal(ROI2look)
        titlestring = [titlestring, ': POSITIVE'];
    else
        titlestring = [titlestring, ': NEGATIVE'];
    end

    set(ht, 'String', titlestring);
    
    % Get new ROI to check
    ROI2look = input('Next ROI to look: ');
end

end