function ROIResults = sbxSlideColocal(ROIsignal, GenInfo, varargin)
% sbxSlideColocal determines whether the two chanels are colocalized or not
% ROIResults = sbxSlideColocal(ROIsignal, GenInfo, varargin)
% Stephen Zhang 2019/05/21

% Tolerance
p = inputParser;
addOptional(p, 'Tolerance', 1, @isnumeric); 
addOptional(p, 'HMMode', 'FirstPointBelowHM'); 
parse(p, varargin{:});
p = p.Results;

% Extract the data into matrices
Signal_Ch1_mat = cell2mat({ROIsignal.Signal_Ch1});
Signal_Ch2_mat = cell2mat({ROIsignal.Signal_Ch2});

% Generate background curves
Background_Ch1_mat = nan(size(Signal_Ch1_mat));
Background_Ch2_mat = nan(size(Signal_Ch2_mat));

% Loop through the ROIs
for i = 1 : size(ROIsignal,1)
    % Figure out the start and end indices
    ind_left = round(ROIsignal(i).bbox(1));
    ind_right = ind_left + ROIsignal(i).bbox(3) - 1;
    
    % Average the relevant segments of the background curve to get a
    % general curve for that ROI.
    Background_Ch1_mat(:,i) = mean(GenInfo.Ch1_curve(ind_left:ind_right, :), 1);
    Background_Ch2_mat(:,i) = mean(GenInfo.Ch2_curve(ind_left:ind_right, :), 1);
end

% Find the half-max indices
HMinds_Ch1 = sbxSlideLocalHMIndex(Signal_Ch1_mat - Background_Ch1_mat, 'HMMode', p.HMMode);
[HMinds_Ch2, ~] = sbxSlideLocalHMIndex(Signal_Ch2_mat - Background_Ch2_mat, 'HMMode', p.HMMode);

% Check left hand
Left_check = (HMinds_Ch1(:,1) - p.Tolerance) <= HMinds_Ch2(:,1);

% Check right hand
Right_check = (HMinds_Ch1(:,2) + p.Tolerance) >= HMinds_Ch2(:,2);

% Output
ROIResults.Colocal = Left_check & Right_check;
ROIResults.Tolerance = p.Tolerance;
end