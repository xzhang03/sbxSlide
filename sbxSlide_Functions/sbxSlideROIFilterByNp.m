function ROI_data = sbxSlideROIFilterByNp(ROI_data,varargin)
% sbxSlideROIFilterByNP removes bad ROIs by signal and neuropil signal
% ROI_data = sbxSlideROIFilterByNp(ROI_data,varargin)
% By Stephen Zhang 3/7/2020

%% Parse Input
% Parse input
p = inputParser;

addOptional(p, 'SubtractNeuropil', true);
addOptional(p, 'Channel', 1);
addOptional(p, 'Method', 'ExtremaToMedian'); % Can be ExtremaToMedian, PeakToMedian, PeakToMean
addOptional(p, 'Threshold', 0); % Threshold value after applying the method

parse(p, varargin{:});
p = p.Results;

%% Initialize
% N cells
ncells = ROI_data.ncells;

% A vector to record data
datavec = zeros(ncells, 1);


% Fields
SigField = sprintf('Signal_Ch%i', p.Channel);
NpSigField = sprintf('Np_Signal_Ch%i', p.Channel);

%% Loop to apply method
for ROI_ind = 1 : ncells
    % Grab signal
    if p.SubtractNeuropil
        CurrSignal = ROI_data.signals(ROI_ind).(SigField) - ...
            ROI_data.signals(ROI_ind).(NpSigField);
    else
        CurrSignal = ROI_data.signals(ROI_ind).(SigField);
    end
    
    % Apply method
    switch p.Method
        case 'PeakToMedian'
            % Max - Median
            datavec(ROI_ind) = max(CurrSignal) - median(CurrSignal);
        case 'PeakToMean'
            % Max - Mean
            datavec(ROI_ind) = max(CurrSignal) - mean(CurrSignal);
        case 'ExtremaToMedian'
            % [Max - Median] or [Min - Median], whichever one with a
            % greater numerical value
            pos = max(CurrSignal) - median(CurrSignal);
            neg = min(CurrSignal) - median(CurrSignal);
            
            if pos > -neg
                datavec(ROI_ind) = pos;
            else
                datavec(ROI_ind) = neg;
            end
        case 'LocalExtremaToMedian'
            % Find local extrema
            [localmax, Pmax] = islocalmax(CurrSignal, 'MaxNumExtrema', 1);
            [localmin, Pmin] = islocalmin(CurrSignal, 'MaxNumExtrema', 1);
            
            pos = CurrSignal(localmax) - median(CurrSignal);
            neg = CurrSignal(localmin) - median(CurrSignal);
            
            if pos > -neg && pos > 0 % Peak above median and Peak amplitude is greater than valley depth
                datavec(ROI_ind) = Pmax(localmax);
            else
                datavec(ROI_ind) = -Pmin(localmin);
            end
    end
end

%% Pass and fail
% A vector to note pass or fail
SigPass = datavec > p.Threshold;

%% Clean up ROIs
% Initialize
ROIs_SigPass = ROI_data.ROIs;
ncells_vec_SigPass = ROI_data.ncells_vec;
cell_counter = 0;

% Loop through
for j = 1 : length(ncells_vec_SigPass)
    % Pass counter
    Curr_pass_counter = 0;
    
    % Current frame
    Curr_frame = ROIs_SigPass(:,:,j);
    
    for i = 1 : ncells_vec_SigPass(j)
        % Propagate counter
        cell_counter = cell_counter + 1;
        
        % If Pass
        if SigPass(cell_counter)
            Curr_pass_counter = Curr_pass_counter + 1;
        else % Fail
            Curr_frame(Curr_frame == i) = 0;
        end
    end
    
    % Update frame
    ROIs_SigPass(:,:,j) = Curr_frame;
    
    % Update cell counter
    ncells_vec_SigPass(j) = Curr_pass_counter;
end

% Update RGB stack
rgb_stack_SigPass = ROI_data.rgbstack;
rgb_stack_SigPass(:,:,1) = edge(sum(ROIs_SigPass,3) > 0);

%% Save 
ROI_data.SigPass.rgbstack = rgb_stack_SigPass;
ROI_data.SigPass.ROIs = ROIs_SigPass;
ROI_data.SigPass.ncells_vec = ncells_vec_SigPass;
ROI_data.SigPass.Pass = SigPass;
ROI_data.SigPass.ncell = sum(ncells_vec_SigPass);
ROI_data.SigPass.par = p;
end