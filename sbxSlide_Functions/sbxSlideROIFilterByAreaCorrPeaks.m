function ROIdata = sbxSlideROIFilterByAreaCorrPeaks(ROIdata, varargin)
% sbxSlideROIFilterByAreaPeaks filters out ROIs based on 1) area overlap, 
% 2) corrlation, as well as 3) the number of peaks under area. An ROI is
% thrown out if both criteria are met.
% ROIdata = sbxSlideROIFilterByAreaCorr(ROIdata, varargin)
% By Stephen Zhang 3/9/2020

p = inputParser;
addOptional(p, 'inputType', 'SigPass'); % 'SigPass', 'Raw'

% Area parameters
addOptional(p, 'AreaOverlap', 0.9); % Overlap area is at least X% of either ROI

% Correlation parameters
addOptional(p, 'Corr', 0.9); % Correlation threshold


% Peak parameters
addOptional(p, 'SmoothWindow', 5); % Window size for smoothing, 0 to turn off
addOptional(p, 'PeakProminenceThreshold', 500); % Threshold for peak prominence.

parse(p, varargin{:});
p = p.Results;

%% IO
switch p.inputType
    case 'Raw'
        signals = ROIdata.signals;
        ROIs = ROIdata.ROIs;
        ncells_vec = ROIdata.ncells_vec;
        rgb_stack = ROIdata.rgbstack;
    case 'SigPass'
        signals = ROIdata.signals(ROIdata.SigPass.Pass);
        ROIs = ROIdata.SigPass.ROIs;
        ncells_vec = ROIdata.SigPass.ncells_vec;
        rgb_stack = ROIdata.SigPass.rgbstack;
end

%% Check area overlaps
%
% Number of comparisons (as an upper limit)
ncomps = sum(cumsum(ncells_vec(2 : end)) .* ncells_vec(1 : end -1));

% Data mat to store calculation results 
% [Frame_1 Cell_1 Frame_2 Cell_2 AreaOv Corr nPeaks]
Data_mat = zeros(ncomps, 7);

% counter
counter = 0;

% Waitbar
hwait = waitbar(0, 'Processing');

% Loop through
% Frame 1
for F1 = 1 : length(ncells_vec) - 1
    
    % Update waitbar
    waitbar(F1/(length(ncells_vec) - 1), hwait, sprintf('Processing %i/%i', F1, length(ncells_vec) - 1));
    
    % Grab frame
    Frame_1 = ROIs(:,:,F1);
    
    % Frame 2
    for F2 = (F1+1) : length(ncells_vec)
        % Grab frame
        Frame_2 = ROIs(:,:,F2);

        % Grab overlap between frames
        Ov_F = (Frame_1 + Frame_2) > 1;
        
        if sum(Ov_F(:)) > 0
            % If there is overlap at all
            % Find out potential participating cells
            % F1
            F1_cellvec = Frame_1(Ov_F);
            F1_cellvec = unique(F1_cellvec);
            F1_cellvec = F1_cellvec(F1_cellvec >0);

            % F2
            F2_cellvec = Frame_2(Ov_F);
            F2_cellvec = unique(F2_cellvec);
            F2_cellvec = F2_cellvec(F2_cellvec >0);
            
            % Loop through cells in F1
            for C1 = 1 : length(F1_cellvec)
                % Grab Cell image and size
                CellId_1 = F1_cellvec(C1);
                CellIm_1 = Frame_1 == CellId_1;
                CellSz_1 = sum(CellIm_1(:));

                for C2 = 1 : length(F2_cellvec)
                    % Grab Cell image and size
                    CellId_2 = F2_cellvec(C2);
                    CellIm_2 = Frame_2 == CellId_2;
                    CellSz_2 = sum(CellIm_2(:));

                    % Grab cell overlaps
                    Ov_c = CellIm_1 & CellIm_2;
                    Ovsz = sum(Ov_c(:));
                    
                    if Ovsz > 0
                        % Propagate counter
                        counter = counter + 1;
                    
                        % Filling
                        Data_mat(counter,1) = F1;
                        Data_mat(counter,2) = CellId_1;
                        Data_mat(counter,3) = F2;
                        Data_mat(counter,4) = CellId_2;
                        Data_mat(counter,5) = max(Ovsz/CellSz_1, Ovsz/CellSz_2);
                    end
                end
            end
        end
    end
end

close(hwait)

% Downsize mat
Data_mat = Data_mat(1:counter, :);


%}
%% Check the number of peaks
% Comparisons to check
% Data_mat = ROIdata.AreaCorrPass.Data_mat; % Debugging code
compstocheck = find(Data_mat(:,5) >= p.AreaOverlap);

% Vectors for locating signals
Sec_vec = [signals(:).Section];
ROI_ind_vec = [signals(:).ROI_index];

% Loop through
for i = 1 : length(compstocheck)
    % Index
    compind = compstocheck(i);
    
    % Signal indices
    Sig1_ind = ...
        (Sec_vec == Data_mat(compind,1)) & (ROI_ind_vec == Data_mat(compind,2));
    Sig2_ind = ...
        (Sec_vec == Data_mat(compind,3)) & (ROI_ind_vec == Data_mat(compind,4));
    
    % Signals
    Sig1 = signals(Sig1_ind).Signal_Ch1;
    Sig2 = signals(Sig2_ind).Signal_Ch1;
    
    % Correlation
    SigCorr = corr(Sig1, Sig2);
    Data_mat(compind,6) = SigCorr;
    
    if SigCorr >= p.Corr
        % Smooth
        if p.SmoothWindow > 0
            Sig1 = smooth(Sig1, p.SmoothWindow);
            Sig2 = smooth(Sig2, p.SmoothWindow);
        end

        % Prominence of Sig1
        [~, P1] = islocalmax(Sig1);
        [~, P2] = islocalmax(Sig2);

        % Save the number of prominent peaks
        Data_mat(compind,7) =...
            max(sum(P1 >= p.PeakProminenceThreshold), sum(P2 >= p.PeakProminenceThreshold));
    end
    
%     Debug
%     plot(1:length(Sig1), Sig1, 1:length(Sig2), Sig2);
%     title(num2str(Data_mat(compind,7)));
end

% Grab pass and fail
PassOrFail = (Data_mat(:,6) < p.Corr) | (Data_mat(:,7) > 1);


%% Clean up the images
% A matrix to label which cells are to be removed
Rem_mat = Data_mat(~PassOrFail, 3:4);

for i = 1 : size(Rem_mat,1)
    % Frame and cell
    F = Rem_mat(i,1);
    C = Rem_mat(i,2);
    
    % Current frame
    CurrFrame = ROIs(:,:,F);
    
    % Remove cell
    CurrFrame(CurrFrame == C) = 0;
    ncells_vec(F) = ncells_vec(F) - 1;
    ROIs(:,:,F) = CurrFrame;
end

%% Output
% Numbers
ROIdata.AreaCorrPass.PassOrFail = PassOrFail;
ROIdata.AreaCorrPass.Data_mat = Data_mat;
ROIdata.AreaCorrPass.ROIs = ROIs;
ROIdata.AreaCorrPass.ncells_vec = ncells_vec;
ROIdata.AreaCorrPass.par = p;
ROIdata.AreaCorrPass.ncells = sum(ncells_vec);

% RGB stack
rgb_stack(:,:,1) = edge(sum(ROIs,3) > 0);
ROIdata.AreaCorrPass.rgb_stack = rgb_stack;
end