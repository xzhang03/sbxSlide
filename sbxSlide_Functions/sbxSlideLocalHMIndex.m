function [inds, HMwidth] = sbxSlideLocalHMIndex(signal, varargin)
% sbxLocalHMIndex finds the half max indices of one peak. Each column of
% the signal matrix is its own sweep
% inds, HMwidth] = sbxLocalHMIndex(signal, varargin)
% By Stephen Zhang, 2019/05/20

% Mode
p = inputParser;
addOptional(p, 'HMMode', 'FirstPointBelowHM'); 
parse(p, varargin{:});
p = p.Results;


% Find peaks (one per column)
[max_vals, max_inds] = max(signal, [], 1);

% Find cases where the max values are on the extreme left or right hand
% side
left_peak_sweeps = max_inds == 1;
right_peak_sweeps = max_inds == size(signal, 1);

% Initialize index output
inds = nan(size(signal,2),2);

% Determine how the indices are calculated: last point above the HM
% (starting from peak) or the first point that dips below the HM (startin
% from the peak).

switch p.HMMode
    case 'FirstPointBelowHM'
        % Loop through
        for i = 1 : size(signal, 2)
            if ~left_peak_sweeps(i) && ~right_peak_sweeps(i)
                % Find where the minimals are in both directions
                [left_min_val, ~] = min(signal(1 : max_inds(i), i), [], 1);
                [right_min_val, ~] = min(signal(max_inds(i) : end, i), [], 1);

                % Find out where the half-max values are
                HMleft = left_min_val + (max_vals(i) - left_min_val) / 2;
                HMright = right_min_val + (max_vals(i) - right_min_val) / 2;
        
                % Find how when the signals cross the half-max boundaries
                inds(i,1) = find(signal(1 : max_inds(i), i) <= HMleft, 1, 'last');
                inds(i,2) = find(signal(max_inds(i) : end, i) <= HMright, 1, 'first') + max_inds(i) - 1;
            end
        end
    case 'LastPointAboveHM'
        % Loop through
        for i = 1 : size(signal, 2)
            if ~left_peak_sweeps(i) && ~right_peak_sweeps(i)
                % Find where the minimals are in both directions
                [left_min_val, ~] = min(signal(1 : max_inds(i), i), [], 1);
                [right_min_val, ~] = min(signal(max_inds(i) : end, i), [], 1);

                % Find out where the half-max values are
                HMleft = left_min_val + (max_vals(i) - left_min_val) / 2;
                HMright = right_min_val + (max_vals(i) - right_min_val) / 2;
                
                % Find how when the signals cross the half-max boundaries
                inds(i,1) = find(signal(1 : max_inds(i), i) >= HMleft, 1, 'first');
                inds(i,2) = find(signal(max_inds(i) : end, i) >= HMright, 1, 'last') + max_inds(i) - 1;
            end
        end
end

% Calculate width
HMwidth = inds(:,2) - inds(:,1) + 1;
end