function outputtable = sbxSlideResultTable(ROIResults, ROIsignal)
% sbxSlideResultTable summarizes results into a table
% outputtable = sbxSlideResultTable(ROIResults, ROIsignal)
% Stephen Zhang 2019/05/21

% Region assignments
Region_vec = [ROIsignal.Region];
Signal_vec = ROIResults.Colocal;

% Number of regions
n_regions = max(Region_vec);

% Initialize output
outputtable = nan(n_regions + 1, 4);

% Loop through
for i = 0 : n_regions
    % Write region info
    outputtable(i + 1, 1) = i;
    
    % Write the number of cells
    currentregion_ind = Region_vec == i;
    outputtable(i + 1, 2) = sum(currentregion_ind);
    outputtable(i + 1, 3) = sum(Signal_vec(currentregion_ind));
    
    % Percent of colocalization
    outputtable(i + 1, 4) = outputtable(i + 1, 3)/outputtable(i + 1, 2);
end

% make table
outputtable = array2table(outputtable, 'VariableNames', {'Region', 'N_Cells', 'N_Colocalizations', 'Fraction_Colocalizations'});
end