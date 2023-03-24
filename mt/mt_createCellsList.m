% Create list of cells
ori = pwd;
folders = dir(ori);
cells_list = [];
for folder = 1:size(folders,1)
    cellFlag = [folders(folder).isdir] && ~isempty(regexp(folders(folder).name, 'cell[0-9]'));
    if cellFlag % Check if this folder is a cell
        % cells_list = [cells_list, str2double(folders(folder).name(5:end))];
        cells_list = [cells_list; convertCharsToStrings(folders(folder).name)];
    end
end

writematrix(cells_list, 'cells_list.txt');
