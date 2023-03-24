% Create folder namelist of cell blocks
ori = pwd;
folders = dir(ori);
cellBlockList = [];
for folder = 1:size(folders,1)
    cellFlag = [folders(folder).isdir] && ~isempty(regexp(folders(folder).name, 'cellBlock[0-9]'));
    if cellFlag % Check if this folder is a cell
        % cells_list = [cells_list, str2double(folders(folder).name(5:end))];
        cellBlockList = [cellBlockList; convertCharsToStrings(folders(folder).name)];
    end
end

writematrix(cellBlockList, 'cellBlockList.txt');

