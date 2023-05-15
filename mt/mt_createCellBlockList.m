% Create folder namelist of cell blocks
ori2 = pwd;
folders = dir(ori2);
cellBlockList = [];
for folder = 1:size(folders,1)
    cellFlag = [folders(folder).isdir] && ~isempty(regexp(folders(folder).name, 'cellBlock[0-9]'));
    if cellFlag % Check if this folder is a cell
        % cells_list = [cells_list, str2double(folders(folder).name(5:end))];
        cellBlockList = [cellBlockList; convertCharsToStrings(folders(folder).name)];
    end
end

writematrix(cellBlockList, 'cellBlockList.txt');

cd ..
[~,curr_folder,~] = fileparts(pwd);
cd(ori2);
flag = regexp(curr_folder,'rest','match');
if string(regexp(curr_folder,'rest','match')) == "rest"
    writematrix('rest', 'vel_threshold_handler.txt');
elseif string(regexp(curr_folder,'zero','match')) == "zero"
    writematrix('zero', 'vel_threshold_handler.txt');
else
    writematrix('normal', 'vel_threshold_handler.txt');
end


