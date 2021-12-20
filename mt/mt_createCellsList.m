    % Control script is called from 'cells' directory;
    % cd into each 'cell' directory (e.g. cell1) to call mtcell object
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
    
%     for cell_idx = 1:length(cells_list)
%         
% %         cell = cells_list(cell_idx);
% %         cd(strcat(ori, '/cell', num2str(cell)));
%         
%         cell_no = cells_list(cell_idx);
%         cd(strcat(ori, '/', cell_no));
% 
%         fprintf('Calls mtcell for cell %d \n', cell); % Replace with code to call mtcell
%         
%         cd(ori);
%     end