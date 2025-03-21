% Control script is called from 'cells' directory;
% cd into each 'cell' directory (e.g. cell0001) to call mtcell object
ori = pwd;
cells_list = string(importdata('cells_list.txt'));

parfor cell_idx = 1:length(cells_list)
   
    fprintf('Cell %d \n', cell_idx);
    
%     if mod(cell_idx,20) == 0 || cell_idx == length(cells_list)
%         writematrix(cell_idx, 'processingCellNo.txt')
%     end
    
    cell_no = cells_list(cell_idx);
    cd(strcat(ori, '/', cell_no));
    
    try
        cellData = mtcell('auto');
        if cellData.data.Alpha ~= 50
        mtcell('Auto', 'Save', 'Redo');
        end
    catch
    end
    
    cd(ori);
end