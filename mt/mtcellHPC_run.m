% Control script is called from 'cells' directory;
% cd into each 'cell' directory (e.g. cell0001) to call mtcell object
ori = pwd;

for cell_idx = 1:length(cells_list)
    
    cell_no = cells_list(cell_idx);
    cd(strcat(ori, '/', cell_no));
    
    if mod(cell_idx,10) == 0
        toc
    end
    
    try
        mtcell('Auto', 'Save', 'Redo');
    catch
    end
    
    cd(ori);
end