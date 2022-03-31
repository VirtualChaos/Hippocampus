function [obj, varargout] = mtneuraldata(varargin)
%@mtneuraldata Constructor function for mtneuraldata class
%   OBJ = mtneuraldata(varargin)
%
%   OBJ = mtneuraldata('auto') attempts to create a mtneuraldata object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on mtneuraldata %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%example [as, Args] = mtneuraldata('save','redo')
%
%dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
				'ObjectLevel','Session','RequiredFile','cells_list.txt', 'NumericArguments', [], ...
				'BinSize',100, 'ThresVel',0);
            
Args.flags = {'Auto','ArgsOnly'};
% Specify which arguments should be checked when comparing saved objects
% to objects that are being asked for. Only arguments that affect the data
% saved in objects should be listed here.
Args.DataCheckArgs = {'BinSize', 'ThresVel'};                         

[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'mtneuraldata';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'mtneuraldata';

% To decide the method to create or load the object
[command,robj] = checkObjCreate('ArgsC',Args,'narginC',nargin,'firstVarargin',varargin);

if(strcmp(command,'createEmptyObjArgs'))
    varargout{1} = {'Args',Args};
    obj = createEmptyObject(Args);
elseif(strcmp(command,'createEmptyObj'))
    obj = createEmptyObject(Args);
elseif(strcmp(command,'passedObj'))
    obj = varargin{1};
elseif(strcmp(command,'loadObj'))
    % l = load(Args.matname);
    % obj = eval(['l.' Args.matvarname]);
	obj = robj;
elseif(strcmp(command,'createObj'))
    % IMPORTANT NOTICE!!! 
    % If there is additional requirements for creating the object, add
    % whatever needed here
    obj = createObject(Args,modvarargin{:});
end

function obj = createObject(Args,varargin)

% check if the right conditions were met to create object
if(~isempty(dir(Args.RequiredFile)))

    ori = pwd;
    data.origin = {pwd};
    foldername = pwd;
    
    cd ..
    sessionData = mtsess('Auto', varargin{:});
    sessionData = sessionData.data;
    cd(ori);
    
    cells_list = string(importdata(Args.RequiredFile));
    
    % placecellStats_combined = zeros(length(cells_list), 4, 10);
    % for cycle = 1:10
    isplacecell = zeros(length(cells_list),5);
    shuffled_sic_values = zeros(length(cells_list), 10000);
    mtplacefield_failed = [];
    
    tic
    for cell_idx = 1:length(cells_list)
        
        cell_no = cells_list(cell_idx);
        cd(strcat(ori, '/', cell_no));
        
        if mod(cell_idx,50) == 0
            fprintf("Current Cell: %d / %d; Estimated Time Remaining: %.2f seconds\n", cell_idx, length(cells_list), (toc/cell_idx) * length(cells_list) - toc);
        end
        
        try
            cellData = load('mtcell.mat');
            cellData = cellData.mtcell.data;
        catch
            cellData = mtcell('auto','save','redo').data;
        end
        
        try
            mtplacefield('auto', 'save', 'redo');
        catch
            mtplacefield_failed = cat(1,mtplacefield_failed, cell_idx);
        end
        
        cell_no_ = char(cell_no);
        cell_no_ = str2double(cell_no_(5:end));
        isplacecell(cell_idx,1) = cell_no_;
        isplacecell(cell_idx,2) = cellData.SIC;
        isplacecell(cell_idx,3) = prctile(cellData.SICsh,95);
        isplacecell(cell_idx,4) = cellData.SIC >= prctile(cellData.SICsh,95);
        shuffled_sic_values(cell_idx,:) = cellData.SICsh;
        
        cellData = rmfield(cellData, {'maps_adsmsh', 'dur_adsmsh', 'radiish'});
        data.cellData.("n" + sprintf('%04d',cell_no_)) = cellData;
        
        data.mtplacefield_failed = mtplacefield_failed;
        
        try
            placefieldData = load('mtplacefield.mat');
            placefieldData = placefieldData.mtplacefield.data;
            spiketimes = load('spiketimes.mat');
            spiketimes = spiketimes.spiketimes;
            spiketrain = load('spiketrain.mat');
            spiketrain = spiketrain.spiketrain;
        catch
            placefieldData = {};
        end
        data.placefieldData.("n" + sprintf('%04d',cell_no_)) = placefieldData;
        data.spiketimes.("n" + sprintf('%04d',cell_no_)) = spiketimes;
        data.spiketrain.("n" + sprintf('%04d',cell_no_)) = spiketrain;
        
        cd(ori);
    end
    isplacecell(:,5) = isplacecell(:,2) > prctile(shuffled_sic_values,95,'all');
    
    data.isplacecell = sortrows(isplacecell);
    data.shuffled_sic_values = shuffled_sic_values;
    data.cellData = orderfields(data.cellData);
    data.placefieldData = orderfields(data.placefieldData);
    
    rmse = zeros(length(cells_list), 2);
    nrmse_mean = zeros(length(cells_list), 2);
    nrmse_stdev = zeros(length(cells_list), 2);
    auc = zeros(length(cells_list), 2);
    placecellStats = zeros(length(cells_list), 4);
    placefieldStats = [];
    fns = fieldnames(data.placefieldData);
    for i = 1:length(cells_list)
        rmse(i,1) = i;
        nrmse_mean(i,1) = i;
        nrmse_stdev(i,1) = i;
        auc(i,1) = i;
        placecellStats(i,1) = i;
        if ~isempty(data.placefieldData.(fns{i})) && logical(data.isplacecell(i,5))
            rmse(i,2) = data.placefieldData.(fns{i}).rmse;
            nrmse_mean(i,2) = data.placefieldData.(fns{i}).nrmse_mean;
            nrmse_stdev(i,2) = data.placefieldData.(fns{i}).nrmse_stdev;
            auc(i,2) = data.placefieldData.(fns{i}).auc;
            placecellStats(i,2) = size(data.placefieldData.(fns{i}).GMM,1); % No. of place fields
            placecellStats(i,3) = mean(data.cellData.(fns{i}).maps_adsm,'omitnan'); % Mean firing rate
            placecellStats(i,4) = std(data.cellData.(fns{i}).maps_adsm,'omitnan'); % Stdev firing rate
            if ~isempty(data.placefieldData.(fns{i}).GMM)
                cell_no_ = char(cells_list(i));
                cell_no_ = str2double(cell_no_(5:end));
                temp = repmat(cell_no_,size(data.placefieldData.(fns{i}).GMM,1),1);
                placefieldStats = cat(1, placefieldStats, [temp data.placefieldData.(fns{i}).GMM]);
            end
        else
            rmse(i,2) = NaN;
            nrmse_mean(i,2) = NaN;
            nrmse_stdev(i,2) = NaN;
            auc(i,2) = NaN;
            placecellStats(i,2) = NaN;
            placecellStats(i,3) = NaN;
            placecellStats(i,4) = NaN;
        end
    end
    
%     rmse_sort = sortrows(rmse,2);
%     nrmse_mean_sort = sortrows(nrmse_mean,2);
%     nrmse_stdev_sort = sortrows(nrmse_stdev,2);
%     auc_sort = sortrows(auc,2);
    
    data.placecellStats = placecellStats;
    data.placefieldStats = placefieldStats;
    data.nTrials = sessionData.nTrials;
    data.nNeuron = sessionData.nNeuron;
    
    % create nptdata so we can inherit from it
    data.numSets = 0;
    data.Args = Args;
    n = nptdata(1,0,pwd);
    d.data = data;
    obj = class(d,Args.classname,n);
    saveObject(obj,'ArgsC',Args);

else
	% create empty object
	obj = createEmptyObject(Args);
end

function obj = createEmptyObject(Args)

% these are object specific fields
data.dlist = [];
data.setIndex = [];

% create nptdata so we can inherit from it
% useful fields for most objects
data.numSets = 0;
data.Args = Args;
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);
