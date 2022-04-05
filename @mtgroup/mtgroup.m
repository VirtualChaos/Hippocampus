function [obj, varargout] = mtgroup(varargin)
%@mtgroup Constructor function for mtgroup class
%   OBJ = mtgroup(varargin)
%
%   OBJ = mtgroup('auto') attempts to create a mtgroup object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on mtgroup %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%example [as, Args] = mtgroup('save','redo')
%
%dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
				'ObjectLevel','Mice','RequiredFile','', 'NumericArguments', [], ...
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
Args.classname = 'mtgroup';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'mtgroup';

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
if(true) % ~isempty(dir(Args.RequiredFile))

    ori = pwd;
    data.origin = {pwd};
    foldername = pwd;
    
    files = dir;
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    
    for i = 1:size(subFolders,1)
        if  subFolders(i).name(1) == 'I' && subFolders(i).name(2) == 'D'
            cd(subFolders(i).name);
            disp(subFolders(i).name)
            
            if isfile('mtmice.mat')
                temp = mtmice('auto').data;
                data.miceData.(subFolders(i).name) = temp.neuralCombined;
                
                fns = fieldnames(temp.sessionCombined);
                for j = 1:length(fns)
                    field = ["F", "Fc", "spikes_corrected", "dF_F0_corrected", "session_data_raw", "acceleration_averaged", "session_data_exclude_zero_trials"];
                    temp.sessionCombined.(fns{j}) = rmfield(temp.sessionCombined.(fns{j}), field);
                end
                data.miceData_sessionData.(subFolders(i).name) = temp.sessionCombined;
            end
                
            cd(ori)
        end
    end

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
