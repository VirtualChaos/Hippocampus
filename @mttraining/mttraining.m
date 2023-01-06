function [obj, varargout] = mtgroup(varargin)
%@mtgroup Constructor function for mttraining class
%   OBJ = mttraining(varargin)
%
%   OBJ = mttraining('auto') attempts to create a mttraining object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on mttraining %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%example [as, Args] = mttraining('save','redo')
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
Args.classname = 'mttraining';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'mttraining';

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
    
    files = dir;
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    
    nMice = 0;
    for i = 1:size(subFolders,1)
        if  subFolders(i).name(1) == 'I' && subFolders(i).name(2) == 'D'
            cd(subFolders(i).name);
            disp(subFolders(i).name)
            
            mice_no = str2num(subFolders(i).name(3:4));
            
            if isfile('Stats/training_stats.mat')
            
                if ismember(mice_no,[45 49 50 59 60 61])
                    group_id = 'Y_Ctrl';
                    
                elseif ismember(mice_no,[33 40 41 47 48 58 81 83 84 85])
                    group_id = 'Y_APP';
                    
                elseif ismember(mice_no,[38 39 55 56 64 65 68 69 73 74 77])
                    group_id = 'O_Ctrl';
                    
                elseif ismember(mice_no,[35 37 53 54 62 63 67 70 72 80])
                    group_id = 'O_APP';
                end
                
                try
                    data.sessionCombined.(group_id).(subFolders(i).name(1:4)) = load('Stats/training_stats.mat').training_stats;
                catch
                    miceData = mtmice('auto','Training');
                    data.sessionCombined.(group_id).(subFolders(i).name(1:4)) = load('Stats/training_stats.mat').training_stats;
                end
                nMice = nMice + 1;
                
            end
                                
            cd(ori)
        end
    end
    
    data.nMice = nMice;

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
