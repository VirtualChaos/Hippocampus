function [obj, varargout] = mtmice(varargin)
%@mtmice Constructor function for mtmice class
%   OBJ = mtmice(varargin)
%
%   OBJ = mtmice('auto') attempts to create a mtmice object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on mtmice %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%example [as, Args] = mtmice('save','redo')
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
Args.classname = 'mtmice';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'mtmice';

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
    sessionCombined = {};
    neuralCombined = {};
    session_no = 0;
    
    for i = 1:size(subFolders,1)
        current_folder = subFolders(i).name;
        if  current_folder(1) == '2' && current_folder(2) == '0' && current_folder(3) == '2'
            
            cd(current_folder);
            
            if isfile('mtsess.mat')
                tic
                session_no = session_no + 1;
                fprintf("Session %d: %s\n", session_no, current_folder);
                sessionData = mtsess('auto').data;
                sessionData.date = current_folder;
                sessionCombined.("s" + session_no) = sessionData;
                
                cd cells
                neuralData = mtneuraldata('auto','save').data;
                neuralData.date = current_folder;
                neuralCombined.("s" + session_no) = neuralData;
                cd ..
                toc
            end
            
            cd(ori);
        end
    end
    
    data.sessionCombined = sessionCombined;
    data.neuralCombined = neuralCombined;

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
