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
				'BinSize',100, 'ThresVel',0, 'Training',0);
            
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
        if Args.Training
            foldername_check = current_folder(1) == 'I' && current_folder(2) == 'D';
        else
            foldername_check = current_folder(1) == '2' && current_folder(2) == '0' && current_folder(3) == '2';
        end
        if  foldername_check
            
            cd(current_folder);
            
            tic
            session_no = session_no + 1;
            fprintf("Session %d: %s\n", session_no, current_folder);
            if Args.Training
                if ~isfile('mtsess.mat') || (Args.RedoLevels)
                    sessionData = mtsess('auto','redo','save','Spikes',0).data;
                else
                    sessionData = mtsess('auto','Spikes',0).data;
                end
            else
                if ~isfile('mtsess.mat') || (Args.RedoLevels)
                    sessionData = mtsess('auto','redo','save').data;
                else
                    sessionData = mtsess('auto').data;
                end
            end
            
            if Args.Training
                temp_numbers = regexp(current_folder,'[0-9]','match');
                sessionData.date = strjoin(temp_numbers(3:10),'');
            else
                sessionData.date = current_folder;
            end
            sessionCombined.("s" + session_no) = sessionData;
            
            if ~Args.Training
                
                cd cells
                neuralData = mtneuraldata('auto','redo','save').data;
                neuralData.date = current_folder;
                neuralCombined.("s" + session_no) = neuralData;
                cd ..
                
            end
            toc
            
            cd(ori);
        end
    end
    
    sessionDays = zeros(size(fieldnames(sessionCombined),1),1);
    datetime_origin = datetime(str2num(sessionCombined.s1.date(1:4)), str2num(sessionCombined.s1.date(5:6)), str2num(sessionCombined.s1.date(7:8)));;
    sessionDays(1) = 1;
    for sess_no = 2:size(fieldnames(sessionCombined),1)
        sessionData = sessionCombined.("s" + sess_no);
        datetime_current = datetime(str2num(sessionData.date(1:4)), str2num(sessionData.date(5:6)), str2num(sessionData.date(7:8)));
        sessionDays(sess_no) = days(datetime_current - datetime_origin) + 1;
    end
    
    data.sessionCombined = sessionCombined;
    data.neuralCombined = neuralCombined;
    data.sessionDays = sessionDays;
    
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
