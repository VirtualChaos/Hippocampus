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
				'BinSize',100, 'ThresVel',0, 'Training',0, 'RedoSess',0, 'RedoNeural',0, 'Match',0, 'Slice',[0,0]);
            
Args.flags = {'Auto','ArgsOnly'};
% Specify which arguments should be checked when comparing saved objects
% to objects that are being asked for. Only arguments that affect the data
% saved in objects should be listed here.
Args.DataCheckArgs = {'BinSize', 'ThresVel','Training',};                         

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
    
    if Args.Match
        match_file = load("Matches/ID68_match_" + Args.Slice(1) + "_" + Args.Slice(2) + ".mat");
    end
    
    ori = pwd;
    data.origin = {pwd};
    
    files = dir;
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    sessionCombined = {};
    neuralCombined = {};
    session_no = 0;
    
    if ~isstring(Args.Slice)
        valid_flag = 1;
    else
        valid_flag = 0;
    end
    
    for i = 1:size(subFolders,1)
        current_folder = subFolders(i).name;
        if Args.Training
            foldername_check = current_folder(1) == 'I' && current_folder(2) == 'D';
        else
            foldername_check = current_folder(1) == '2' && current_folder(2) == '0' && current_folder(3) == '2';
        end
        
        if Args.Training & foldername_check
            temp_numbers = regexp(current_folder,'[0-9]','match');
            training_current_folder = strjoin(temp_numbers(3:10),'');
        end
               
        if foldername_check && isstring(Args.Slice)
            %'Training',1,'Slice',["20220421","20220422"]
            if ~Args.Training
                current_date = datetime(str2num(current_folder(1:4)), str2num(current_folder(5:6)), str2num(current_folder(7:8)));
            else
                current_date = datetime(str2num(training_current_folder(1:4)), str2num(training_current_folder(5:6)), str2num(training_current_folder(7:8)));
            end
            slice_start = convertStringsToChars(Args.Slice(1));
            slice_start_date = datetime(str2num(slice_start(1:4)), str2num(slice_start(5:6)), str2num(slice_start(7:8)));
            slice_end = convertStringsToChars(Args.Slice(2));
            slice_end_date = datetime(str2num(slice_end(1:4)), str2num(slice_end(5:6)), str2num(slice_end(7:8)));
            
            if days(current_date - slice_start_date) == 0
                valid_flag = 1;
            elseif days(current_date - slice_end_date) > 0
                valid_flag = 0;
            end
        end
        
        if  foldername_check && valid_flag
                                 
            cd(current_folder);
            
            if ~isempty(dir('ID*.mat'))
                
                tic
                session_no = session_no + 1;
                fprintf("Session %d: %s\n", session_no, current_folder);
                if Args.Training
                    if ~isfile('mtsess.mat')
                        sessionData = mtsess('auto','redo','save','Spikes',0).data;
                    else
                        if Args.RedoSess
                            sessionData = mtsess('auto','redo','save','Spikes',0).data;
                        else
                            sessionData = mtsess('auto','Spikes',0).data;
                        end
                    end
                else
                    if Args.RedoSess
                        sessionData = mtsess('auto','redo','save').data;
                    else
                        sessionData = mtsess('auto','save').data;
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
                    if Args.RedoNeural
                        neuralData = mtneuraldata('auto','redo','save').data;
                    else
                        if isfile('mtneuraldata.mat')
                            neuralData = mtneuraldata('auto').data;
                        else
                            neuralData = mtneuraldata('auto','save').data;
                        end
                    end                    
                    neuralData.date = current_folder;
                    neuralCombined.("s" + session_no) = neuralData;
                    cd ..
                    
                end
                toc
                                
            end
                
        end
                        
            cd(ori);
    end

    if Args.Match
        
        valid = false(size(match_file.roiMatchData.allSessionMapping,1),size(match_file.roiMatchData.allSessionMapping,2));
        for sess_no = 1:size(match_file.roiMatchData.allSessionMapping,2)
            valid(:,sess_no) = match_file.roiMatchData.allSessionMapping(:,sess_no) < neuralCombined.("s" + sess_no).nNeuron;
        end
        
        match_file.roiMatchData.allSessionMapping = match_file.roiMatchData.allSessionMapping(~any(valid == 0,2),:);
        
        data.roiMatchData = match_file.roiMatchData;
        
    end
    
    sessionDays = zeros(size(fieldnames(sessionCombined),1),1);
    datetime_origin = datetime(str2num(sessionCombined.s1.date(1:4)), str2num(sessionCombined.s1.date(5:6)), str2num(sessionCombined.s1.date(7:8)));
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
