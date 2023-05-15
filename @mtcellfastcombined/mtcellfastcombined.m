function [obj, varargout] = mtcellfastcombined(varargin)
%@mtcellfastcombined Constructor function for mtcellfastcombined class
%   OBJ = mtcellfastcombined(varargin)
%
%   OBJ = mtcellfastcombined('auto') attempts to create a mtcellfastcombined object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on mtcellfastcombined %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%example [as, Args] = mtcellfastcombined('save','redo')
%
%dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
				'ObjectLevel','Cell','RequiredFile','spiketrain.mat', ...
				'BinSize',100, 'ThresVel',1, 'NumShuffles',10000, 'ShuffleLimits',[0.1 0.9], ...
                'AdaptiveSmooth',1, 'Alpha',50);
            
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
Args.classname = 'mtcellfastcombined';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'mtcellfastcombined';

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Cell class
    ori = pwd;
    data.origin = {pwd};
    sessionData = mtsess('auto', varargin{:});
    sessionData = sessionData.data;
    toKeep = {'session_data_exclude_zero_trials', 'sessionMidpoint', 'actual_start_time', ...
                'velocity_averaged', 'nTrials', 'nNeuron', 'tsFindex', 'tsF', 'session_data_raw'};
    f = fieldnames(sessionData);
    toRemove = f(~ismember(f,toKeep));
    sessionData = rmfield(sessionData,[toRemove]);
    cd('cells');
    vel_threshold_handler = string(importdata("vel_threshold_handler.txt"));
    Args.VelThresholdHandler = vel_threshold_handler;
    cd(ori);
    
    spiketrain_all = load(Args.RequiredFile);
    spiketrain_all = spiketrain_all.spiketrain;
    
    %% Downsample treadmill data
    session_data_raw_dsp = sessionData.session_data_raw(sessionData.tsFindex,:);
    trialNo0                         =    session_data_raw_dsp(:,2);%trialNo(tsFindex);
    pos_bin_downsample0              =    session_data_raw_dsp(:,3);%round(pos_downsample0*100);
    spd_downsample0                  =    session_data_raw_dsp(:,10);%spd(tsFindex);
    
    %% Define the trials and speed thresh for calculating firingmap
    ok_trial                    =   session_data_raw_dsp(:,2)  > 0 ;
    if vel_threshold_handler == "rest"
        ok_spd = spd_downsample0 <= Args.ThresVel;
    elseif vel_threshold_handler == "zero"
        ok_spd = spd_downsample0 == 0;
    else
        ok_spd = spd_downsample0 > Args.ThresVel;
    end
    ok_firemap                  =   ok_trial & ok_spd;
    tsF_firemap                 =   sessionData.tsF(ok_firemap);
    spikes_firemap              =   spiketrain_all(:,ok_firemap);
    trialNo_dsp_firemap         =   trialNo0(ok_firemap);
    pos_trial_bin_dsp_firemap   =   pos_bin_downsample0(ok_firemap);
    %spd_dsp_firemap             =   spd_downsample0(ok_firemap);
    
    F_trial = [];
    pos_trial_downsample_trial = [];
    tsF_trial = [];
    for i = 1:size(spikes_firemap,1)
        for j = 1:max(trialNo_dsp_firemap)
            id                            = find(trialNo_dsp_firemap == j);
            Ftmp                          = spikes_firemap(i,id);
            F_trial{i}{j}                 = Ftmp;
            pos_trial_downsample_trial{j} = pos_trial_bin_dsp_firemap(id);
            tsF_trial{j}                  = tsF_firemap(id);
            for ilocation = 1:100
                if isempty(find(pos_trial_downsample_trial{j} == ilocation))
                    spatialoccupancy_trial(i,j,ilocation) = 0;
                    Firingrate_trial(i,j,ilocation)       = 0;
                else
                    spatialoccupancy_trial(i,j,ilocation) = length(find(pos_trial_downsample_trial{j} == ilocation))/length(tsF_trial{j});
                    Firingrate_trial(i,j,ilocation)       = sum(Ftmp(find(pos_trial_downsample_trial{j} == ilocation)))./(length(find(pos_trial_downsample_trial{j}  == ilocation)));
                end
            end
            Firingrate_trial_sm0(i,j,:) = smooth(squeeze(Firingrate_trial(i,j,:)),3);
            
        end
        firingtmp                        = Firingrate_trial_sm0(i,:,:);
        Firingrate_trial_sm(i,:,:)       = firingtmp;
        
        firingavgtmp                     = squeeze(mean(Firingrate_trial(i,:,:),2));
        firingavgtmp                     = smooth(firingavgtmp,3);
        Firingrate_average_sm(i,:)       = firingavgtmp;
    end
    
    data.binFiringRate = Firingrate_trial;

    data.nTrials = sessionData.nTrials;
    data.nNeuron = sessionData.nNeuron;
    data.spiketrain = spikes_firemap;
    
    % create nptdata so we can inherit from it
    data.numSets = 1;
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
