function [obj, varargout] = mtsess(varargin)
%@mtsess Constructor function for mtsess class
%   OBJ = mtsess(varargin)
%
%   OBJ = mtsess('auto') attempts to create a mtsess object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on mtsess %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%example [as, Args] = mtsess('save','redo')
%
%dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
				'ObjectLevel','Session','RequiredFile','ID*.mat', ...
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
Args.classname = 'mtsess';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'mtsess';

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
    matfilelist    = dir(fullfile(foldername,'ID*.mat'));
    Treadmill_Data = load(fullfile(foldername,matfilelist(1).name));
    Treadmill_Data = Treadmill_Data.Data;
    cd(ori);

%% Fluorescence 
F_raw                 =  double(readNPY(fullfile(foldername,'F.npy')));
spikes_raw            =  double(readNPY(fullfile(foldername,'spks.npy')));
iscell                =  double(readNPY(fullfile(foldername,'iscell.npy')));
Fc_raw                =  double(readNPY(fullfile(foldername,'Fc.npy')));
Fbase_raw             =  double(readNPY(fullfile(foldername,'Fbase.npy')));
iscellno              =  find(iscell(:,1)==1);
F0                    =  F_raw(iscellno,:);
Fc0                   =  Fc_raw(iscellno,:);
Fbase0                =  Fbase_raw(iscellno,:);
dF_F0                 =  double(Fc0./Fbase0);
spikes0               =  spikes_raw(iscellno,:);  
[nNeuron, nImg]       =  size(F0);
spikes0_corrected = load(fullfile(foldername,'spikes0_corrected.mat'));
spikes0_corrected = spikes0_corrected.spikes0_corrected;
dF_F0_corrected_ = load(fullfile(foldername,'dF_F0_corrected.mat'));
dF_F0_corrected_ = dF_F0_corrected_.dF_F0_corrected;
%% Treadmill data
water0                =  Treadmill_Data.Water_real;
pos_cum0              =  Treadmill_Data.Distance_real;
lick                  =  Treadmill_Data.Lick > 3.3; % Use 3.3 as threshold
watertimes            =  squeeze(find(diff(water0)>3));
trialNo               =  zeros(size(pos_cum0));
pos_trial0            =  zeros(size(pos_cum0));
pos_trial0_normalized =  zeros(size(pos_cum0));
for      i     =  1 : (size(watertimes)-1)   
    trialNo((watertimes(i)+1):watertimes(i+1)) = i;
    pos_tmp    =  pos_cum0((watertimes(i)+1):watertimes(i+1)) - pos_cum0(watertimes(i)); 
    pos_trial0((watertimes(i)+1):watertimes(i+1))  =  pos_tmp;
    pos_trial0_normalized((watertimes(i)+1):watertimes(i+1)) = (pos_tmp-min(pos_tmp))/(max(pos_tmp)-min(pos_tmp)); 
end
%% Match Timestamp_treadmill with Timestamp_res
Timestamp_treadmill   =  Treadmill_Data.Timestamp_real;
Timestamp_res         =  Treadmill_Data.Timestamp_resonant_syn;
ResStr                =  Treadmill_Data.Resonant_start;
ResEnd                =  Treadmill_Data.Resonant_End;
ResStr1               =  circshift(ResStr,-1);
ResEnd1               =  circshift(ResEnd,-1);
% Finding timestamps of peaks
s0                    =  find(ResStr<2 & ResStr1 >2);
e0                    =  find(ResEnd<2 & ResEnd1 >2);
s                     =  s0(1:length(e0));
e                     =  e0;
Imgstate0             = zeros(length(Timestamp_res),1);
for i = 1:length(e0)
    Imgstate0(s(i):e(i)) = 5;
end
Imgstate                         =    interp1(Timestamp_res,  Imgstate0,   Timestamp_treadmill, 'linear');
Imgstate1                        =    circshift(Imgstate,-1);
s1                               =    find(Imgstate<2 & Imgstate1 >2);
Imgstate2                        =    circshift(Imgstate,1);
s2                               =    find(Imgstate<2 & Imgstate2 >2);

tsFindex0                        =    round (s1 + 1/2 * (s2-s1)); %Time of Images
tsFindex1                        =    tsFindex0; %(1:nImg);  
usedtrialNo                      =    length(watertimes);
ok                               =    (tsFindex1>=(watertimes(1)+1) & tsFindex1<=watertimes(usedtrialNo) ); % Filter out first/last half-trial 
tsFindex                         =    tsFindex1( ok );

tsF                              =    Timestamp_treadmill(tsFindex);
spikes                           =    spikes0(:,ok);
spikes_corrected                 =    spikes0_corrected(:,ok);
dF_F0_corrected                  =    dF_F0_corrected_(:,ok);
F                                =    F0(:,ok);
Fc                               =    Fc0(:,ok);

%% Session Data

% nTrials = trialNo_dsp(end);
nTrials = max(trialNo);
BinSize = Args.BinSize;

% Trial information (per trial): Trial no., Start time, End time, Duration,
%                                Duration spent licking, Mean running speed (cm/s)
%                                Mean acceleration (cm/s^2), Bin Start Idx, Bin End Idx
data_trial = zeros(nTrials, 9);

% Bin information (per bin): Bin id, Trial no., Bin no., Start time, End time, Duration,
%                            Duration spent licking, Mean running speed (cm/s), 
%                            Mean acceleration (cm/s^2), Trial Start Idx, Trial End Idx
data_bin = zeros(BinSize * nTrials, 11);

% get duration spent at each bin
pos_trial_bin = ceil(pos_trial0_normalized*100);
pos_trial_bin(pos_trial_bin == 0) = 1;
new_binNo = ((trialNo - 1) * BinSize) + pos_trial_bin;
water = zeros(size(pos_cum0));
for i = 1:size(watertimes, 1)
    water(watertimes(i)) = 1;
end
session_data_raw = cat(2, Timestamp_treadmill, trialNo, pos_trial_bin, new_binNo, water, lick, pos_trial0_normalized);
exclude_zero_trials_idx = find(session_data_raw(:,2) ~= 0);
session_data_exclude_zero_trials = session_data_raw(exclude_zero_trials_idx, :);
actual_start_time = session_data_exclude_zero_trials(1,1);
session_data_exclude_zero_trials(:,1) = session_data_exclude_zero_trials(:,1) - actual_start_time; % Adjust start time

% data_bin
for i = 1:(nTrials * BinSize)
    data_bin(i,1) = i; % New bin no.
    data_bin(i,2) = fix((i-1)/BinSize)+1; % Trial no.
    if mod(i,BinSize) == 0
        data_bin(i,3) = BinSize;
    else
        data_bin(i,3) = mod(i,BinSize); % New bin no.
    end
    bin_start_idx = find(session_data_raw(:,4) == i, 1, 'first');
    bin_end_idx = find(session_data_raw(:,4) == i, 1, 'last');
    data_bin(i,4) = session_data_raw(bin_start_idx, 1); % Start time
    data_bin(i,5) = session_data_raw(bin_end_idx, 1); % End time
    data_bin(i,6) = data_bin(i,5) - data_bin(i,4); % Duration
    lick_count = nnz(session_data_raw(bin_start_idx:bin_end_idx, 5));
    data_bin(i,7) = lick_count / (bin_end_idx - bin_start_idx + 1); % Duration spent licking
    data_bin(i,8) = 2.2 / data_bin(i,6); % Mean running speed (Bin length = 2.2cm)
    if i == 1
        data_bin(i,9) = 0;
    else
        data_bin(i,9) = data_bin(i,8) / data_bin(i,6); % Mean acceleration
    end
    data_bin(i,10) = bin_start_idx;
    data_bin(i,11) = bin_end_idx;
    
    if mod(i,100) == 0
        fprintf('Processing bin %i\n',i);
    end
end
data_bin(:,4) = data_bin(:,4) - actual_start_time; % Adjust start time
data_bin(:,5) = data_bin(:,5) - actual_start_time; % Adjust end time

% compute trial durations from the timestamps
trial_start = [];
trial_end = [];
for i = 1:nTrials
    trial_start{i} = data_bin(find(data_bin(:,2) == i, 1, 'first'), 4);
    trial_end{i} = data_bin(find(data_bin(:,2) == i, 1, 'last'), 5);
end
trial_start = cell2mat(trial_start);
trial_end = cell2mat(trial_end);
trial_dur = trial_end - trial_start; % in seconds

% data_trial
data_trial(:,2) = trial_start; % Start time
data_trial(:,3) = trial_end; % End time
data_trial(:,4) = trial_dur; % Duration
for i = 1:nTrials
    data_trial(i,1) = i; % Trial no.
    trial_start_idx = find(session_data_exclude_zero_trials(:,2) == i, 1, 'first');
    trial_end_idx = find(session_data_exclude_zero_trials(:,2) == i, 1, 'last');
    lick_count = nnz(session_data_exclude_zero_trials(trial_start_idx:trial_end_idx, 5));
    data_trial(i,5) = lick_count / (trial_end_idx - trial_start_idx + 1); % Duration spent licking
    data_trial(i,6) = 220 / data_trial(i,4); % Mean running speed (Track length = 220cm)
    if i == 1
        data_trial(i,7) = 0;
    else
        data_trial(i,7) = (data_trial(i,6) - data_trial(i-1,6)) / data_trial(i,4); % Mean acceleration
    end
    data_trial(i,8) = trial_start_idx;
    data_trial(i,9) = trial_end_idx;
end

velocity = [];
acceleration = [];
for i = 2:size(Timestamp_treadmill, 1)
    velocity(i) = (session_data_raw(i,7) - session_data_raw(i-1,7)) / (session_data_raw(i,1) - session_data_raw(i-1,1));
    acceleration(i) = velocity(i) / (session_data_raw(i,1) - session_data_raw(i-1,1));
end
velocity = 220*velocity';
velocity(velocity < -1000) = 0;
acceleration = acceleration';
acceleration(acceleration < -1000) = 0;
session_data_raw = cat(2, session_data_raw, velocity, acceleration);

TrialTime_all_exclude_zero_trials = session_data_raw(exclude_zero_trials_idx, (1:3));
TrialTime_idx = zeros(nTrials, 2);
for i = 1:nTrials
    TrialTime_idx(i,1) = find(TrialTime_all_exclude_zero_trials(:,2) == i, 1, 'first');
    TrialTime_idx(i,2) = find(TrialTime_all_exclude_zero_trials(:,2) == i, 1, 'last');
end

% Velocity and Acceleration
% trial = 1;
% plot_start_idx = TrialTime_idx(trial,1);
% plot_end_idx = TrialTime_idx(trial,2);
% figure;
% subplot(2,1,1);
% plot(session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,7)*220); title('Distance'); xlabel('Time (s)'); ylabel('Distance (cm)'); xline(data_trial(1,3));
% subplot(2,1,2);
% plot(session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,8)); title('Velocity'); xlabel('Time (s)'); ylabel('Velocity (cm/s)'); xline(data_trial(1,3));

% Sliding Filter
windowSize = 501; 
padSize = floor(windowSize/2);
velocity_padded = padarray(velocity, padSize);
velocity_averaged = zeros(size(velocity,1), 1);
acceleration_padded = padarray(acceleration, padSize);
acceleration_averaged = zeros(size(acceleration,1), 1);
for i = 1+padSize:size(velocity, 1)+padSize
    velocity_averaged(i - padSize) = mean(velocity_padded((i - padSize):(i + padSize)));
    acceleration_averaged(i - padSize) = mean(acceleration_padded((i - padSize):(i + padSize)));
end
% figure;
% subplot(2,1,1);
% plot(session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,7)*220); title('Distance'); xlabel('Time (s)'); ylabel('Distance (cm)'); xline(data_trial(1,3));
% subplot(2,1,2);
% plot(session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), velocity_averaged(plot_start_idx:plot_end_idx)); title('Velocity'); xlabel('Time (s)'); ylabel('Velocity (cm/s)'); xline(data_trial(1,3));
% subplot(3,1,3);
% plot(session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), acceleration_averaged(plot_start_idx:plot_end_idx)); title('Acceleration'); xlabel('Time (s)'); ylabel('Acceleration (cm/s^2)'); xline(data_trial(1,3));

velocity_averaged_ = cat(2, session_data_raw(:,1:3), velocity_averaged);
velocity_averaged_ = velocity_averaged_(exclude_zero_trials_idx, :);
velocity_averaged_filt = velocity_averaged_;
velocity_averaged_filt(velocity_averaged_filt(:,4) < Args.ThresVel, :) = NaN;
% figure; plot(velocity_averaged_sliced(plot_start_idx:plot_end_idx,1), velocity_averaged_sliced(plot_start_idx:plot_end_idx, 4)); title('Velocity'); xlabel('Time (s)'); ylabel('Velocity (cm/s)'); xline(data_trial(1,3));
% figure; plot(velocity_averaged_positive(plot_start_idx:plot_end_idx,1), velocity_averaged_positive(plot_start_idx:plot_end_idx, 4)); title('Velocity'); xlabel('Time (s)'); ylabel('Velocity (cm/s)'); xline(data_trial(1,3));

session_data_exclude_zero_trials = cat(2, session_data_exclude_zero_trials, velocity_averaged_(:,4), velocity_averaged_filt(:,4));

velocity_binned = zeros(BinSize, 1);
for i = 1:BinSize
    temp_arr = velocity_averaged_filt(find(velocity_averaged_filt(:,3) == i),4);
    velocity_binned(i) = mean(temp_arr, 'omitnan'); % Positive velocity values only
end

% figure; plot(1:size(velocity_binned(:,1),1), velocity_binned(:,1), 'b'); title('Mean Velocity along Track Distance'); xlabel('Distance (Bin)'); ylabel('Mean Velocity (cm/s)');

% figure; histogram(velocity_averaged_(:,4)); title('Velocity Distribution'); xlabel('Velocity (cm/s)'); ylabel('Count');

% Licking and Water
% figure; plot(session_data_exclude_zero_trials(:,1), session_data_exclude_zero_trials(:,5)*2, 'b'); xlabel('Time (s)');% xline(data_trial(:,3));
% hold on
% plot(session_data_exclude_zero_trials(:,1), session_data_exclude_zero_trials(:,6), 'r'); title('Licking'); xlabel('Time (s)');% xline(data_trial(:,3));

% Time distribution of licks (Method 1 - Bin 50 to Bin 50)
mid_bin_timestamps = zeros(nTrials, 2);
for i = 1:nTrials
    temp_idx = find(session_data_exclude_zero_trials(:,2) == i & session_data_exclude_zero_trials(:,3) == 50, 1, 'first');
    mid_bin_timestamps(i,1) = temp_idx;
    mid_bin_timestamps(i,2) = session_data_exclude_zero_trials(temp_idx,1);
end
water_timestamps = Timestamp_treadmill(watertimes) - actual_start_time;

lick_count_idx = find(diff(session_data_exclude_zero_trials(:,6)) == 1) + 1;
lick_count = session_data_exclude_zero_trials(:,6) < 0; % Resetting array to 0
lick_count(lick_count_idx) = 1;
lick_count = cat(2, session_data_exclude_zero_trials(:,1:3), lick_count);
lick_count_vel_filt = lick_count;
lick_count_vel_filt(find(velocity_averaged_(:,4) < Args.ThresVel),4) = 0; % Ignore lick when smoothed velocity <= 5 cm/s
session_data_exclude_zero_trials = cat(2, session_data_exclude_zero_trials, lick_count(:,4), lick_count_vel_filt(:,4));

lick_timestamps = session_data_exclude_zero_trials(logical(lick_count_vel_filt(:,4)),[1:4]);
lick_timestamps_spliced = lick_timestamps(lick_timestamps(:,1) > mid_bin_timestamps(1,2) & lick_timestamps(:,1) < mid_bin_timestamps(end,2), :); % splice licks from before 1st midpoint and after last midpoint
lick_timestamps_adjusted = lick_timestamps_spliced(:,1); 

for i = 1:(nTrials-1)
    temp_condition = find(lick_timestamps_spliced(:,1) >= mid_bin_timestamps(i,2) & lick_timestamps_spliced(:,1) <= mid_bin_timestamps(i+1,2));
    lick_timestamps_adjusted(temp_condition) = lick_timestamps_spliced(temp_condition) - water_timestamps(i+1);
end
edges = floor(min(lick_timestamps_adjusted)):0.5:ceil(max(lick_timestamps_adjusted));
% edges = -2:0.01:2;
% figure; histogram(lick_timestamps_adjusted, edges); title('Time distribution of licks'); xlabel('Time (s)');

% Bin distribution of licks
% figure; histogram(lick_timestamps_spliced(:,2), 1:101); title('Bin distribution of licks'); xlabel('Bin No.');

lick_count_binned = zeros(BinSize*nTrials,1);
for i =  1:(BinSize*nTrials)
    temp_arr = lick_timestamps(find(lick_timestamps(:,4) == i)); %Count licks per bin
    lick_count_binned(i,1) = nnz(temp_arr);
end
lick_count_binned = reshape(lick_count_binned,[BinSize,nTrials]);

lick_binned = zeros(BinSize, 3);
for i =  1:BinSize
    temp_arr = lick_timestamps(find(lick_timestamps(:,3) == i)); %Count licks across each bin
    lick_binned(i,1) = nnz(temp_arr);
    temp_arr2 = data_bin(find(data_bin(:,3) == i), 6); %Sum duration of all bin instances
    lick_binned(i,2) = sum(temp_arr2);
end
lick_binned(:,3) = lick_binned(:,1) ./ lick_binned(:,2); % Calculate lick rate
% figure; plot(1:size(lick_binned,1), lick_binned(:,3), 'b'); title('Lick Rate along Track Distance'); xlabel('Distance (Bin)'); ylabel('Lick Rate (s^-1)');

% figure;
% plot(session_data_exclude_zero_trials(:,1), session_data_exclude_zero_trials(:,7)*220, 'b'); title('Distance'); xlabel('Time (s)'); ylabel('Distance (cm)'); xline(data_trial(:,3), 'r');
% hold on
% plot(session_data_exclude_zero_trials(:,1), lick_count_vel_filt(:,4)*250, 'g');

sessionMidpoint = find(session_data_exclude_zero_trials(:,2) == ceil(nTrials/2)+1, 1, 'first');

data.Args = Args;
data.nTrials = nTrials;
data.nNeuron = nNeuron;
data.data_bin = data_bin;
data.data_trial = data_trial;
data.session_data_raw = session_data_raw;
data.session_data_exclude_zero_trials = session_data_exclude_zero_trials;
data.sessionMidpoint = sessionMidpoint;
data.actual_start_time = actual_start_time;
data.velocity_binned = velocity_binned;
data.lick_timestamps = lick_timestamps;
data.lick_timestamps_spliced = lick_timestamps_spliced;
data.lick_timestamps_adjusted = lick_timestamps_adjusted;
data.lick_count_binned = lick_count_binned;
data.lick_binned = lick_binned;
data.tsF = tsF;
data.spikes = spikes;
data.spikes_corrected = spikes_corrected;
data.dF_F0 = dF_F0;
data.dF_F0_corrected = dF_F0_corrected;
data.F = F;
data.Fc = Fc;

% create nptdata so we can inherit from it
data.numSets = 0;
data.Args = Args;
n = nptdata(1,0,pwd);
d.data = data;
obj = class(d,Args.classname,n);
saveObject(obj,'ArgsC',Args);

% Make cell directories and save spiketimes
mkdir cells
cd cells

ori2 = pwd;

for i = 1:5 %nNeuron
    cell_folderName = strcat('cell', num2str(i));
    mkdir(cell_folderName);
    cd(cell_folderName);
    spiketrain = spikes_corrected(i,:);
    save('spiketrain.mat', 'spiketrain', '-v7.3');
    
    timestamp_dsp = zeros(size(tsFindex));
    for ii = 1:size(tsFindex)
       timestamp_dsp(ii) = Timestamp_treadmill(tsFindex(ii)) - actual_start_time;
    end
    spiketimes = timestamp_dsp(find(spiketrain)); % Get timestamp of spikes idx
    save('spiketimes.mat', 'spiketimes', '-v7.3');
    
    cd(ori2)
end

cd(ori)

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
