function [obj, varargout] = mtneural(varargin)
%@mtneural Constructor function for mtneural class
%   OBJ = mtneural(varargin)
%
%   OBJ = mtneural('auto') attempts to create a mtneural object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on mtneural %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%example [as, Args] = mtneural('save','redo')
%
%dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
				'ObjectLevel','Session','RequiredFile','cells_list.txt', 'NumericArguments', [], ...
				'BinSize',100, 'ThresVel',1, 'FiltFrac',0.6, 'PeakThreshold',3, 'TrialPeakThreshold',1, 'TrialPeakCheckFrac',2/3);
            
Args.flags = {'Auto','ArgsOnly'};
% Specify which arguments should be checked when comparing saved objects
% to objects that are being asked for. Only arguments that affect the data
% saved in objects should be listed here.
Args.DataCheckArgs = {'BinSize', 'ThresVel', 'FiltFrac', 'PeakThreshold', 'TrialPeakThreshold', 'TrialPeakCheckFrac'};                         

[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'mtneural';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'mtneural';

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
    
    sessionData = mtsess('Auto', varargin{:});
    sessionData = sessionData.data;
    cellcombinedData = mtsess('Auto', varargin{:});
    cellcombinedData = cellcombinedData.data;
    
    
    cells_list = string(importdata(Args.RequiredFile));
    
    isplacecell = zeros(sessionData.nTrials,5);
    shuffled_sic_values = zeros(sessionData.nTrials, 10000);
    spiketrain_combined = zeros(sessionData.nTrials, size(sessionData.F,2));
    mtplacefield_failed = [];
    binFiringRate_combined = zeros(sessionData.nTrials, sessionData.nTrials, sessionData.Args.BinSize);
    
    % Place field identification
    for cell_idx = 1:sessionData.nTrials
        if mod(cell_idx,50) == 0
            fprintf("Current Cell: %d / %d; Estimated Time Remaining: %.2f seconds\n", cell_idx, sessionData.nTrials, (toc/cell_idx) * sessionData.nTrials - toc);
        end
        
        
        
    end
    
    % Organising data to store
    tic
    for cell_idx = 1:sessionData.nTrials
        
        cell_no = cells_list(cell_idx);
        cd(strcat(ori, '/', cell_no));
        

        
        try
            cellData = load('mtcell.mat');
            cellData = cellData.mtcell.data;
        catch
            cellData = mtcell('auto','save','redo').data;
        end
        
        try
            mtplacefield('auto','save','redo', ... 
                            'FiltFrac',Args.FiltFrac,'PeakThreshold',Args.PeakThreshold,'TrialPeakThreshold',Args.TrialPeakThreshold,'TrialPeakCheckFrac',Args.TrialPeakCheckFrac);
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
        
        binFiringRate_combined(cell_idx,:,:) = cellData.binFiringRate;
        
        try
            placefieldData = mtplacefield('auto','save','redo', ... 
                            'FiltFrac',Args.FiltFrac,'PeakThreshold',Args.PeakThreshold,'TrialPeakThreshold',Args.TrialPeakThreshold,'TrialPeakCheckFrac',Args.TrialPeakCheckFrac).data;
        catch
            mtplacefield_failed = cat(1,mtplacefield_failed, cell_idx);
            placefieldData = {};
        end
        
        spiketimes = load('spiketimes.mat');
        spiketimes = spiketimes.spiketimes;
        spiketrain = load('spiketrain.mat');
        spiketrain = spiketrain.spiketrain;
        
        data.placefieldData.("n" + sprintf('%04d',cell_no_)) = placefieldData;
        data.spiketimes.("n" + sprintf('%04d',cell_no_)) = spiketimes;
%         try
%             data.spiketimes.("n" + sprintf('%04d',cell_no_)) = spiketimes;
%         catch
%             data.spiketimes.("n" + sprintf('%04d',cell_no_)) = [];
%         end
        %data.spiketrain.("n" + sprintf('%04d',cell_no_)) = spiketrain;
        spiketrain_combined(cell_idx,:) = spiketrain;
               
        cd(ori);
    end
    data.mtplacefield_failed = mtplacefield_failed;
    isplacecell(:,5) = isplacecell(:,2) > prctile(shuffled_sic_values,95,'all');
    data.spiketrain = spiketrain_combined;
    data.binFiringRate = binFiringRate_combined;
    
    data.isplacecell = sortrows(isplacecell);
    data.shuffled_sic_values = shuffled_sic_values;
    data.cellData = orderfields(data.cellData);
    data.placefieldData = orderfields(data.placefieldData);
    
    data.sessData.data_trial = sessionData.data_trial;
    data.sessData.tsF = sessionData.tsF;
    data.sessData.dF_F0_corrected = sessionData.dF_F0_corrected;
    
    rmse = zeros(length(cells_list), 2);
    nrmse_mean = zeros(length(cells_list), 2);
    nrmse_stdev = zeros(length(cells_list), 2);
    auc = zeros(length(cells_list), 2);
    log_likelihood = zeros(length(cells_list), 2);
    aic = zeros(length(cells_list), 2);
    bic = zeros(length(cells_list), 2);
    threshold = zeros(length(cells_list), 5);
    placecellStats = zeros(length(cells_list), 4);
    placefieldStats = [];
    fns = fieldnames(data.placefieldData);
    for i = 1:length(cells_list)
        rmse(i,1) = i;
        nrmse_mean(i,1) = i;
        nrmse_stdev(i,1) = i;
        auc(i,1) = i;
        log_likelihood(i,1) = i;
        aic(i,1) = i;
        bic(i,1) = i;
        threshold(i,1) = i;
        placecellStats(i,1) = i;
        if ~isempty(data.placefieldData.(fns{i})) && logical(data.isplacecell(i,4))
            rmse(i,2) = data.placefieldData.(fns{i}).rmse;
            nrmse_mean(i,2) = data.placefieldData.(fns{i}).nrmse_mean;
            nrmse_stdev(i,2) = data.placefieldData.(fns{i}).nrmse_stdev;
            auc(i,2) = data.placefieldData.(fns{i}).auc;
            try
                log_likelihood(i,2) = data.placefieldData.(fns{i}).log_likelihood;
            catch
                data.placefieldData.(fns{i})
            end
            aic(i,2) = data.placefieldData.(fns{i}).aic;
            bic(i,2) = data.placefieldData.(fns{i}).bic;
            threshold(i,2) = data.placefieldData.(fns{i}).mean_threshold;
            threshold(i,3) = data.placefieldData.(fns{i}).std_threshold;
            threshold(i,4) = data.placefieldData.(fns{i}).cv_threshold;
            threshold(i,5) = data.placefieldData.(fns{i}).vmr_threshold;
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
            log_likelihood(i,2) = NaN;
            aic(i,2) = NaN;
            bic(i,2) = NaN;
            threshold(i,2:end) = NaN;
            placecellStats(i,2) = NaN;
            placecellStats(i,3) = mean(data.cellData.(fns{i}).maps_adsm,'omitnan'); % Mean firing rate
            placecellStats(i,4) = std(data.cellData.(fns{i}).maps_adsm,'omitnan'); % Stdev firing rate
        end
    end
    
%     rmse_sort = sortrows(rmse,2);
%     nrmse_mean_sort = sortrows(nrmse_mean,2);
%     nrmse_stdev_sort = sortrows(nrmse_stdev,2);
%     auc_sort = sortrows(auc,2);
    
    data.rmse = rmse;
    data.nrmse_mean = nrmse_mean;
    data.nrmse_stdev = nrmse_stdev;
    data.auc = auc;
    data.log_likelihood = log_likelihood;
    data.aic = aic;
    data.bic = bic;
    data.threshold = threshold;
    data.placecellStats = placecellStats;
    data.placefieldStats = placefieldStats;
    data.nTrials = sessionData.nTrials;
    data.nNeuron = sessionData.nNeuron;
    
    % binFiringRate_combined = neuralData.data.binFiringRate;
    % spiketrain_combined = neuralData.data.spiketrain;
    
    % Bayesian Decoding
%     tsFindex = sessionData.tsFindex;
%     tsF = sessionData.tsF;
%     nTrials = sessionData.nTrials;
%     velThreshold = 1; % 1 cm/s
%     
%     tsFindex_first_idx_exclude_zero_trials = find(sessionData.session_data_raw(:,2) == 1,1,'first');
%     tsFindex_last_idx_exclude_zero_trials = find(sessionData.session_data_raw(:,2) == nTrials,1,'last');
%     ok_ezt = (tsFindex > tsFindex_first_idx_exclude_zero_trials) & (tsFindex < tsFindex_last_idx_exclude_zero_trials);
%     tsFindex_ezt = tsFindex - tsFindex_first_idx_exclude_zero_trials + 1; % Adjust idx spdwhen shifting from session_raw to session_ezt 
%     tsFindex_ezt = tsFindex_ezt(ok_ezt);
%     vel_ave_downsample = sessionData.session_data_exclude_zero_trials(tsFindex_ezt,8);
%     tsFindex_ezt_velthreshold = tsFindex_ezt(vel_ave_downsample > velThreshold);
%     ok_ezt_velthreshold = ok_ezt(vel_ave_downsample > velThreshold);
%         
%     tsF_downsample = sessionData.session_data_exclude_zero_trials(tsFindex_ezt_velthreshold,1) + sessionData.actual_start_time;
%     trialNo_downsample = sessionData.session_data_exclude_zero_trials(tsFindex_ezt_velthreshold,2);
%     pos_downsample = sessionData.session_data_exclude_zero_trials(tsFindex_ezt_velthreshold,7);
%     pos_bin_downsample = sessionData.session_data_exclude_zero_trials(tsFindex_ezt_velthreshold,3);
%     vel_ave_downsample = sessionData.session_data_exclude_zero_trials(tsFindex_ezt_velthreshold,8);
%     vel_ave_filt_downsample = sessionData.session_data_exclude_zero_trials(tsFindex_ezt_velthreshold,9);
%     lick_downsample = sessionData.session_data_exclude_zero_trials(tsFindex_ezt_velthreshold,10);
%     lick_burst_downsample = sessionData.lick_burst(tsFindex_ezt_velthreshold);
%     
%     spiketrain_combined_ezt = spiketrain_combined(:,ok_ezt_velthreshold);
% 

       
%     %% Firingmap calculation
%     %  Firingrate_trial      :  neuron * trialNo * firingrate (real spks number at per time point)
%     %  Firingrate_trial_sm   :  neuron * trialNo * firingrate (normalized and smoothed)
%     %  Firingrate_average_sm :  neuron * firingrate           (average each trial, normalized and smoothed)
%     
%     F_trial = [];
%     pos_trial_downsample_trial = [];
%     tsF_trial = [];
%     for i = 1:size(spiketrain_combined_ezt,1)
%         for j = 1:max(trialNo_downsample)
%             id                            = find(trialNo_downsample == j);
%             Ftmp                          = spiketrain_combined_ezt(i,id);
%             F_trial{i}{j}                 = Ftmp;
%             pos_trial_downsample_trial{j} = pos_bin_downsample(id);
%             tsF_trial{j}                  = tsF_downsample(id);
%             for ilocation = 1:100
%                 if isempty(find(pos_trial_downsample_trial{j} == ilocation));
%                     spatialoccupancy_trial(i,j,ilocation) = 0;
%                     Firingrate_trial(i,j,ilocation)       = 0;
%                 else
%                     spatialoccupancy_trial(i,j,ilocation) = length(find(pos_trial_downsample_trial{j} == ilocation))/length(tsF_trial{j});
%                     Firingrate_trial(i,j,ilocation)       = sum(Ftmp(find(pos_trial_downsample_trial{j} == ilocation)))./(length(find(pos_trial_downsample_trial{j}  == ilocation)));
%                 end
%             end
%             Firingrate_trial_sm0(i,j,:) = smooth(squeeze(Firingrate_trial(i,j,:)),3);
%                         
%         end
%         firingtmp                        = Firingrate_trial_sm0(i,:,:);
%         Firingrate_trial_sm(i,:,:)       = firingtmp;
%         
%         firingavgtmp                     = squeeze(mean(Firingrate_trial(i,:,:),2));
%         firingavgtmp                     = smooth(firingavgtmp,3);
%         Firingrate_average_sm(i,:)       = firingavgtmp;
%     end
%     
%     %% Prior information
%     % Prior: P_pos
%     for ilocation = 1:100
%         ok_iloc = (pos_bin_downsample == ilocation);
%         P_pos(ilocation) = sum(ok_iloc) / length(pos_bin_downsample);
%     end
% 
%     % Prior: meanF_bin_odd: spikes in each bin at each trial/timepoints at each trial; Average odd trials
%     ok_prior_trial     = 1:2:nTrials;
%     meanF_bin_odd      =  mean(Firingrate_trial(1:end,ok_prior_trial,:),2);  % use all neurons and only odd trials for training
%     meanF_bin_odd      =  (squeeze(meanF_bin_odd))' / 100; % Mean firing rates across odd trials
% 
%     tau                 = 0.5; % Duration of time window(second), do 1 decoding at each time window
%     ts_interval         = tsF(2) - tsF(1);
%     tau_bin             = round(tau./ts_interval);
% 
%     ts_decoded     = [];
%     pos_decoded    = [];
%     trial_decoded  = [];
% 
%     for i_Trial = 1:nTrials
% 
%         ok_trial      =  trialNo_downsample == i_Trial;
%         tsF_trial     =  tsF_downsample(ok_trial);
%         spks_trial    =  spiketrain_combined_ezt(:, ok_trial);
%         Trial_bin     =  sum(ok_trial);
%         TwinNo        =  ceil(Trial_bin/tau_bin);
% 
%         for i_Twin = 1:TwinNo
%             st    =  (i_Twin - 1) * tau_bin + 1;
%             ed    =  i_Twin * tau_bin;
%             if ed > Trial_bin
%                 ed =  Trial_bin;
%             end
%             if (tsF_trial(ed) - tsF_trial(st)) < 1 % Time window span not exceed 1second
%                 ts_decoded     =  cat(1, ts_decoded, mean(tsF_trial(st:ed)));
%                 nowF           =  mean(spks_trial(1:end, st:ed),2)';
%                 nowF_mat       =  repmat(nowF, sessionData.Args.BinSize, 1);
% 
%                 % Key formula of Bayesian decoding: P(pos|nowf) -> Pspatial * (meanF(pos))^nowf * exp(-tau*f(pos))
%                 P_pos_spks     =   sum(log((meanF_bin_odd+eps).^nowF_mat),2) - tau * sum(meanF_bin_odd, 2);
%                 %figure(1);plot(P_pos_spks)
%                 [~, idpeak]    =  max(P_pos_spks);
%                 pos_decoded    =  cat(1, pos_decoded,   idpeak);
%                 trial_decoded  =  cat(1, trial_decoded, i_Trial * ones(length(idpeak),1));
%             end
%         end
%     end
% 
%     %% Decoded error at even trials calculation
%     % Real position at each decoded timepoint
%     pos_true = interp1(tsF_downsample, pos_bin_downsample, ts_decoded, 'linear') *2.2; % cm
%     pos_decoded = pos_decoded *2.2; % cm
% 
%     % Decoding error: mean absolute difference
%     ok_even = rem(trial_decoded,2) == 0;
%     ok_odd = rem(trial_decoded,2) == 1;
%     decoded_error_even_matrix = abs(pos_decoded(ok_even) - pos_true(ok_even));
%     pos_err_even_mean  = mean(decoded_error_even_matrix);
%     pos_err_even_median  = median(decoded_error_even_matrix);
% 
%     decoded_error_even_matrix_circ = 220-abs(pos_decoded(ok_even) -pos_true(ok_even));
%     derror = min([decoded_error_even_matrix decoded_error_even_matrix_circ],[],2);
%     pos_err_even_mean_circ  = mean(derror);
%     pos_err_even_median_circ  = median(derror);
% 
%     pos_err = [pos_err_even_mean,pos_err_even_median,pos_err_even_mean_circ,pos_err_even_median_circ];
% 
%     figure;
%     subplot(2,1,1)
%     plot(ts_decoded(:),pos_true(:));
%     plot(sessionData.session_data_exclude_zero_trials(:,1) + sessionData.actual_start_time, sessionData.session_data_exclude_zero_trials(:,7)*220,'k','LineWidth',2);hold on;
%     plot(ts_decoded,pos_decoded,'o','MarkerFaceColor','r','MarkerSize',10);
%     hold on;
% 
%     subplot(2,1,2)
%     histogram(derror,'Normalization','probability'); title('Decoding error distribution'); xlabel('Distance (cm)'); ylabel('Count');


%% Bayesian decoding (original code)
cd ..

spdthresh = 1;
dF_F0_corrected = load('dF_F0_corrected.mat').dF_F0_corrected;
spikes0_corrected = load('spikes0_corrected.mat').spikes0_corrected;
spikes0_corrected = abs(spikes0_corrected);
matfilelist    = dir('ID*.mat');
Treadmill_Data = load(matfilelist(1).name);
Treadmill_Data = Treadmill_Data.Data;

% Fluorescence: spikes0
% F_raw                 =  double(readNPY(fullfile(foldername,'F.npy')));
% spikes_raw            =  double(readNPY(fullfile (foldername,'spks.npy')));
iscell                =  double(readNPY('iscell.npy'));
iscellno              =  find(iscell(:,1)==1);
F0                    =  dF_F0_corrected; %F_raw(iscellno,:);
spikes0               =  spikes0_corrected; %spikes_raw(iscellno,:);
[nNeuron, nImg]       =  size(F0);

% Lick signal
Lick_raw              =  Treadmill_Data.Lick;
% Match Timestamp_treadmill with Timestamp_res
Timestamp_treadmill   =  Treadmill_Data.Timestamp_real;
Timestamp_res         =  Treadmill_Data.Timestamp_resonant_syn;
ResStr                =  Treadmill_Data.Resonant_start;
ResEnd                =  Treadmill_Data.Resonant_End;
ResStr1               =  circshift(ResStr,-1);
ResEnd1               =  circshift(ResEnd,-1);
s0                    =  find(ResStr<2 & ResStr1 >2);
e0                    =  find(ResEnd<2 & ResEnd1 >2);
s                     =  s0(1:length(e0));
e                     =  e0;
Imgstate0             =  zeros(length(Timestamp_res),1);
for i = 1:length(e0)
    Imgstate0(s(i):e(i)) = 5;
end
Imgstate              =  interp1(Timestamp_res,  Imgstate0,   Timestamp_treadmill, 'linear');
Imgstate1             =  circshift(Imgstate,-1);
s1                    =  find(Imgstate<2 & Imgstate1 >2);
Imgstate2             =  circshift(Imgstate,1);
s2                    =  find(Imgstate<2 & Imgstate2 >2);
tsFindex_raw          =  round (s1 + 1/2 * (s2-s1));

%% Treadmill data:pos_trial0; pos_trial0_normalized
water0                =  Treadmill_Data.Water_real;

%Miss water point
% T_misswater           = ;%20210508_ID38:1127.77;20210607_ID45:2038.82; %20210604_ID45:1348;
% ok_misswater          = (Timestamp_treadmill > T_misswater) & (Timestamp_treadmill <(T_misswater+1));
% water0(ok_misswater)  = 5;

pos_cum0              =  Treadmill_Data.Distance_real;
watertimes            =  squeeze(find(diff(water0)>2.25));
trialNo               =  zeros(size(pos_cum0));
pos_trial0            =  zeros(size(pos_cum0));
pos_trial0_normalized =  zeros(size(pos_cum0));
for      i     =  1 : (size(watertimes)-1);
    trialNo((watertimes(i)+1):watertimes(i+1)) =i;
    pos_tmp                =  pos_cum0((watertimes(i)+1):watertimes(i+1)) - pos_cum0(watertimes(i));
    pos_trial0 ((watertimes(i)+1):watertimes(i+1)) ...
        =  pos_tmp;
    pos_trial0_normalized((watertimes(i)+1):watertimes(i+1)) ...
        = (pos_tmp-min(pos_tmp))/(max(pos_tmp)-min(pos_tmp));
end
pos_tmp_start                              = pos_cum0(1:watertimes(1)) + max(pos_trial0) - max(pos_cum0(1:watertimes(1)));
pos_trial0(1:watertimes(1))                = pos_tmp_start;
pos_trial0_normalized(1:watertimes(1))     =(pos_tmp_start-0)/(max(pos_trial0)-0);
pos_tmp_end                                = pos_cum0(watertimes(end):end) - pos_cum0(watertimes(end));
pos_trial0(watertimes(end):end)            = pos_tmp_end;
pos_trial0_normalized(watertimes(end):end) =(pos_tmp_end-0)/(max(pos_trial0)-0);

%% Treadmill data: spd
step = 500;
pos_cum0_1 = [zeros(1,step),pos_cum0',zeros(1,step)];
spd = zeros(length(Timestamp_treadmill),1);
for i = (step/2 +1): (length(Timestamp_treadmill) - step/2)
    spd(i) =  pos_cum0_1(i+step/2) - pos_cum0_1(i-step/2);
end
spd = spd'.*220/1300/(step/5000); % cm/s

%% Downsample pos and spd to F
tsFindex0                        =    tsFindex_raw(1:nImg);
tsF0                             =    Timestamp_treadmill(tsFindex0);
trialNo0                         =    trialNo(tsFindex0);
nTrial0                          =    max(trialNo0);
pos_downsample0                  =    pos_trial0_normalized(tsFindex0);
pos_bin_downsample0              =    round(pos_downsample0*100);
spd_downsample0                  =    spd(tsFindex0);
Lick0                            =    Lick_raw(tsFindex0);
%%
ok_rest = ones(1,length(spd_downsample0));
for i = 1:100
    spd_shift_1(i,:)= [spd_downsample0(i+1:end) zeros(1,i)];
    spd_shift_2(i,:) = [zeros(1,i) spd_downsample0(1:end-i)];
    ok_rest = ok_rest & spd_downsample0 == 0 & spd_shift_1(i,:)==0 & spd_shift_2(i,:) == 0;
end

%% Define the trials and speed thresh for calculating firingmap
ok_trial                    =   trialNo0  > 0 ;
ok_spd                      =   spd_downsample0 > spdthresh;
ok_firemap                  =   ok_trial & ok_spd';
tsF_firemap                 =   tsF0(ok_firemap);
spikes_firemap              =   spikes0(:,ok_firemap);
trialNo_dsp_firemap         =   trialNo0(ok_firemap);
pos_trial_bin_dsp_firemap   =   pos_bin_downsample0(ok_firemap);
spd_dsp_firemap             =   spd_downsample0(ok_firemap);

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
            if isempty(find(pos_trial_downsample_trial{j} == ilocation));                                              
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

% Prior:
% P_pos         : nBin * 1,       probability of mouse show up at each bin
% meanF_bin_odd : nBin * nNeuron, mean firing rate of each neuron at each bin
% nowF_mat      : nBin * nNeuron, current firing rate (nowF)

% Decoding formula:
% P(pos|nowf) -> P_pos * (meanF_bin_odd)^nowF_mat * exp(-tau*meanF_bin_odd)

ok = trialNo_dsp_firemap <=  max(trialNo_dsp_firemap);

tsF_used                  =  tsF_firemap(ok);
spikes_used               =  spikes_firemap(:,ok);
trialNo_dsp_used          =  trialNo_dsp_firemap(ok);
pos_trial_bin_dsp_used    =  pos_trial_bin_dsp_firemap(ok);
spd_dsp_used              =  spd_dsp_firemap(ok);

[nNeuron,~]               =  size(spikes_used);
nTrial                    =  max(trialNo_dsp_used);
Neuron_list               =  1:nNeuron;
nBin                      =  max(pos_trial_bin_dsp_used);

%% Prior information
% Prior: P_pos
for ilocation = 1:100
    ok_iloc                = pos_trial_bin_dsp_used == ilocation;
    P_pos(ilocation)   = sum(ok_iloc)/length(pos_trial_bin_dsp_firemap);
end
%P_pos = ones(100,1);
% Prior: meanF_bin_odd: spikes in each bin at each trial/timepoints at each trial; Average odd trials
ok_prior_trial     = 1:2:nTrial;
meanF_bin_odd      =  mean(Firingrate_trial(Neuron_list,ok_prior_trial,:),2);  % use all neurons and only odd trials for training
meanF_bin_odd      =  (squeeze(meanF_bin_odd))';

tau                 = 0.5;                       % Duration of time window(second), do 1 decoding at each time
ts_interval         = tsF0(2) - tsF0(1);
tau_bin             = round(tau ./ts_interval);

ts_decoded     = [];
pos_decoded    = [];
trial_decoded  = [];
%% Decoding process: each trial -> each 0.5s time window
for i_Trial = 1:nTrial
    
    ok_trial      =  trialNo_dsp_used == i_Trial;
    tsF_trial     =  tsF_used(ok_trial);
    spks_trial    =  spikes_used(:, ok_trial);
    Trial_bin     =  sum(ok_trial);
    TwinNo        =  ceil(Trial_bin/tau_bin);
    
    for i_Twin = 1:TwinNo
        st    =  (i_Twin - 1) * tau_bin + 1;
        ed    =  i_Twin * tau_bin;
        if ed > Trial_bin
            ed =  Trial_bin;
        end
        if (tsF_trial(ed) - tsF_trial(st)) < 1 % Time window span not exceed 1second
            ts_decoded     =  cat(1, ts_decoded, mean(tsF_trial(st:ed)));
            nowF           =  mean(spks_trial(:, st:ed),2)';
            nowF_mat       =  repmat(nowF, nBin, 1);
            
            % Key formula of Bayesian decoding: P(pos|nowf) -> Pspatial * (meanF(pos))^nowf * exp(-tau*f(pos))
            P_pos_spks     =   sum(log((meanF_bin_odd+eps).^nowF_mat),2)- tau * sum(meanF_bin_odd, 2);
            %figure(1);plot(P_pos_spks)
            [~, idpeak]    =  max(P_pos_spks);
            pos_decoded    =  cat(1, pos_decoded,   idpeak);
            trial_decoded  =  cat(1, trial_decoded, i_Trial * ones(length(idpeak),1));
        end
    end
end

%% Decoded error at even trials calculation
% Real position at each decoded timepoint
pos_true = interp1(tsF_used, pos_trial_bin_dsp_used, ts_decoded, 'linear') *2.2; % cm
pos_decoded = pos_decoded *2.2; % cm

% Decoding error: mean absolute difference
ok_even = rem(trial_decoded,2) == 0;
ok_odd = rem(trial_decoded,2) == 1;
decoded_error_even_matrix = abs(pos_decoded(ok_even) - pos_true(ok_even));
pos_err_even_mean  = mean(decoded_error_even_matrix);
pos_err_even_median  = median(decoded_error_even_matrix);

decoded_error_even_matrix_circ = 220-abs(pos_decoded(ok_even) -pos_true(ok_even));
derror = min([decoded_error_even_matrix decoded_error_even_matrix_circ],[],2);
pos_err_even_mean_circ  = mean(derror);
pos_err_even_median_circ  = median(derror);

pos_err = [pos_err_even_mean,pos_err_even_median,pos_err_even_mean_circ,pos_err_even_median_circ];

% figure;
% subplot(2,1,1)
% plot(ts_decoded(:),pos_true(:));
% plot(Timestamp_treadmill, pos_trial0_normalized*220,'k','LineWidth',2);hold on;
% plot(ts_decoded,pos_decoded,'o','MarkerFaceColor','r','MarkerSize',10);
% hold on;
% 
% subplot(2,1,2)
% histogram(derror,'Normalization','probability'); title('Decoding error distribution'); xlabel('Distance (cm)'); ylabel('Count');

data.pos_err = pos_err;
data.derror = derror;


    cd(ori);
    
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
