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
				'BinSize',100, 'ThresVel',5, 'NumShuffles',1000, 'ShuffleLimits',[0.1 0.9], ...
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
            Firingrate_trial_sm0(i,j,:) = Smooth(squeeze(Firingrate_trial(i,j,:)),3);
            
        end
        firingtmp                        = Firingrate_trial_sm0(i,:,:);
        Firingrate_trial_sm(i,:,:)       = firingtmp;
        
        firingavgtmp                     = squeeze(mean(Firingrate_trial(i,:,:),2));
        firingavgtmp                     = Smooth(firingavgtmp,3);
        Firingrate_average_sm(i,:)       = firingavgtmp;
    end
    
    data.binFiringRate = Firingrate_trial;

    data.nTrials = sessionData.nTrials;
    data.nNeuron = sessionData.nNeuron;
    data.spiketrain = spikes_firemap;
    
    %% Select place cell by SIC
    nshuffle = Args.NumShuffles;
    [~,nImg_firemap] = size(spikes_firemap);
    shuffle = randi([round(nImg_firemap*0) round(nImg_firemap*1)],nshuffle,1);     
    for ilocation = 1:Args.BinSize
        IDbin{ilocation} = find(pos_trial_bin_dsp_firemap ==ilocation);
        spatialoccupancy(ilocation) = length(IDbin{ilocation}) / length(pos_trial_bin_dsp_firemap);
    end
    j = 0;
    SICVec = zeros(1,sessionData.nNeuron);
    SICVecPos = zeros(sessionData.nNeuron,Args.BinSize);
    for i = 1:sessionData.nNeuron
        spiketmp          = spikes_firemap(i,:);
        spikes_repeat     = repmat(spiketmp,nshuffle,1);
        spikes_shuffle    = reshape(cell2mat(arrayfun(@(k) circshift(spikes_repeat(k,:),shuffle(k)),1:nshuffle,'uni',0)),nImg_firemap,nshuffle)';
%         spikes_shuffle    = spikes_shuffle;
        for ilocation = 1:Args.BinSize
             Firingrate(i,ilocation)         = sum(spikes_firemap(i,IDbin{ilocation} ))./length(IDbin{ilocation} );    
             Firingrate_shuffle(:,ilocation) = sum(spikes_shuffle(:,IDbin{ilocation}),2)./(size(spikes_shuffle(:,IDbin{ilocation}),2) );
        end         
        [SICVec(i), SICVecPos(i,:)] = SpatialInfo(spatialoccupancy',Firingrate(i,:));
        SICmatrix = cell2mat(arrayfun(@(k) SpatialInfo(spatialoccupancy',Firingrate_shuffle(k,:)),1:nshuffle,'uni',0));  
        if length(find(SICmatrix<SICVec(i)))  > nshuffle *0.95  
            j = j+1;
            placecell_passSIC(j) = i;        
        end
        SICmatrix_all{i} = SICmatrix;
    end
    SIC.SICVec              = SICVec;
    SIC.SICmatrix_all       = SICmatrix_all;
    SIC.placecell_passSIC   = placecell_passSIC;
    SIC.Firingrate          = Firingrate;
    
    %% Select place cell
    placecellPeakAmp   = 0.3; % Place cell criteria: peak amplitude must > 0.3 of maximum
    placecellFieldMin  = 0.5; % Place cell criteria: field region > 0.5 of peak
    placecellInOut     = 3;   % Place cell criteria: In out activity ratio is 3
    placecellGoodtrial = 1/2; % Place cell criteria: Peak inside place field more that half of trials

    PeakAmp         = placecellPeakAmp;
    Fieldthreshold  = placecellFieldMin;
    In_out_ratio    = placecellInOut;
    GoodtrialNo     = ceil(sessionData.nTrials*placecellGoodtrial);
    
    Field = {};
    FieldPeakID = {};
    FieldSize = {};
    for iNeuron =  1:sessionData.nNeuron
        if ~ismember(iNeuron,placecell_passSIC)
            FieldNo(iNeuron)       = 0;
        else
            FiringNeuron = Firingrate_average_sm(iNeuron,:);
            Firingtmp  = Firingrate_average_sm(iNeuron,:);
            posmaptrial = squeeze(Firingrate_trial(iNeuron,:,:));
            iField = 1;
            while true
                if iField >= 4 % criteria
                    break;
                end
                [peak,idx] = max(Firingtmp);
                if peak <= PeakAmp*max(FiringNeuron) % criteria
                    break;
                end
                [~,Peakx] = ind2sub(size(Firingtmp),idx);
                field     = FindField(Firingtmp,Peakx,1,peak*Fieldthreshold,true,0);   % criteria
                [fieldboundary,~]          = FieldBoundaries(field,true,0);
                fieldXmin                  = fieldboundary(1,1);
                fieldXmax                  = fieldboundary(1,2);
                if fieldXmax > fieldXmin
                    fieldbins = [fieldXmin:fieldXmax];
                else
                    fieldXmin = fieldXmin + 1;
                    fieldXmax = fieldXmax - 1;
                    fieldbins = [fieldXmin:100 1:fieldXmax];
                end
                
                activityinfield            = mean(Firingtmp(fieldbins),'omitnan');
                %activityoutfield           = mean(Firingtmp([1:fieldXmin fieldXmax:100]));
                activityoutfield           = mean(Firingtmp(find(~ismember(1:100,fieldbins))),'omitnan');
                SingleFieldinoutratio      = activityinfield/activityoutfield;
                if SingleFieldinoutratio <= In_out_ratio % criteria
                    break
                end
                %[~,maxind]  = max(posmaptrial,[],2);
%                 Goodtrial = find(maxind<=fieldXmax & maxind>=fieldXmin);
                [~,maxind]  = maxk(posmaptrial,5,2);
                Goodtrial = find(any(ismember(maxind,fieldbins),2));
                if (length(Goodtrial) <= GoodtrialNo)   % criteria
                    break
                end
                % Check field overlap
                bin_buffer = 3;
                fieldbins_check1 = 0;
                fieldbins_check2 = 0;
                if iField >= 2
                    [fieldboundary,~] = FieldBoundaries(Field{iNeuron,iField-1},true,0);
                    fieldXmin                  = fieldboundary(1,1) - bin_buffer;
                    fieldXmax                  = fieldboundary(1,2) + bin_buffer;
                    if fieldXmin <= 0
                        fieldXmin = fieldXmin + 100;
                    end
                    if fieldXmax >= 101
                        fieldXmax = fieldXmax - 100;
                    end
                    if fieldXmax > fieldXmin
                        fieldbins_check1 = [fieldXmin:fieldXmax];
                    else
%                         fieldXmin = fieldXmin + 1;
%                         fieldXmax = fieldXmax - 1;
                        fieldbins_check1 = [fieldXmin:100 1:fieldXmax];
                    end
                end
                if iField == 3
                    [fieldboundary,~] = FieldBoundaries(Field{iNeuron,iField-2},true,0);
                    fieldXmin                  = fieldboundary(1,1) - bin_buffer;
                    fieldXmax                  = fieldboundary(1,2) + bin_buffer;
                    if fieldXmin <= 0
                        fieldXmin = fieldXmin + 100;
                    end
                    if fieldXmax >= 101
                        fieldXmax = fieldXmax - 100;
                    end
                    if fieldXmax > fieldXmin
                        fieldbins_check2 = [fieldXmin:fieldXmax];
                    else
%                         fieldXmin = fieldXmin + 1;
%                         fieldXmax = fieldXmax - 1;
                        fieldbins_check2 = [fieldXmin:100 1:fieldXmax];
                    end
                end
                if ~ismember(Peakx,fieldbins_check1) & ~ismember(Peakx,fieldbins_check2)
                    Field      {iNeuron,iField} = field;
                    FieldPeakID{iNeuron,iField} = Peakx;
                    FieldSize  {iNeuron,iField} = sum(field(:));
                    iField = iField +1;
                end
                Firingtmp(field) = NaN;
                posmaptrial(:,field) = NaN;
                if all(isnan(Firingtmp))
                    break;
                end
            end
            FieldNo(iNeuron)       = iField-1;
        end
    end
    %clearvars -except spikes0 F0 nNeuron ntrial Field FieldPeakID FieldSize FieldNo Firingrate_trial_sm Firingrate_average_sm
    %% placecell fraction
    PCNo = 0;
    for iNeuron = 1:sessionData.nNeuron
        if FieldNo(iNeuron) >0
            PCNo = PCNo +1;
            ID_PC  (PCNo) = iNeuron;
            Peak_PC(PCNo) = FieldPeakID{iNeuron,1};
        end
    end
    if PCNo ~= 0
        [~,IDtmp]     = sort(Peak_PC);
        ID_PC_sorted  = ID_PC(IDtmp);
        PCNo = 0;
        PC1No = 0;
        PC2No = 0;
        PC3No = 0;
        FdNo  = 0;
        for iNeuron = 1:sessionData.nNeuron
            if FieldNo(iNeuron) == 1
                PCNo         = PCNo +1;
                PC_ID(PCNo)  = iNeuron;
                PC_FdNo(PCNo)= 1;
                PC1No        = PC1No +1;
                PC1_ID(PC1No)= iNeuron;
                FdNo         = FdNo +1;
                Fdsize(FdNo) = FieldSize{iNeuron,1};
            elseif FieldNo(iNeuron) == 2
                PCNo         = PCNo +1;
                PC_ID(PCNo)  = iNeuron;
                PC_FdNo(PCNo)= 2;
                PC2No        = PC2No +1;
                PC2_ID(PC2No)= iNeuron;
                FdNo         = FdNo +1;
                Fdsize(FdNo) = FieldSize{iNeuron,1};
                FdNo         = FdNo +1;
                Fdsize(FdNo) = FieldSize{iNeuron,2};
            elseif FieldNo(iNeuron) == 3
                PCNo         = PCNo +1;
                PC_ID(PCNo)  = iNeuron;
                PC3No        = PC3No +1;
                PC_FdNo(PCNo)= 3;
                PC3No        = PC3No +1;
                PC3_ID(PC3No)= iNeuron;
                FdNo         = FdNo +1;
                Fdsize(FdNo) = FieldSize{iNeuron,1};
                FdNo         = FdNo +1;
                Fdsize(FdNo) = FieldSize{iNeuron,2};
                FdNo         = FdNo +1;
                Fdsize(FdNo) = FieldSize{iNeuron,3};
            end
        end
        PCfraction(1) = PCNo/sessionData.nNeuron;
        PCfraction(2) = PC1No/sessionData.nNeuron;
        PCfraction(3) = PC2No/sessionData.nNeuron;
        PCfraction(4) = PC3No/sessionData.nNeuron;
        PCfraction(5) = PC1No/PCNo;
        PCfraction(6) = PC2No/PCNo;
        PCfraction(7) = PC3No/PCNo;
    else
        PCfraction   = zeros(1,7);
        PC_ID        = [];
        FieldSize    = [];
        FieldNo      = [];
        Fdsize       = [];
        ID_PC_sorted = [];
        FieldPeakID = [];
    end
    %% Output results
    data.PCfraction   = PCfraction;   % place cell fraction
    data.PC_ID        = PC_ID;        % place cell ID
    data.ID_PC_sorted = ID_PC_sorted; % place cell ID sorted by 1st peak position
    data.FieldNo      = FieldNo;      % place cell field number
    data.FieldSize    = FieldSize;    % place cell field size
    data.Fdsize       = Fdsize;       % all field size
    data.SIC          = SIC;
    data.FieldPeakID  = FieldPeakID;
    
    %% Bayesian Decoding Shuffle
    nNeuron_Baye   = 100; % fixed neuron number in each session
    %% Prior information: P_pos; meanF_bin_odd(spikes in each bin at each trial/timepoints at each trial; Average odd trials)
    for ilocation = 1:100
        ok_iloc           = pos_trial_bin_dsp_firemap == ilocation;
        P_pos(ilocation)  = sum(ok_iloc)/length(pos_trial_bin_dsp_firemap);
    end
    % P_pos = ones(100,1);
    %% Define decoding time points
    tau = 0.5; % Duration of time window(second), do 1 decoding at each time window
    tsF_interval = sessionData.tsF(2) - sessionData.tsF(1);
    tau_bin = round(tau./tsF_interval);
    %% Shuffle process
    for ishuffle = 1:100
        if mod(ishuffle,10) == 0
            fprintf('Processing Shuffle %i\n',ishuffle);
        end
        norder = randperm(sessionData.nNeuron);
        n_select = norder(1:100);
        spikes_firemap = spiketrain_all(n_select,ok_firemap);
        
        pos_trial_downsample_trial = cell(sessionData.nTrials);
        tsF_trial                  = cell(sessionData.nTrials);
        spatialoccupancy_trial     = zeros(nNeuron_Baye,sessionData.nTrials,Args.BinSize);
        Firingrate_trial      = zeros(nNeuron_Baye,sessionData.nTrials,Args.BinSize);
        Firingrate_trial_sm   = zeros(nNeuron_Baye,sessionData.nTrials,Args.BinSize);
        Firingrate_average    = zeros(nNeuron_Baye,Args.BinSize);
        Firingrate_average_sm = zeros(nNeuron_Baye,Args.BinSize);
        
        for i = 1:nNeuron_Baye
            for j = 1:sessionData.nTrials
                id                            = find(trialNo_dsp_firemap == j);
                Ftmp                          = spikes_firemap(i,id);
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
                Firingrate_trial_sm(i,j,:) = Smooth(squeeze(Firingrate_trial(i,j,:)),3);
            end
            Firingrate_average(i,:)        = squeeze(mean(Firingrate_trial(i,:,:),2));
            Firingrate_average_sm(i,:)     = Smooth(Firingrate_average(i,:),3);
        end
        
        %% Calculate meanF_bin_odd
        ok_prior_trial =  1:2:sessionData.nTrials;
        meanF_bin_odd  =  mean(Firingrate_trial_sm(1:nNeuron_Baye,ok_prior_trial,:),2);  % use all neurons and only odd trials for training
        meanF_bin_odd  =  (squeeze(meanF_bin_odd))';
        %% Decoding process: each trial -> each 0.5s time window
        ts_decoded     = [];
        pos_decoded    = [];
        trial_decoded  = [];
        bin_actual     = [];
        for i_Trial = 1:sessionData.nTrials
            ok_trial      =  trialNo_dsp_firemap == i_Trial;
            tsF_trial     =  tsF_firemap(ok_trial);
            bin_trial     =  pos_trial_bin_dsp_firemap(ok_trial);
            spks_trial    =  spikes_firemap(:,ok_trial);
            Trial_bin     =  sum(ok_trial);
            TimewinNo     =  ceil(Trial_bin/tau_bin);
            for i_Timewin = 1:TimewinNo
                st    =  (i_Timewin - 1) * tau_bin + 1;
                ed    =  i_Timewin * tau_bin;
                if ed > Trial_bin
                    ed =  Trial_bin;
                end
                if (tsF_trial(ed) - tsF_trial(st)) < 1 % Time window span not exceed 1second
                    ts_decoded     =  cat(1, ts_decoded, mean(tsF_trial(st:ed)));
                    bin_actual     =  cat(1, bin_actual, mode(bin_trial(st:ed)));
                    nowF           =  mean(spks_trial(:, st:ed),2)';
                    nowF_mat       =  repmat(nowF, Args.BinSize, 1);
                    % Key formula of Bayesian decoding: P(pos|nowf) -> Pspatial * (meanF(pos))^nowf * exp(-tau*f(pos))
                    P_pos_spks     =  sum(log((meanF_bin_odd+eps).^nowF_mat),2)- tau * sum(meanF_bin_odd, 2);
                    [~, idpeak]    =  max(P_pos_spks);
                    pos_decoded    =  cat(1, pos_decoded,   idpeak);
                    trial_decoded  =  cat(1, trial_decoded, i_Trial * ones(length(idpeak),1));
                end
            end
        end
        %% Decoded error at even trials calculation
        pos_true    = 2.2 * interp1(tsF_firemap, pos_trial_bin_dsp_firemap, ts_decoded, 'linear'); % Real position at each decoded timepoint (cm)
        pos_decoded_bin = pos_decoded;
        pos_decoded = 2.2 * pos_decoded; % cm
        ok_even = rem(trial_decoded,2) == 0;
        ok_odd  = rem(trial_decoded,2) == 1;
        decoding_error_even       = abs(pos_decoded(ok_even) - pos_true(ok_even));
        decoding_error_even_oppo  = 220 - decoding_error_even;
        decoding_error_even_circ  = min([decoding_error_even decoding_error_even_oppo],[],2);
        pos_err_even_mean         = mean  (decoding_error_even);
        pos_err_even_median       = median(decoding_error_even);
        pos_err_even_mean_circ    = mean  (decoding_error_even_circ);
        pos_err_even_median_circ  = median(decoding_error_even_circ);
        pos_err = [pos_err_even_mean,pos_err_even_median,pos_err_even_mean_circ,pos_err_even_median_circ];
        pos_err_shuffle(ishuffle,:) = pos_err;
        
        pos_decoded_bin_shuffle(ishuffle,:) = pos_decoded_bin;
    end
    
    BayeDecoding_shuffle.ts_decoded                    = ts_decoded;
    BayeDecoding_shuffle.bin_actual                    = bin_actual;
    BayeDecoding_shuffle.pos_decoded_bin_shuffle       = pos_decoded_bin;
    BayeDecoding_shuffle.pos_err_shuffle               = pos_err_shuffle;
    BayeDecoding_shuffle.pos_err_shuffle_mean          = mean(pos_err_shuffle,1);
    
    data.BayeDecoding_shuffle = BayeDecoding_shuffle;

    %%
    % Pairwise correlation analysis
    trialNo       = session_data_raw_dsp(:,2);
    speed_downsample         = spd_downsample0';

    ok_rest = ones(1,length(speed_downsample));
    spd_shift_1 = zeros(50,length(speed_downsample));
    spd_shift_2 = zeros(50,length(speed_downsample));
    for i = 1:50
        spd_shift_1(i,:) = [speed_downsample(i+1:end) zeros(1,i)];
        spd_shift_2(i,:) = [zeros(1,i) speed_downsample(1:end-i)];
        ok_rest          = ok_rest & speed_downsample == 0 & spd_shift_1(i,:)==0 & spd_shift_2(i,:)==0;
    end
    ok_trial       =   trialNo  > 0 ;
    ok_rest        =   ok_rest & ok_trial';
    ok_spd         =   speed_downsample > Args.ThresVel;
    ok_run         =   ok_trial & ok_spd' ;
    
    spikes_rest    =   spiketrain_all(:,ok_rest);
    spikes_run     =   spiketrain_all(:,ok_run);
    
    meanspks_all        = mean(spiketrain_all,2);
    meanspks_res        = mean(spikes_rest,2);
    meanspks_run        = mean(spikes_run,2);
    
    cr_all     = corrcoef(spiketrain_all');
    cr_res     = corrcoef(spikes_rest');
    cr_run     = corrcoef(spikes_run');
    
    ok = triu(true(size(cr_all)),1);
    cr_all_vector     = cr_all(ok);
    cr_res_vector     = cr_res(ok);
    cr_run_vector     = cr_run(ok);
    
    MeanCor     = zeros(3,1);
    MedianCor   = zeros(3,1);
    IDpeak_all  = zeros(3,1);
    Xpeak_all   = zeros(3,1);
    FWHM_all    = zeros(3,1);
    fgauL_sym   = cell(3,1);
    gofgauL_sym = cell(3,1);
    fgauR_sym   = cell(3,1);
    gofgauR_sym = cell(3,1);
    Cor_now     = [cr_all_vector cr_res_vector cr_run_vector]';
    for k = 1:3
        MeanCor(k)    = mean(Cor_now(k,:),'omitnan');
        MedianCor(k)  = median(Cor_now(k,:),'omitnan');
        [N,~,~] = histcounts(Cor_now(k,:),'Normalization','probability','BinWidth',0.0005,'BinLimits',[-1,1]);
        X = -0.9995:0.0005:1;
        IDpeak0         = find(N  == max(N));
        IDpeak          = IDpeak0(1);
        IDpeak_all(k)   = IDpeak;
        Xpeak_all(k)  = X(IDpeak);
        xLeft           = X(1:IDpeak);
        xRight          = X(IDpeak:end);
        Left            = N(1:IDpeak);
        Right           = N(IDpeak:end);
        xLeft_sym       = X(1):0.0005:(2*X(IDpeak)-X(1));
        Left_sym        = [N(1:IDpeak) flip(N(2:IDpeak))];
        xRight_sym      = (2*X(IDpeak)-X(end)):0.0005:X(end);
        Right_sym       = [flip(N(IDpeak:end)) N((IDpeak+1):end)];
        a1              = N(IDpeak);
        b1              = X(IDpeak);
        expression      = strcat(string(a1),'*exp(-((x-',string(b1),')/c1)^2)');
        ft              = fittype( expression, 'independent', 'x', 'dependent', 'y' );
        opts            = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display    = 'Off';
        opts.StartPoint = 0.141886338627215;
        [fgauL_sym{k}, gofgauL_sym{k}]  = fit(xLeft_sym',Left_sym', ft, opts );
        [fgauR_sym{k}, gofgauR_sym{k}]  = fit(xRight_sym',Right_sym', ft, opts );
        c1L             = fgauL_sym{k}.c1;
        muL             = c1L/sqrt(2);
        FWHML           = muL*2.355;
        c1R             = fgauR_sym{k}.c1;
        muR             = c1R/sqrt(2);
        FWHMR           = muR*2.355;
        FWHM_all(k)   = FWHML + FWHMR;
        % plot fit
        %{
    Leftfit_sym  = a1.*exp(-((xLeft_sym-b1)./c1L).^2);
    Rightfit_sym = a1.*exp(-((xRight_sym-b1)./c1R).^2);
    figure(1);
    plot(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'r');hold on;
    plot(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'r');hold on;
    plot(X,N,'k');xlim([-0.1,0.1])
        %}
    end
    if ok_rest == 0
        Xpeak_all(2) = NaN;
        FWHM_all(2) = NaN;
        
    end
    
    data.Cor_Mspks.meanspikes   = [meanspks_all  meanspks_res  meanspks_run ]';
    data.Cor_Mspks.correlations = [cr_all_vector cr_res_vector cr_run_vector]';
    data.Cor_Mspks.Cor_GaussFit = [MeanCor; MedianCor; Xpeak_all; FWHM_all]';
    
    idx = 1:size(cr_res,1)+1:numel(cr_res);
    if sum(isnan(cr_res),'all') ~= 0
        temp_arr = reshape(cr_res.',1,[]);
        nan_list = find(isnan(temp_arr(idx)));
        cr_res(nan_list,:) = 0;
        cr_res(:,nan_list) = 0;
    end
    if sum(isnan(cr_run),'all') ~= 0
        temp_arr = reshape(cr_run.',1,[]);
        nan_list = find(isnan(temp_arr(idx)));
        cr_run(nan_list,:) = 0;
        cr_run(:,nan_list) = 0;
    end
    cr_res(idx) = 0;
    cr_run(idx) = 0;
    
    data.Cor_Mspks.cr_res = cr_res;
    data.Cor_Mspks.cr_run = cr_run;
    
    % Place-Non Place cell connectivity
    temp_cell_list = 1:sessionData.nNeuron;
    pc_list = PC_ID;
    non_pc_list = temp_cell_list(~ismember(temp_cell_list,pc_list));
    
    % SIC list
    SIC_list = [[1:size(SICVec,2)]' SICVec'];
    SIC_list = sortrows(SIC_list,2);
    
    pc_cr_res = cr_res(pc_list,:);
    nonpc_cr_res = cr_res(non_pc_list,:);
    
    pc_pc_cr_res = cr_res(pc_list,pc_list);
    pc_nonpc_cr_res = cr_res(pc_list,non_pc_list);
    
    nonpc_nonpc_cr_res = cr_res(non_pc_list,non_pc_list);
    nonpc_pc_cr_res = cr_res(non_pc_list,pc_list);
    
    high_SIC_list = SIC_list(SIC_list(:,2) > prctile(SIC_list(:,2),75),:);
    low_SIC_list = SIC_list(SIC_list(:,2) < prctile(SIC_list(:,2),25),:);
    
    high_SIC_cr_res = cr_res(high_SIC_list(:,1),:);
    low_SIC_cr_res = cr_res(low_SIC_list(:,1),:);
    
    high_high_SIC_cr_res = cr_res(high_SIC_list(:,1),high_SIC_list(:,1));
    high_low_SIC_cr_res = cr_res(high_SIC_list(:,1),low_SIC_list(:,1));
    
    low_low_SIC_cr_res = cr_res(low_SIC_list(:,1),low_SIC_list(:,1));
    low_high_SIC_cr_res = cr_res(low_SIC_list(:,1),high_SIC_list(:,1));
    
%     temp = proportion_pos_neg(proportion_pos_neg(:,3) > 0.5,1);
    
    MeanCor_all = zeros(1,12);
    MedianCor_all = zeros(1,12);
    FWHM_all = zeros(1,12);
    
    % PC to all; Non-PC to all; High SIC to all; Low SIC to all
    % PC to PC; PC to Non-PC; High SIC to High SIC; High SIC to Low SIC
    % Non-PC to Non-PC; Non-PC to PC; Low SIC to Low SIC; Low SIC to High SIC
    
    if true
%     figure;
%     subplot(3,2,1)
    Cor_now = reshape(pc_cr_res.',1,[]); % Red
    [MeanCor,MedianCor,FWHM,IDpeak,xLeft_sym,Leftfit_sym,xRight_sym,Rightfit_sym] = FWHM_compute(Cor_now);
    MeanCor_all(1) = MeanCor;
    MedianCor_all(1) = MedianCor;
    FWHM_all(1) = FWHM;
    %         plot(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'Color',[0.4660 0.6740 0.1880]);hold on;
    %         plot(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'Color',[0.4660 0.6740 0.1880]);hold on;
%     area(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'EdgeColor',[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840], ...
%         'FaceAlpha',0.5,'LineWidth',2);
%     hold on
%     area(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'EdgeColor',[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840], ...
%         'FaceAlpha',0.5,'LineWidth',2);
%     hold on
%     xline(MedianCor,'Color',[0.6350 0.0780 0.1840],'LineStyle','--');
%     %     plot(X,N,'k');xlim([-0.1,0.1])
%     hold on
    Cor_now = reshape(nonpc_cr_res.',1,[]); % Blue
    [MeanCor,MedianCor,FWHM,IDpeak,xLeft_sym,Leftfit_sym,xRight_sym,Rightfit_sym] = FWHM_compute(Cor_now); % Blue
    MeanCor_all(2) = MeanCor;
    MedianCor_all(2) = MedianCor;
    FWHM_all(2) = FWHM;
    %         plot(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'Color',[0.9290 0.6940 0.1250]);hold on;
    %         plot(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'Color',[0.9290 0.6940 0.1250]);hold on;
%     area(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak), ...
%         'EdgeColor',[0.3010 0.7450 0.9330],'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.5,'LineWidth',2);
%     hold on
%     area(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end), ...
%         'EdgeColor',[0.3010 0.7450 0.9330],'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.5,'LineWidth',2);
%     %     plot(X,N,'k');xlim([-0.1,0.1])
%     hold on
%     xline(MedianCor,'Color',[0.3010 0.7450 0.9330],'LineStyle','--');
%     xlim([-0.1,0.1])
%     xlabel("Correlation Coefficient",'FontSize',18); ylabel("Normalised density",'FontSize',18);
%     title("PC-Non PC")
%     hold off
    
%     subplot(3,2,2)
    Cor_now = reshape(high_SIC_cr_res.',1,[]); % Red
    [MeanCor,MedianCor,FWHM,IDpeak,xLeft_sym,Leftfit_sym,xRight_sym,Rightfit_sym] = FWHM_compute(Cor_now);
    MeanCor_all(3) = MeanCor;
    MedianCor_all(3) = MedianCor;
    FWHM_all(3) = FWHM;
    %         plot(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'Color',[0.4660 0.6740 0.1880]);hold on;
    %         plot(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'Color',[0.4660 0.6740 0.1880]);hold on;
%     area(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'EdgeColor',[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840], ...
%         'FaceAlpha',0.5,'LineWidth',2);
%     hold on
%     area(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'EdgeColor',[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840], ...
%         'FaceAlpha',0.5,'LineWidth',2);
%     hold on
%     xline(MedianCor,'Color',[0.6350 0.0780 0.1840],'LineStyle','--');
%     %     plot(X,N,'k');xlim([-0.1,0.1])
%     hold on
    Cor_now = reshape(low_SIC_cr_res.',1,[]); % Blue
    [MeanCor,MedianCor,FWHM,IDpeak,xLeft_sym,Leftfit_sym,xRight_sym,Rightfit_sym] = FWHM_compute(Cor_now); % Blue
    MeanCor_all(4) = MeanCor;
    MedianCor_all(4) = MedianCor;
    FWHM_all(4) = FWHM;
    %         plot(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'Color',[0.9290 0.6940 0.1250]);hold on;
    %         plot(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'Color',[0.9290 0.6940 0.1250]);hold on;
%     area(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak), ...
%         'EdgeColor',[0.3010 0.7450 0.9330],'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.5,'LineWidth',2);
%     hold on
%     area(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end), ...
%         'EdgeColor',[0.3010 0.7450 0.9330],'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.5,'LineWidth',2);
%     %     plot(X,N,'k');xlim([-0.1,0.1])
%     hold on
%     xline(MedianCor,'Color',[0.3010 0.7450 0.9330],'LineStyle','--');
%     xlim([-0.1,0.1])
%     xlabel("Correlation Coefficient",'FontSize',18); ylabel("Normalised density",'FontSize',18);
%     title("SIC")
%     hold off
    
%     subplot(3,2,3)
    Cor_now = reshape(pc_pc_cr_res.',1,[]); % Red
    [MeanCor,MedianCor,FWHM,IDpeak,xLeft_sym,Leftfit_sym,xRight_sym,Rightfit_sym] = FWHM_compute(Cor_now);
    MeanCor_all(5) = MeanCor;
    MedianCor_all(5) = MedianCor;
    FWHM_all(5) = FWHM;
    %         plot(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'Color',[0.4660 0.6740 0.1880]);hold on;
    %         plot(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'Color',[0.4660 0.6740 0.1880]);hold on;
%     area(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'EdgeColor',[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840], ...
%         'FaceAlpha',0.5,'LineWidth',2);
%     hold on
%     area(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'EdgeColor',[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840], ...
%         'FaceAlpha',0.5,'LineWidth',2);
%     hold on
%     xline(MedianCor,'Color',[0.6350 0.0780 0.1840],'LineStyle','--');
%     %     plot(X,N,'k');xlim([-0.1,0.1])
%     hold on
    Cor_now = reshape(pc_nonpc_cr_res.',1,[]); % Blue
    [MeanCor,MedianCor,FWHM,IDpeak,xLeft_sym,Leftfit_sym,xRight_sym,Rightfit_sym] = FWHM_compute(Cor_now); % Blue
    MeanCor_all(6) = MeanCor;
    MedianCor_all(6) = MedianCor;
    FWHM_all(6) = FWHM;
    %         plot(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'Color',[0.9290 0.6940 0.1250]);hold on;
    %         plot(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'Color',[0.9290 0.6940 0.1250]);hold on;
%     area(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak), ...
%         'EdgeColor',[0.3010 0.7450 0.9330],'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.5,'LineWidth',2);
%     hold on
%     area(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end), ...
%         'EdgeColor',[0.3010 0.7450 0.9330],'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.5,'LineWidth',2);
%     %     plot(X,N,'k');xlim([-0.1,0.1])
%     hold on
%     xline(MedianCor,'Color',[0.3010 0.7450 0.9330],'LineStyle','--');
%     xlim([-0.1,0.1])
%     xlabel("Correlation Coefficient",'FontSize',18); ylabel("Normalised density",'FontSize',18);
%     title("PC")
%     hold off
    
%     subplot(3,2,4)
    Cor_now = reshape(high_high_SIC_cr_res.',1,[]); % Red
    [MeanCor,MedianCor,FWHM,IDpeak,xLeft_sym,Leftfit_sym,xRight_sym,Rightfit_sym] = FWHM_compute(Cor_now);
    MeanCor_all(7) = MeanCor;
    MedianCor_all(7) = MedianCor;
    FWHM_all(7) = FWHM;
    %         plot(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'Color',[0.4660 0.6740 0.1880]);hold on;
    %         plot(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'Color',[0.4660 0.6740 0.1880]);hold on;
%     area(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'EdgeColor',[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840], ...
%         'FaceAlpha',0.5,'LineWidth',2);
%     hold on
%     area(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'EdgeColor',[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840], ...
%         'FaceAlpha',0.5,'LineWidth',2);
%     hold on
%     xline(MedianCor,'Color',[0.6350 0.0780 0.1840],'LineStyle','--');
%     %     plot(X,N,'k');xlim([-0.1,0.1])
%     hold on
    Cor_now = reshape(high_low_SIC_cr_res.',1,[]); % Blue
    [MeanCor,MedianCor,FWHM,IDpeak,xLeft_sym,Leftfit_sym,xRight_sym,Rightfit_sym] = FWHM_compute(Cor_now); % Blue
    MeanCor_all(8) = MeanCor;
    MedianCor_all(8) = MedianCor;
    FWHM_all(8) = FWHM;
    %         plot(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'Color',[0.9290 0.6940 0.1250]);hold on;
    %         plot(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'Color',[0.9290 0.6940 0.1250]);hold on;
%     area(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak), ...
%         'EdgeColor',[0.3010 0.7450 0.9330],'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.5,'LineWidth',2);
%     hold on
%     area(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end), ...
%         'EdgeColor',[0.3010 0.7450 0.9330],'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.5,'LineWidth',2);
%     %     plot(X,N,'k');xlim([-0.1,0.1])
%     hold on
%     xline(MedianCor,'Color',[0.3010 0.7450 0.9330],'LineStyle','--');
%     xlim([-0.1,0.1])
%     xlabel("Correlation Coefficient",'FontSize',18); ylabel("Normalised density",'FontSize',18);
%     title("High SIC")
%     hold off

%     subplot(3,2,5)
    Cor_now = reshape(nonpc_nonpc_cr_res.',1,[]); % Red
    [MeanCor,MedianCor,FWHM,IDpeak,xLeft_sym,Leftfit_sym,xRight_sym,Rightfit_sym] = FWHM_compute(Cor_now);
    MeanCor_all(9) = MeanCor;
    MedianCor_all(9) = MedianCor;
    FWHM_all(9) = FWHM;
    %         plot(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'Color',[0.4660 0.6740 0.1880]);hold on;
    %         plot(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'Color',[0.4660 0.6740 0.1880]);hold on;
%     area(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'EdgeColor',[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840], ...
%         'FaceAlpha',0.5,'LineWidth',2);
%     hold on
%     area(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'EdgeColor',[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840], ...
%         'FaceAlpha',0.5,'LineWidth',2);
%     hold on
%     xline(MedianCor,'Color',[0.6350 0.0780 0.1840],'LineStyle','--');
%     %     plot(X,N,'k');xlim([-0.1,0.1])
%     hold on
    Cor_now = reshape(nonpc_pc_cr_res.',1,[]); % Blue
    [MeanCor,MedianCor,FWHM,IDpeak,xLeft_sym,Leftfit_sym,xRight_sym,Rightfit_sym] = FWHM_compute(Cor_now); % Blue
    MeanCor_all(10) = MeanCor;
    MedianCor_all(10) = MedianCor;
    FWHM_all(10) = FWHM;
    %         plot(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'Color',[0.9290 0.6940 0.1250]);hold on;
    %         plot(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'Color',[0.9290 0.6940 0.1250]);hold on;
%     area(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak), ...
%         'EdgeColor',[0.3010 0.7450 0.9330],'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.5,'LineWidth',2);
%     hold on
%     area(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end), ...
%         'EdgeColor',[0.3010 0.7450 0.9330],'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.5,'LineWidth',2);
%     %     plot(X,N,'k');xlim([-0.1,0.1])
%     hold on
%     xline(MedianCor,'Color',[0.3010 0.7450 0.9330],'LineStyle','--');
%     xlim([-0.1,0.1])
%     xlabel("Correlation Coefficient",'FontSize',18); ylabel("Normalised density",'FontSize',18);
%     title("Non PC")
%     hold off
    
%     subplot(3,2,6)
    Cor_now = reshape(low_low_SIC_cr_res.',1,[]); % Red
    [MeanCor,MedianCor,FWHM,IDpeak,xLeft_sym,Leftfit_sym,xRight_sym,Rightfit_sym] = FWHM_compute(Cor_now);
    MeanCor_all(11) = MeanCor;
    MedianCor_all(11) = MedianCor;
    FWHM_all(11) = FWHM;
    %         plot(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'Color',[0.4660 0.6740 0.1880]);hold on;
    %         plot(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'Color',[0.4660 0.6740 0.1880]);hold on;
%     area(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'EdgeColor',[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840], ...
%         'FaceAlpha',0.5,'LineWidth',2);
%     hold on
%     area(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'EdgeColor',[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840], ...
%         'FaceAlpha',0.5,'LineWidth',2);
%     hold on
%     xline(MedianCor,'Color',[0.6350 0.0780 0.1840],'LineStyle','--');
%     %     plot(X,N,'k');xlim([-0.1,0.1])
%     hold on
    Cor_now = reshape(low_high_SIC_cr_res.',1,[]); % Blue
    [MeanCor,MedianCor,FWHM,IDpeak,xLeft_sym,Leftfit_sym,xRight_sym,Rightfit_sym] = FWHM_compute(Cor_now); % Blue
    MeanCor_all(12) = MeanCor;
    MedianCor_all(12) = MedianCor;
    FWHM_all(12) = FWHM;
    %         plot(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'Color',[0.9290 0.6940 0.1250]);hold on;
    %         plot(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'Color',[0.9290 0.6940 0.1250]);hold on;
%     area(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak), ...
%         'EdgeColor',[0.3010 0.7450 0.9330],'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.5,'LineWidth',2);
%     hold on
%     area(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end), ...
%         'EdgeColor',[0.3010 0.7450 0.9330],'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.5,'LineWidth',2);
%     %     plot(X,N,'k');xlim([-0.1,0.1])
%     hold on
%     xline(MedianCor,'Color',[0.3010 0.7450 0.9330],'LineStyle','--');
%     xlim([-0.1,0.1])
%     xlabel("Correlation Coefficient",'FontSize',18); ylabel("Normalised density",'FontSize',18);
%     title("Low SIC")
%     hold off
    end
    
    data.Cor_Mspks.MeanCor_all = MeanCor_all;
    data.Cor_Mspks.MedianCor_all = MedianCor_all;
    data.Cor_Mspks.FWHM_all = FWHM_all;
    
    for loop_pc_non_pc = 1:5
        
        switch loop_pc_non_pc
            case 1 % All cells
                cr_res_temp = cr_res;
                cr_run_temp = cr_run;
            case 2 % PC-PC
                cr_res_temp = cr_res(pc_list,pc_list);
                cr_run_temp = cr_run(pc_list,pc_list);
            case 3 % non PC-non PC
                cr_res_temp = cr_res(non_pc_list,non_pc_list);
                cr_run_temp = cr_run(non_pc_list,non_pc_list);
            case 4 % PC-non PC
                cr_res_temp = cr_res(pc_list,non_pc_list);
                cr_run_temp = cr_run(pc_list,non_pc_list);
            case 5 % non PC-PC
                cr_res_temp = cr_res(non_pc_list,pc_list);
                cr_run_temp = cr_run(non_pc_list,pc_list);
        end
        
        k_values = zeros(2,6);
        % figure;
        for loop_links = 1:6
%             subplot(3,2,loop_links)
            if loop_links == 1 | loop_links == 3 | loop_links == 5
                network = cr_res_temp;
                % link_dist = sortrows([[1:size(network,1)]' sum(network)'],2);
                proportion_pos_neg = zeros(size(network,1),6);
                for i = 1:size(network,1)
                    proportion_pos_neg(i,1) = i;
                    proportion_pos_neg(i,2) = sum(network(i,:));% Sum of all links
                    proportion_pos_neg(i,3) = sum(network(i,:) > 0) / size(network,2); % Proportion of positive links
                    proportion_pos_neg(i,4) = sum(network(i,network(i,:) > 0)); % Sum of positive links
                    proportion_pos_neg(i,5) = sum(network(i,:) < 0) / size(network,2); % Proportion of negative links
                    proportion_pos_neg(i,6) = sum(network(i,network(i,:) < 0)); % Sum of positive links
                end   
                rest_proportion_pos_neg = proportion_pos_neg;
                
            elseif loop_links == 2 | loop_links == 4 | loop_links == 6
                network = cr_run_temp;
                % link_dist = sortrows([[1:size(network,1)]' sum(network)'],2);
                proportion_pos_neg = zeros(size(network,1),6);
                for i = 1:size(network,1)
                    proportion_pos_neg(i,1) = i;
                    proportion_pos_neg(i,2) = sum(network(i,:));% Sum of all links
                    proportion_pos_neg(i,3) = sum(network(i,:) > 0) / size(network,2); % Proportion of positive links
                    proportion_pos_neg(i,4) = sum(network(i,network(i,:) > 0)); % Sum of positive links
                    proportion_pos_neg(i,5) = sum(network(i,:) < 0) / size(network,2); % Proportion of negative links
                    proportion_pos_neg(i,6) = sum(network(i,network(i,:) < 0)); % Sum of positive links
                end
                run_proportion_pos_neg = proportion_pos_neg;
            end
                        
            if loop_links == 1 | loop_links == 2 % All
                link_dist = sortrows(proportion_pos_neg(:,[1 2]),2);
            elseif loop_links == 3 | loop_links == 4 % Positive
                link_dist = sortrows(proportion_pos_neg(:,[1 4]),2);
            elseif loop_links == 5 | loop_links == 6 % Negative
                link_dist = sortrows(proportion_pos_neg(:,[1 6]),2);
                link_dist(:,2) = link_dist(:,2) * -1; % If negative
            end
            [N,edges] = histcounts(link_dist(:,2),0:0.2:max(link_dist(:,2))+0.2,'Normalization','probability');
            k_mean = mean(link_dist(:,2),'omitnan');
            k_median = median(link_dist(:,2),'omitnan');
            k_values(1,loop_links) = k_mean;
            k_values(2,loop_links) = k_median;
            % [N,edges] = histcounts(link_dist(:,2),0:5:max(link_dist(:,2)+5),'Normalization','probability'); %1:max(link_dist(:,2))+1
            degree_dist = [[0:0.2:max(link_dist(:,2))]' N'];
            
%             bar(edges(1:end-1),degree_dist(:,2));
%             xline(k_mean,'k--');
%             xline(k_median,'r--');
            %     loglog(degree_dist(:,2),'x');
%             if loop_links == 1
%                 title("Resting (All links)")
%             elseif loop_links == 2
%                 title("Running (All links)")
%             elseif loop_links == 3
%                 title("Resting (Positive links)")
%             elseif loop_links == 4
%                 title("Running (Positive links)")
%             elseif loop_links == 5
%                 title("Resting (Negative links)")
%             elseif loop_links == 6
%                 title("Running (Negative links)")
%             end
        end
        
        switch loop_pc_non_pc
            case 1 % All cells
                data.Cor_Mspks.links.all.rest_proportion_pos_neg = rest_proportion_pos_neg;
                data.Cor_Mspks.links.all.run_proportion_pos_neg = run_proportion_pos_neg;
                data.Cor_Mspks.links.all.k_values = k_values;
            case 2 % PC-PC
                data.Cor_Mspks.links.pc_pc.rest_proportion_pos_neg = rest_proportion_pos_neg;
                data.Cor_Mspks.links.pc_pc.run_proportion_pos_neg = run_proportion_pos_neg;
                data.Cor_Mspks.links.pc_pc.k_values = k_values;
            case 3 % non PC-non PC
                data.Cor_Mspks.links.nonpc_nonpc.rest_proportion_pos_neg = rest_proportion_pos_neg;
                data.Cor_Mspks.links.nonpc_nonpc.run_proportion_pos_neg = run_proportion_pos_neg;
                data.Cor_Mspks.links.nonpc_nonpc.k_values = k_values;
            case 4 % PC-non PC
                data.Cor_Mspks.links.pc_nonpc.rest_proportion_pos_neg = rest_proportion_pos_neg;
                data.Cor_Mspks.links.pc_nonpc.run_proportion_pos_neg = run_proportion_pos_neg;
                data.Cor_Mspks.links.pc_nonpc.k_values = k_values;
            case 5 % non PC-PC
                data.Cor_Mspks.links.nonpc_pc.rest_proportion_pos_neg = rest_proportion_pos_neg;
                data.Cor_Mspks.links.nonpc_pc.run_proportion_pos_neg = run_proportion_pos_neg;
                data.Cor_Mspks.links.nonpc_pc.k_values = k_values;
        end
        
    end
    
    
    %%
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
