function [obj, varargout] = mtcellcombined(varargin)
%@mtcellcombined Constructor function for mtcellcombined class
%   OBJ = mtcellcombined(varargin)
%
%   OBJ = mtcellcombined('auto') attempts to create a mtcellcombined object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on mtcellcombined %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%example [as, Args] = mtcellcombined('save','redo')
%
%dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
				'ObjectLevel','Cell','RequiredFile','spiketimes.mat', ...
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
Args.classname = 'mtcellcombined';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'mtcellcombined';

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
    spiketimes_all = load(Args.RequiredFile);
    spiketimes_all = spiketimes_all.spiketimes;
    
    for cell_idx = 1:sessionData.nNeuron
        
        if mod(cell_idx,50) == 0
            fprintf('Processing cell %i\n',cell_idx);
        end
        
        spiketimes = spiketimes_all.("n" + sprintf('%04d',cell_idx));
        
        if isempty(spiketimes)
            spiketimes = 0;
        end
        
        for repeat = 1:3 % 1 = full trial, 2 = 1st half, 3 = 2nd half
            
            % selecting rows from sessionTimeC
            stc = sessionData.session_data_exclude_zero_trials(:,[1,4,3]);
            stc(:,1) = stc(:,1);
            stc(:,4) = [diff(stc(:,1)); 0];
            
            halving_markers = ones(length(stc),1);
            % halving_markers(ceil(length(stc)/2):end) = 2;
            halving_markers(sessionData.sessionMidpoint:end) = 2;
            
            conditions0 = ones(size(stc,1),1);
            
            if repeat == 1
                conditions = conditions0;
                NumShuffles = Args.NumShuffles;
            elseif repeat == 2
                conditions = conditions0 & (halving_markers==1);
                NumShuffles = 0;
            elseif repeat == 3
                conditions = conditions0 & (halving_markers==2);
                NumShuffles = 0;
            end
            
            clear halving_markers
            
            % spike shuffling
            
            spiketimes_ = spiketimes';
            maxTime = spiketimes_(end);
            tShifts = [0 ((rand([1,NumShuffles])*diff(Args.ShuffleLimits))+Args.ShuffleLimits(1))*maxTime];
            full_arr = repmat(spiketimes_, NumShuffles+1, 1);
            full_arr = full_arr + tShifts';
            keepers = length(spiketimes_) - sum(full_arr>maxTime, 2);
            for row = 2:size(full_arr,1)
                full_arr(row,:) = [full_arr(row,1+keepers(row):end)-maxTime full_arr(row,1:keepers(row))];
            end
            flat_spiketimes = NaN(2,size(full_arr,1)*size(full_arr,2));
            temp = full_arr';
            flat_spiketimes(1,:) = temp(:);
            flat_spiketimes(2,:) = repelem(1:size(full_arr,1), size(full_arr,2));
            flat_spiketimes = flat_spiketimes';
            flat_spiketimes = sortrows(flat_spiketimes);
            
            flat_spiketimes(flat_spiketimes(:,1) < sessionData.actual_start_time,:) = [];
            
            clear full_arr temp
            
            if Args.ThresVel > 0
                conditions = conditions & (sessionData.velocity_averaged(:,4) > Args.ThresVel);
            end
            
            bins = 1:Args.BinSize;
            
            consol_arr = zeros(Args.BinSize, NumShuffles + 1);
            if repeat == 1
                binArr = zeros(Args.BinSize * sessionData.nTrials, NumShuffles + 1);
            end
            
            consol_arr = zeros(Args.BinSize, NumShuffles+1);
            binArr = zeros(sessionData.nTrials * Args.BinSize, NumShuffles+1);
            
            % C++ code for spike binning
%             A = coder.typeof(0, [inf,2], [1 0]);
%             B = coder.typeof(0, [inf,4], [1 0]);
%             C = coder.typeof(0, [inf,1], [1 0]);
%             D = coder.typeof(0, [100,10001], [0 1]);
%             E = coder.typeof(0, [inf,10001], [1 1]);
            %codegen spikeBinner_matlab -args {A, B, C, 1, D, E} spikeBinner.c spikeBinner.h -report
            
            [consol_arr, binArr] = spikeBinner_matlab_mex(flat_spiketimes, stc, double(conditions), repeat, consol_arr, binArr);
            
            %         interval = 1;
            %
            %         for sp = 1:size(flat_spiketimes,1)
            %
            %             while interval < size(stc,1)
            %                 if flat_spiketimes(sp,1) >= stc(interval) && flat_spiketimes(sp,1) < stc(interval+1)
            %                     break;
            %                 end
            %                 interval = interval + 1;
            %             end
            %
            %             bins_hit = stc(interval,3);
            %             bins_hit = bins_hit(logical(conditions(interval)));
            %
            %             bins_hit(bins_hit<1) = [];
            %
            %             consol_arr(bins_hit,flat_spiketimes(sp,2)) = consol_arr(bins_hit,flat_spiketimes(sp,2)) + 1;
            %
            %             % For individual bins per trial
            %             if repeat == 1
            %                 tbins_hit = stc(interval,2);
            %                 tbins_hit = tbins_hit(logical(conditions(interval)));
            %                 tbins_hit(tbins_hit<1) = [];
            %
            %                 binArr(tbins_hit,flat_spiketimes(sp,2)) = binArr(tbins_hit,flat_spiketimes(sp,2)) + 1;
            %             end
            %         end
            
            firing_counts_full = consol_arr';
            stc_ss = stc(find(conditions==1),[2:4]);
            stc_ss(~(stc_ss(:,2) > 0), :) = [];
            % stc_ss = [stc_ss; [1600 0]];
            gpdur = accumarray(stc_ss(:,2),stc_ss(:,3))';
            if repeat == 1
                binSpikeCount = reshape(binArr(:,1), [Args.BinSize,sessionData.nTrials])';
                binDuration = reshape(accumarray(stc_ss(:,1),stc_ss(:,3))', [Args.BinSize,sessionData.nTrials])';
            end
            clear conditions0 conditions binArr flat_spiketimes stc stc_ss
            
            if Args.AdaptiveSmooth
                
                firing_rates_full_raw = firing_counts_full./repmat(gpdur,size(firing_counts_full,1),1);
                to_save = NaN(1,Args.BinSize);
                to_save(bins) = firing_rates_full_raw(1,bins);
                
                if repeat == 1
                    data.("n" + sprintf('%04d',cell_idx)).maps_raw = to_save;
                elseif repeat == 2
                    data.("n" + sprintf('%04d',cell_idx)).maps_raw1 = to_save;
                elseif repeat == 3
                    data.("n" + sprintf('%04d',cell_idx)).maps_raw2 = to_save;
                end
                nan_track = isnan(to_save);
                
                % smoothing part starts here
                
                wip = ones(NumShuffles+1,1);
                
                gpdur1 = zeros(1,Args.BinSize);
                gpdur1(bins) = gpdur(bins);
                
                gpdur1 = repmat(gpdur1,NumShuffles + 1,1);
                gpdur1 = permute(gpdur1,[2,1]);
                firing_counts_full1 = zeros(NumShuffles + 1, Args.BinSize);
                firing_counts_full1(:,bins) = firing_counts_full(:,bins);
                firing_counts_full1 = permute(firing_counts_full1,[2,1]);
                
                to_compute = 1:0.5:Args.BinSize/2;
                possible = NaN(length(to_compute),2,Args.BinSize,NumShuffles + 1);
                to_fill = NaN(size(possible,3), size(possible,4));
                to_fill_time = NaN(size(possible,3), size(possible,4));
                to_fill_radius = NaN(size(possible,3), size(possible,4));
                
                for idx = 1:length(to_compute)
                    
                    %                     f=fspecial('disk',to_compute(idx));
                    %                     f(f>=(max(max(f))/3))=1;
                    %                     f(f~=1)=0;
                    
                    f = ones(1, to_compute(idx)*2+1)';
                    % f = ones(1, to_compute(idx)*2+1) ./ (to_compute(idx)*2+1);
                    
                    possible(idx,1,:,:) = repmat(imfilter(gpdur1(:,1), f, 'conv'), 1, NumShuffles+1);   %./scaler;
                    possible(idx,2,:,find(wip)) = imfilter(firing_counts_full1(:,find(wip)), f, 'conv');   %./scaler;
                    
                    logic1 = squeeze(Args.Alpha./(possible(idx,1,:,:).*sqrt(possible(idx,2,:,:))) <= to_compute(idx));
                    slice1 = squeeze(possible(idx,1,:,:));
                    slice2 = squeeze(possible(idx,2,:,:));
                    
                    to_fill(logic1 & isnan(to_fill)) = slice2(logic1 & isnan(to_fill))./slice1(logic1 & isnan(to_fill));
                    to_fill_time(logic1 & isnan(to_fill_time)) = slice1(logic1 & isnan(to_fill_time));
                    to_fill_radius(logic1 & isnan(to_fill_radius)) = to_compute(idx);
                    
                    %                     disp('smoothed with kernel size:');
                    %                     disp(to_compute(idx));
                    %                     disp('grids left');
                    %                     disp(sum(sum(isnan(to_fill(:,:)))));
                    
                    check = squeeze(sum(isnan(to_fill),1));
                    wip(check==0) = 0;
                    
                    %                     if sum(sum(isnan(to_fill(:,:)))) == 0
                    %                         disp('breaking');
                    %                         break;
                    %                     end
                end
                
                clear possible
                
                to_fill(isnan(to_fill)) = 0;
                to_fill = permute(to_fill, [2 1]);
                to_fill = to_fill(:,bins);
                
                to_fill_time(isnan(to_fill_time)) = 0;
                to_fill_time = permute(to_fill_time, [2 1]);
                to_fill_time = to_fill_time(:,bins);
                
                to_fill_radius = permute(to_fill_radius, [2 1]);
                
                % firing_rates_full = to_fill;
                firing_rates_full = to_fill ./ to_fill_radius;
                
                % smoothing part ends
                
                if repeat == 1
                    to_save = NaN(1,Args.BinSize);
                    to_save(bins) = firing_rates_full(1,:);
                    to_save(find(nan_track==1)) = nan;
                    data.("n" + sprintf('%04d',cell_idx)).maps_adsm = to_save;
                    to_save = NaN(size(firing_rates_full,1)-1,Args.BinSize);
                    to_save(:,bins) = firing_rates_full(2:end,:);
                    data.("n" + sprintf('%04d',cell_idx)).maps_adsmsh = to_save;
                    to_save = NaN(size(firing_rates_full,1)-1,Args.BinSize);
                    to_save(:,bins) = to_fill_time(2:end,:);
                    data.("n" + sprintf('%04d',cell_idx)).dur_adsmsh = to_save;
                    to_save = NaN(1,Args.BinSize);
                    to_save(bins) = to_fill_time(1,:);
                    data.("n" + sprintf('%04d',cell_idx)).dur_adsm = to_save;
                    data.("n" + sprintf('%04d',cell_idx)).radii = to_fill_radius(1,:);
                    data.("n" + sprintf('%04d',cell_idx)).radiish = to_fill_radius(2:end,:);
                    data.("n" + sprintf('%04d',cell_idx)).binSpikeCount = binSpikeCount;
                    data.("n" + sprintf('%04d',cell_idx)).binDuration = binDuration;
                    data.("n" + sprintf('%04d',cell_idx)).binFiringRate = binSpikeCount./binDuration;
                    
                elseif repeat == 2
                    to_save = NaN(1,Args.BinSize);
                    to_save(bins) = firing_rates_full(1,:);
                    to_save(find(nan_track==1)) = nan;
                    data.("n" + sprintf('%04d',cell_idx)).maps_adsm1 = to_save;
                    %to_save = NaN(size(firing_rates_full,1)-1,Args.BinSize);
                    %to_save(:,bins) = firing_rates_full(2:end,:);
                    %data.("n" + sprintf('%04d',cell_idx)).maps_adsmsh1 = to_save;
                    %to_save = NaN(size(firing_rates_full,1)-1,Args.BinSize);
                    %to_save(:,bins) = to_fill_time(2:end,:);
                    %data.("n" + sprintf('%04d',cell_idx)).dur_adsmsh1 = to_save;
                    to_save = NaN(1,Args.BinSize);
                    to_save(bins) = to_fill_time(1,:);
                    data.("n" + sprintf('%04d',cell_idx)).dur_adsm1 = to_save;
                    data.("n" + sprintf('%04d',cell_idx)).radii1 = to_fill_radius;
                elseif repeat == 3
                    to_save = NaN(1,Args.BinSize);
                    to_save(bins) = firing_rates_full(1,:);
                    to_save(find(nan_track==1)) = nan;
                    data.("n" + sprintf('%04d',cell_idx)).maps_adsm2 = to_save;
                    %to_save = NaN(size(firing_rates_full,1)-1,Args.BinSize);
                    %to_save(:,bins) = firing_rates_full(2:end,:);
                    %data.("n" + sprintf('%04d',cell_idx)).maps_adsmsh2 = to_save;
                    %to_save = NaN(size(firing_rates_full,1)-1,Args.BinSize);
                    %to_save(:,bins) = to_fill_time(2:end,:);
                    %data.("n" + sprintf('%04d',cell_idx)).dur_adsmsh2 = to_save;
                    to_save = NaN(1,Args.BinSize);
                    to_save(bins) = to_fill_time(1,:);
                    data.("n" + sprintf('%04d',cell_idx)).dur_adsm2 = to_save;
                    data.("n" + sprintf('%04d',cell_idx)).radii2 = to_fill_radius;
                end
                
            else
                firing_rates_full = firing_counts_full./repmat(gpdur,size(firing_counts_full,1),1);
                to_fill_time = repmat(gpdur,size(firing_counts_full,1),1); % HM added
                
                to_save = NaN(1,Args.BinSize);
                to_save(bins) = firing_rates_full(1,:);
                if repeat == 1
                    data.("n" + sprintf('%04d',cell_idx)).maps_raw = to_save;
                elseif repeat == 2
                    data.("n" + sprintf('%04d',cell_idx)).maps_raw1 = to_save;
                elseif repeat == 3
                    data.("n" + sprintf('%04d',cell_idx)).maps_raw2 = to_save;
                end
            end
            
            gpdur1 = zeros(NumShuffles+1,Args.BinSize);
            gpdur1(:,bins) = to_fill_time;
            Pi1 = gpdur1./sum(gpdur1,2); % consider nansum to play safe
            %         Pi1 = repmat(Pi1, NumShuffles+1, 1);
            
            lambda_i = NaN(NumShuffles+1,Args.BinSize);
            lambda_i(:,bins) = firing_rates_full;
            lambda_i(isnan(lambda_i)) = 0;
            lambda_bar = sum(Pi1 .* lambda_i,2);
            % divide firing for each position by the overall mean
            FRratio = lambda_i./repmat(lambda_bar,1,Args.BinSize);
            % compute first term in SIC
            SIC1 = Pi1 .* lambda_i;
            SIC2 = log2(FRratio);
            zeros_placing = SIC1==0;
            
            bits_per_sec = SIC1 .* SIC2 ./ lambda_bar;
            bits_per_sec(zeros_placing) = NaN;
            lambda_bar_ok = lambda_bar>0;
            lambda_bar_bad = ~lambda_bar_ok;
            sic_out = nansum(bits_per_sec, 2);
            sic_out(lambda_bar_bad) = NaN;
            
            %     histogram(sic_out);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if repeat == 1
                data.("n" + sprintf('%04d',cell_idx)).SIC = sic_out(1);
                % shuffle_SIC(shuffle_no) = sic_out(1);
                data.("n" + sprintf('%04d',cell_idx)).SICsh = sic_out(2:end,1);
                %     data.("n" + sprintf('%04d',cell_idx)).median_occ_firings = median_stats';
                %     data.("n" + sprintf('%04d',cell_idx)).variance_occ_firings = var_stats';
                %     data.("n" + sprintf('%04d',cell_idx)).perc_occ_firings = perc_stats';
                %     data.("n" + sprintf('%04d',cell_idx)).occ_data = occ_data;
            elseif repeat == 2
                data.("n" + sprintf('%04d',cell_idx)).SIC1 = sic_out;
            elseif repeat == 3
                data.("n" + sprintf('%04d',cell_idx)).SIC2 = sic_out;
            end
            %     data.("n" + sprintf('%04d',cell_idx)).median_occ_firings = median_stats';
            %     data.("n" + sprintf('%04d',cell_idx)).variance_occ_firings = var_stats';
            %     data.("n" + sprintf('%04d',cell_idx)).perc_occ_firings = perc_stats';
            
            %     data.("n" + sprintf('%04d',cell_idx)).occ_data = occ_data;
        end
        
    end
    
    data.nTrials = sessionData.nTrials;
    data.nNeuron = sessionData.nNeuron;
    
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
