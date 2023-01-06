function [obj, varargout] = mtplacefield(varargin)
%@mtplacefield Constructor function for mtplacefield class
%   OBJ = mtplacefield(varargin)
%
%   OBJ = mtplacefield('auto') attempts to create a mtplacefield object by ...
%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on mtplacefield %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%example [as, Args] = mtplacefield('save','redo')
%
%dependencies:

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
    'ObjectLevel','Cell','RequiredFile','mtcell.mat', ...
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
Args.classname = 'mtplacefield';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'mtplacefield';

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
    %     cd ../..
    %     sessionData = mtsess('Auto', varargin{:});
    %     sessionData = sessionData.data;
    %     cd(ori);
    cellData = load(Args.RequiredFile);
    cellData = cellData.mtcell.data;
    
    % Determine place fields
    % basemapLsm = imgaussfilt(data.maps_adsm, 3);
    basemapLsm = cellData.maps_adsm;
    basemapLrw = cellData.maps_raw;
    
    Fs = 100;            % Sampling frequency
    T = 1/Fs;             % Sampling period
    L = 100;             % Length of signal
    s = (0:L-1)*T;        % Space vector
    
    Y = fft(basemapLsm);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
%         figure;
%         plot(f,P1);
%         title('Single-Sided Amplitude Spectrum of X(t)');
%         xlabel('f (Hz)');
%         ylabel('|P1(f)|');
%         xline(10,'--r');
    
    % figure;
    filt_freq_arr = 0:2:22;
    % Low-Pass Filter
    for filt_freq_idx = 1:12
        Y_filt = Y;
        Y_filt(find(f > filt_freq_arr(filt_freq_idx),1):end-find(f > filt_freq_arr(filt_freq_idx),1)) = 0; % Low-Pass filter at 150Hz
        mapLsm_FFT_Low_Pass = abs(ifft(Y_filt));
        
%         subplot(4,3,filt_freq_idx);
%         plot(basemapLrw, 'g');legends{1} = sprintf('Raw');
%         hold on
%         plot(basemapLsm, 'b');legends{2} = sprintf('AdSmooth');
%         hold on
%         plot(mapLsm_FFT_Low_Pass, 'r');legends{3} = sprintf('AdSm-FFT-Low-Pass');
%         %legend(legends);
%         title("filt_freq =" + filt_freq_arr(filt_freq_idx) + "Hz")
    end
    
    filt_freq = f(find((cumsum(P1) ./ sum(P1) > Args.FiltFrac), 1, 'first')); % 75% threshold
    Y_filt = Y;
    Y_filt(find(f > filt_freq,1):end-find(f > filt_freq,1)) = 0; % Low-Pass filter at 15Hz
    mapLsm_FFT_Low_Pass = abs(ifft(Y_filt));
    
%         figure;
%         plot(basemapLrw, 'g');legends{1} = sprintf('Raw');
%         hold on
%         plot(basemapLsm, 'b');legends{2} = sprintf('AdSmooth');
% %         title('Firing Rate Map'); xlabel('Position Bins'); ylabel('Firing Rate (spike/sec)');
%         hold on
%         plot(mapLsm_FFT_Low_Pass, 'r');legends{3} = sprintf('AdSm-FFT-Low-Pass');
%         legend(legends);
%         title("Firing Rate Map; " + "filt_freq =" + filt_freq + "Hz"); xlabel('Position Bins'); ylabel('Firing Rate (spike/sec)');
%     
    % GMM Iterative Method
    
    % Interpolate
    mapLsm_FFT_Low_Pass_padded = [mapLsm_FFT_Low_Pass(76:100) mapLsm_FFT_Low_Pass mapLsm_FFT_Low_Pass(1:25)];
    
    fit_data = [-24:1:125; mapLsm_FFT_Low_Pass_padded]';
    fit_data_valid = fit_data(26:125,:);
    
    guess = []; % Mean (Centre Bin), Std (Half-Width), Height (Rate)
    fit_data_ = fit_data;
    fit_data_valid_ = fit_data_valid;
    min_peak_height = mean(fit_data_valid_(:,2), 'omitnan');
    max_gauss_width = 25;
    gauss_overlap_threshold = 1;
    mean_threshold = mean(fit_data_valid_(:,2), 'omitnan');
    std_threshold = std(fit_data_valid_(:,2), 'omitnan');
    height_threshold = mean_threshold + std_threshold;
    cv_threshold = std_threshold / mean_threshold;
    vmr_threshold = (std_threshold^2) / mean_threshold;
    
    data.mean_threshold = mean_threshold;
    data.std_threshold = std_threshold;
    data.cv_threshold = cv_threshold;
    data.vmr_threshold = vmr_threshold;
    
    % Finding baseline
%     while size(guess,1) < 3
%         % Find max peak
%         [max_height, max_bin_idx] = max(fit_data_valid_(:,2), [], 'omitnan');
%         max_bin = fit_data_valid_(max_bin_idx);
%         % Stop searching for peaks once height drops below height threshold
%         if max_height <= height_threshold
%             break
%         end
%         
%         % Set the guess parameters for gaussian fitting
%         guess_bin = max_bin; % Mean
%         guess_height = max_height; % Height
%         
%         % Find half height index on each side of the center bin (gaussian mean)
%         half_height = 0.5 * guess_height;
%         max_bin_idx_ = find(fit_data_(:,1) == guess_bin);
%         left_idx = find(fit_data_(1:max_bin_idx_,2) <= half_height, 1, 'last');
%         left_side = fit_data_(left_idx,1);
%         right_idx = find(fit_data_(max_bin_idx_:end,2) <= half_height, 1, 'first') + max_bin_idx_ - 1;
%         right_side = fit_data_(right_idx,1);
%         
%         % Use shorter side
%         if isempty(left_side) && ~isempty(right_side)
%             shorter_width = abs(right_side - max_bin);
%         elseif isempty(right_side) && ~isempty(left_side)
%             shorter_width = abs(left_side - max_bin);
%         elseif isempty(left_side) && isempty(right_side)
%             break
%         else
%             shorter_width = min(abs(left_side - max_bin), abs(right_side - max_bin));
%         end
%         fwhm = shorter_width * 2;
%         guess_std = fwhm / (2 * sqrt(2 * log(2)));
%         
%         % Halt fitting process if peak drops below minimum height
%         if guess_height < min_peak_height
%             break
%         end
%         
%         if guess_std > max_gauss_width
%             break
%         end
%         
%         guess = [guess; [guess_bin, guess_std, guess_height]];
%         
%         ys = zeros(size(fit_data,1),1);
%         ys = ys + guess_height * exp(-(fit_data(:,1)-guess_bin).^2 / (2*(guess_std.^2)));
%         
%         left = floor(guess_bin - guess_std);
%         right = ceil(guess_bin + guess_std);
%         ys(1:left+25-1) = 0;
%         ys(right+25+1:end) = 0;
%         ys(left+25:right+25) = ys(left+25:right+25) - (guess_height - min([ys(left+25) ys(right+25)]));
%         
%         if guess_bin < 50
%             ys(101:125) = ys(1:25);
%         else
%             ys(26:50) = ys(126:150);
%         end
%         
%         fit_data_(:,2) = fit_data_(:,2) - ys;
%         fit_data_valid_ = fit_data_(26:125,:);
%         
% %         if mean(fit_data_valid_(:,2), 'omitnan') < mean(fit_data_valid(:,2), 'omitnan')
% %             break
% %         end
%         
%     end
%     
%     % Check peak heights compared to baseline
%     %std_threshold = Args.PeakThreshold * std(fit_data_valid_(:,2), 'omitnan');
%     baseline_mean_height = mean(fit_data_valid_(:,2), 'omitnan');
%     peak_height_from_baseline = max(fit_data_valid(:,2), [], 'omitnan') - baseline_mean_height;
%     
%     %data.std_threshold = std_threshold;
%     data.baseline_mean_height = baseline_mean_height;
%     data.peak_height_from_baseline = peak_height_from_baseline;
    
%     figure;
%     plot(fit_data_valid(:,1), fit_data_valid(:,2), 'b'); legends{1} = sprintf('Firing Rate Map');
%     hold on
%     plot(fit_data_(26:125, 1), fit_data_(26:125, 2),'r'); legends{2} = sprintf('Estimated Baseline');
%     legend(legends);
%     hold off
    
    guess = []; % Mean (Centre Bin), Std (Half-Width), Height (Rate)
    fit_data_ = fit_data;
    fit_data_valid_ = fit_data_valid;
    
    % . Find peak: Loop through, finding a candidate peak, and fitting with a guess gaussian
    %   Stopping procedures: limit on # of peaks, or relative or absolute height thresholds
    while size(guess,1) < 10
        
        % Find max peak
        [max_height, max_bin_idx] = max(fit_data_valid_(:,2), [], 'omitnan');
        max_bin = fit_data_valid_(max_bin_idx);
        
        % Set the guess parameters for gaussian fitting
        guess_bin = max_bin; % Mean
        guess_height = max_height; % Height
                
        % Find half height index on each side of the center bin (gaussian mean)
        half_height = 0.5 * max_height;
        max_bin_idx_ = find(fit_data_(:,1) == max_bin);
        left_idx = find(fit_data_(1:max_bin_idx_,2) <= half_height, 1, 'last');
        left_side = fit_data_(left_idx,1);
        right_idx = find(fit_data_(max_bin_idx_:end,2) <= half_height, 1, 'first') + max_bin_idx_ - 1;
        right_side = fit_data_(right_idx,1);
        
        % Use shorter side
        if isempty(left_side) && ~isempty(right_side)
            shorter_width = abs(right_side - max_bin);
        elseif isempty(right_side) && ~isempty(left_side)
            shorter_width = abs(left_side - max_bin);
        elseif isempty(left_side) && isempty(right_side)
            break
        else
            shorter_width = min(abs(left_side - max_bin), abs(right_side - max_bin));
        end
        fwhm = shorter_width * 2;
        guess_std = fwhm / (2 * sqrt(2 * log(2)));
        
        if guess_std > max_gauss_width
            break
        end
        
        if max_height <= height_threshold
            break
        end

        ys = zeros(size(fit_data,1),1);
        ys = ys + guess_height * exp(-(fit_data(:,1)-guess_bin).^2 / (2*(guess_std.^2)));
        
        left = floor(guess_bin - guess_std);
        right = ceil(guess_bin + guess_std);
        ys(1:left+25-1) = 0;
        ys(right+25+1:end) = 0;
        ys(left+25:right+25) = ys(left+25:right+25) - (guess_height - min([ys(left+25) ys(right+25)]));
        
        if guess_bin < 50
            ys(101:125) = ys(1:25);
        else
            ys(26:50) = ys(126:150);
        end
        fit_data_(:,2) = fit_data_(:,2) - ys;
        fit_data_valid_ = fit_data_(26:125,:);

        % Stop searching for peaks
%         if max([fit_data(left+25,2) fit_data(right+25,2)]) <= std_threshold | guess_height < min_peak_height %Args.PeakThreshold * std(fit_data_valid_(:,2), 'omitnan')
% %             max([fit_data(left+25,2) fit_data(right+25,2)])
% %             std_threshold
%             break
%         end
        
        guess = [guess; [guess_bin, guess_std, guess_height]];
        
    end
    
    % Check overlaps of gaussians
    for i = 1:2
        if size(guess,1) > 1
            guess = sortrows(guess,1,'ascend');
            gauss_bounds = zeros(size(guess,1),2);
            for peak_no = 1:size(guess,1)
                gauss_bounds(peak_no,1) = (guess(peak_no,1) - (guess(peak_no,2) * gauss_overlap_threshold)); % Left bound
                gauss_bounds(peak_no,2) = (guess(peak_no,1) + (guess(peak_no,2) * gauss_overlap_threshold)); % Right bound
            end
            
            gauss_bounds = [gauss_bounds; gauss_bounds(1,:) + 100]; % Wrap around bin 100 to bin 1 for comparison
            guess_ = [guess; guess(1,:)];
            
            keep_peak_ = true(size(guess_,1), 1);
            for peak_no = 1:size(guess,1)
                if gauss_bounds(peak_no,2) > gauss_bounds(peak_no+1,1) % Check right bound of current peak with left bound of next peak
                    if guess_(peak_no,3) > guess_(peak_no+1,3) % Check height of current peak with height of next peak
                        keep_peak_(peak_no+1) = false; % Drop next peak
                    elseif guess_(peak_no,3) < guess_(peak_no+1,3)
                        keep_peak_(peak_no) = false; % Drop current peak
                    end
                end
            end
            
            keep_peak = keep_peak_(2:end-1);
            keep_peak = [keep_peak_(1) & keep_peak_(end); keep_peak];
            guess = guess(keep_peak, :);
        end
    end
    
    % Check gaussian width
%     if size(guess,1) > 1
%         keep_peak = true(size(guess,1), 1);
%         keep_peak = guess(:,2) < max_gauss_width;
%         guess = guess(keep_peak, :);
%     end
    
%     fit_data_ = fit_data;
%     fit_data_valid_ = fit_data_valid;
%     for i = 1:size(guess,1)
%         ys = zeros(size(fit_data,1),1);
%         ys = ys + guess(i,3) * exp(-(fit_data(:,1)-guess(i,1)).^2 / (2*(guess(i,2).^2)));
%         
%         left = floor(guess_bin - guess_std);
%         right = ceil(guess_bin + guess_std);
%         ys(1:left+25-1) = 0;
%         ys(right+25+1:end) = 0;
%         ys(left+25:right+25) = ys(left+25:right+25) - (guess_height - min([ys(left+25) ys(right+25)]));
%         
%         if guess_bin < 50
%             ys(101:125) = ys(1:25);
%         else
%             ys(26:50) = ys(126:150);
%         end
%         
%         fit_data_(:,2) = fit_data_(:,2) - ys;
%         fit_data_valid_ = fit_data_(26:125,:);
%     end
    
    binFiringRate = cellData.binFiringRate;
        
    for trial_no = 1:size(binFiringRate, 1)
        [trial_fits.("t" + trial_no), trial_threshold.("t" + trial_no)] = fit_trial(binFiringRate(trial_no,:), height_threshold);
    end
    
    data.trial_fits = trial_fits;
    data.trial_threshold = trial_threshold;
    guess_original = guess;
    
    if size(guess,1) > 4
        guess = [];
    end
    
    if vmr_threshold < 0.4
        guess = [];
    end
    
%     figure;
%     hold all
%     title('Firing Rate Map'); xlabel('Position Bins'); ylabel('Firing Rate (spike/sec)');
%     plot(basemapLrw, 'g'); legends{1} = sprintf('Raw');
%     
%     plot(basemapLsm, 'b'); legends{2} = sprintf('AdSmooth');
%     
%     plot(mapLsm_FFT_Low_Pass, 'r'); legends{3} = sprintf('AdSm-FFT-Low-Pass');
    
    if ~isempty(guess)
            
%         keep_guess = false(size(guess,1),1);
%         fns = fieldnames(trial_fits);
%         for i = 1:size(guess,1)
%             trial_peak_check = false(size(binFiringRate,1),1);
%             for j = 1:size(binFiringRate,1)
%                 trial_peak_ind = trial_fits.(fns{j})(trial_fits.(fns{j})(:,2) ~= 0,1);
%                 trial_peak_multi_check = true(size(trial_peak_ind,1),1);
%                 trial_peak_multi_check(trial_peak_ind < (guess(i,1) - guess(i,2))) = false;
%                 trial_peak_multi_check(trial_peak_ind > (guess(i,1) + guess(i,2))) = false;
%                 trial_peak_check(j) = any(trial_peak_multi_check,1);
%             end
%             if sum(trial_peak_check) > size(trial_peak_check,1)*Args.TrialPeakCheckFrac
%                 keep_guess(i) = true;
%             end
%         end
%         guess = guess(keep_guess,:);
        
%         figure;
%         hold all
%         title('Firing Rate Map'); xlabel('Position Bins'); ylabel('Firing Rate (spike/sec)');
%         plot(basemapLrw, 'g'); legends{1} = sprintf('Raw');
%         
%         plot(basemapLsm, 'b'); legends{2} = sprintf('AdSmooth');
%         
%         plot(mapLsm_FFT_Low_Pass, 'r'); legends{3} = sprintf('AdSm-FFT-Low-Pass');
%         
%         mu_lines = [guess(:,1) [1:size(guess,1)]'];
%         mu_lines_labels = cellstr("Field " + num2str(mu_lines(:,2)));
%         xline(mu_lines(:,1), '--', mu_lines_labels, 'Color', '#D95319');
%         ax = gca;
%         for peak_no = 1:size(guess,1)
%             left = guess(peak_no,1) - guess(peak_no,2);
%             right = guess(peak_no,1) + guess(peak_no,2);
%             x = [left right right left];
%             y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];
%             patch(x, y, 'k', 'FaceAlpha', '0.1', 'LineStyle', 'none');
%         end
        
        gaus = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;
        
        y_separated = zeros(size(guess,1),100);
        k = zeros(size(guess,1),100);
        
        for i = 1:size(guess,1)
            x = linspace(-24,125,150);
            mu = guess(i,1);
            sig = guess(i,2);
            amp = guess(i,3);
            vo = 0;
            y = gaus(x,mu,sig,amp,vo);
            if mu < 50
                y(101:125) = y(1:25);
            else
                y(26:50) = y(126:150);
            end
            % Plot gaussian
%             plot(x(26:125), y(26:125), 'k-', 'LineWidth',3);
%             hold on

            % Calculate log-likelihood
            y_separated(i,:) = y(26:125) ./ sum(y(26:125));
            k(i,:) = y(26:125) > 0.01;

        end
%             axis([1 100 -inf inf]);
%             yline((baseline_mean_height + peak_height_from_baseline/2),'r');
%             % yline(Args.PeakThreshold * std(fit_data_valid_(:,2), 'omitnan'),'b');
%             yline(std_threshold,'b');
        
        x = linspace(-24,125,150);        
        y = 0;
        for i = 1:size(guess,1)
            y = gaus(x,guess(i,1),guess(i,2),guess(i,3),vo) + y;
        end
        
%         legends{4} = ''; legends{5} = ''; legends{6} = ''; legends{7} = ''; legends{8} = ''; legends{9} = '';
%         plot(x(26:125), y(26:125), 'k-', 'LineWidth',3); legends{10} = sprintf('GMM');
%         legend(legends, 'Location', 'northeastoutside');
        
%         figure;
%         imagesc(binFiringRate); title('Firing Rates'); xlabel('Position Bins'); ylabel('Trials');
%         
%             mu_lines = [guess(:,1) [1:size(guess,1)]'];
%             mu_lines_labels = cellstr("Field " + num2str(mu_lines(:,2)));
%             xline(mu_lines(:,1), '--', mu_lines_labels, 'Color', '#D95319');
%             ax = gca;
%             for peak_no = 1:size(guess,1)
%                 left = guess(peak_no,1) - guess(peak_no,2);
%                 right = guess(peak_no,1) + guess(peak_no,2);
%                 x_ = [left right right left];
%                 y_ = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];
%                 patch(x_, y_, 'k', 'FaceAlpha', '0.1', 'LineStyle', 'none');
%             end
        
        if ~isempty(guess)
            % Calculate accuracy metrics
            rmse = sqrt(mean((basemapLrw - y(26:125)).^2,'omitnan'));
            nrmse_mean = rmse / mean(basemapLrw, 'omitnan');
            nrmse_stdev = rmse / std(basemapLrw, 'omitnan');
            auc = mean(abs(y ./ mapLsm_FFT_Low_Pass_padded), 'omitnan');
            rss = sum((mapLsm_FFT_Low_Pass - y(26:125)).^2, 'omitnan');
            r2 = 1 - (rss/sum((mapLsm_FFT_Low_Pass - mean(mapLsm_FFT_Low_Pass,'omitnan')).^2, 'omitnan'));
            aic_old = -2*r2*100 + 2*size(guess,1); %100*log(rss/100) + 2*size(guess,1);
            bic_old = -2*r2*100 + size(guess,1)*log(100); %100*log(rss/100) + size(guess,1)*log(100);
            log_likelihood = sum(k .* log(y_separated),'all','omitnan');
            aic = -2*log_likelihood + 2*size(guess,1);
            bic = -2*log_likelihood + size(guess,1)*log(100);
            
            guess(:,4) = (guess(:,3) - mean(fit_data_(26:125, 2),'omitnan')) / std(fit_data_(26:125, 2),'omitnan'); % Height z-score
            
            data.basemapLrw = basemapLrw;
            data.basemapLsm = basemapLsm;
            data.mapLsm_FFT_Low_Pass = mapLsm_FFT_Low_Pass;
            data.GMM = guess;
            data.GMM_original = guess_original;
            data.baseline = fit_data_;
            data.rmse = rmse;
            data.nrmse_mean = nrmse_mean;
            data.nrmse_stdev = nrmse_stdev;
            data.auc = auc;
            data.log_likelihood = log_likelihood;
            data.aic = aic_old;%aic;
            data.bic = bic_old;%bic;
        end
    end
       
    if isempty(guess)
        data.basemapLrw = basemapLrw;
        data.basemapLsm = basemapLsm;
        data.mapLsm_FFT_Low_Pass = mapLsm_FFT_Low_Pass;
        data.GMM = guess;
        data.GMM_original = guess_original;
        data.baseline = fit_data_;
        data.rmse = NaN;
        data.nrmse_mean = NaN;
        data.nrmse_stdev = NaN;
        data.auc = NaN;
        data.log_likelihood = NaN;
        data.aic = NaN;
        data.bic = NaN;
    end
    
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

function [trial_fits, trial_threshold] = fit_trial(ratemap, TrialThreshold)
    
%     for trial_no = 1:size(binFiringRate, 1)
%     ratemap = binFiringRate(trial_no,:);
    mapLsm_FFT_Low_Pass_padded = [ratemap(76:100) ratemap ratemap(1:25)];
    
    fit_data = [-24:1:125; mapLsm_FFT_Low_Pass_padded]';
    fit_data_valid = fit_data(26:125,:);
    
    fit_data_ = fit_data;
    fit_data_valid_ = fit_data_valid;
    %min_peak_height = max(fit_data_valid_(:,2), [], 'omitnan') * 0.5;
  
    % Check peak heights compared to baseline
    %std_threshold = TrialThreshold * std(fit_data_valid_(:,2), 'omitnan');
    
    fit_data_ = fit_data;
    fit_data_valid_ = fit_data_valid;
    %fit_data_valid_(fit_data_valid_(:,2) < std_threshold,2) = 0;
    fit_data_valid_(fit_data_valid_(:,2) < TrialThreshold) = 0;
    %disp(std_threshold);
        
    trial_fits = fit_data_valid_;
    %trial_std_threshold = std_threshold;
    trial_threshold = TrialThreshold;
%     trial_fits.("t" + trial_no) = fit_data_valid_;
%     trial_std_threshold.("t" + trial_no) = std_threshold;
%     end