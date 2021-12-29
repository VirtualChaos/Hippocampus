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
				'BinSize',100, 'ThresVel',1);
            
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
    peakrate_full = nanmax(basemapLsm);
    peakrate_set = peakrate_full;
    
    Fs = 100;            % Sampling frequency
    T = 1/Fs;             % Sampling period
    L = 100;             % Length of signal
    s = (0:L-1)*T;        % Space vector
    
    Y = fft(basemapLsm);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
%     figure;
%     plot(f,P1);
%     title('Single-Sided Amplitude Spectrum of X(t)');
%     xlabel('f (Hz)');
%     ylabel('|P1(f)|');
    
    % figure;
    filt_freq_arr = 0:2:22;
    % Low-Pass Filter
    for filt_freq_idx = 1:12
        Y_filt = Y;
        Y_filt(find(f > filt_freq_arr(filt_freq_idx),1):end-find(f > filt_freq_arr(filt_freq_idx),1)) = 0; % Low-Pass filter at 150Hz
        mapLsm_FFT_Low_Pass = abs(ifft(Y_filt));
        
        % subplot(4,3,filt_freq_idx);
        % plot(basemapLrw, 'g');legends{1} = sprintf('Raw');
        % hold on
        % plot(basemapLsm, 'b');legends{2} = sprintf('AdSmooth');
        % hold on
        % plot(mapLsm_FFT_Low_Pass, 'r');legends{3} = sprintf('AdSm-FFT-Low-Pass');
        % %legend(legends);
        % title("filt_freq =" + filt_freq_arr(filt_freq_idx) + "Hz")
    end
    
    filt_freq = f(find((cumsum(P1) ./ sum(P1) > 0.75), 1, 'first')); % 75% threshold
    Y_filt = Y;
    Y_filt(find(f > filt_freq,1):end-find(f > filt_freq,1)) = 0; % Low-Pass filter at 150Hz
    mapLsm_FFT_Low_Pass = abs(ifft(Y_filt));
    
%     figure;
%     plot(basemapLrw, 'g');legends{1} = sprintf('Raw');
%     hold on
%     plot(basemapLsm, 'b');legends{2} = sprintf('AdSmooth');
%     hold on
%     plot(mapLsm_FFT_Low_Pass, 'r');legends{3} = sprintf('AdSm-FFT-Low-Pass');
%     %legend(legends);
%     title("filt_freq =" + filt_freq + "Hz")
      
    % GMM Iterative Method
    
    % Interpolate
    mapLsm_FFT_Low_Pass_padded = [mapLsm_FFT_Low_Pass(75:100) mapLsm_FFT_Low_Pass mapLsm_FFT_Low_Pass(1:26)];
    % mapLsm_FFT_Low_Pass_padded = [basemapLsm(75:100) basemapLsm basemapLsm(1:26)];
    % vq = interp1(1:152,mapLsm_FFT_Low_Pass_padded,0.1:0.1:152,'linear');
    % vq = vq(11:1510);
    %
    % fit_data = [-24.9:0.1:125; vq]';
    % fit_data_valid = fit_data(251:1250,:);
    
    fit_data = [-25:1:126; mapLsm_FFT_Low_Pass_padded]';
    fit_data_valid = fit_data(27:126,:);
    
    guess = []; % Mean (Centre Bin), Std (Half-Width), Height (Rate)
    fit_data_ = fit_data;
    fit_data_valid_ = fit_data_valid;
    peak_threshold = 3; % Multiples of stdev of height
    min_peak_height = max(fit_data_valid_(:,2), [], 'omitnan') * 0.5;
    min_gauss_std = 0;
    max_gauss_std = 12.5;
    gauss_overlap_threshold = 0.5;
    
    % . Find peak: Loop through, finding a candidate peak, and fitting with a guess gaussian
    %   Stopping procedures: limit on # of peaks, or relative or absolute height thresholds
    while size(guess,1) < 3
        % Find max peak
        [max_height, max_bin_idx] = max(fit_data_valid_(:,2), [], 'omitnan');
        max_bin = fit_data_valid_(max_bin_idx);
        % Stop searching for peaks once height drops below height threshold
        if max_height <= peak_threshold * std(fit_data_valid_(:,2), 'omitnan')
            break
        end
        
        % Set the guess parameters for gaussian fitting
        guess_bin = max_bin; % Mean
        guess_height = max_height; % Height
        
        % Halt fitting process if peak drops below minimum height
        if guess_height < min_peak_height
            break
        end
        
        % Find half height index on each side of the center bin (gaussian mean)
        half_height = 0.5 * max_height;
        max_bin_idx_ = find(fit_data_(:,1) == max_bin);
        left_idx = find(fit_data_(1:max_bin_idx_,2) <= half_height, 1, 'last');
        left_side = fit_data_(left_idx,1);
        right_idx = find(fit_data_(max_bin_idx_:end,2) <= half_height, 1, 'first') + max_bin_idx_ - 1;
        right_side = fit_data_(right_idx,1);
        
        % Use shorter side
        if isempty(left_side)
            shorter_width = abs(right_side - max_bin);
        elseif isempty(right_side)
            shorter_width = abs(left_side - max_bin);
        elseif isempty(left_side) && isempty(right_side)
            shorter_width = 10000;
        else
            shorter_width = min(abs(left_side - max_bin), abs(right_side - max_bin));
        end
        fwhm = shorter_width * 2;
        guess_std = fwhm / (2 * sqrt(2 * log(2)));
        
        % Check that gaussian std are within limits
%         if guess_std < min_gauss_std
%             guess_std = min_gauss_std;
%         end
%         if guess_std > max_gauss_std
%             guess_std = max_gauss_std;
%         end
        
        guess = [guess; [guess_bin, guess_std, guess_height]];
        
        ys = zeros(size(fit_data,1),1);
        ys = ys + guess_height * exp(-(fit_data(:,1)-guess_bin).^2 / (2*(guess_std.^2)));
        fit_data_(:,2) = fit_data_(:,2) - ys;
        %fit_data_valid_ = fit_data_(251:1250,:);
        fit_data_valid_ = fit_data_(27:126,:);
    end
    
    % Check overlaps of gaussians
    if size(guess,1) > 1
        %     guess = sortrows(guess);
        %     gauss_bounds = zeros(size(guess,1),2);
        %     for peak_no = 1:size(guess,1)
        %         gauss_bounds(peak_no,1) = (guess(peak_no,1) - (guess(peak_no,2) * gauss_overlap_threshold)); % Left bound
        %         gauss_bounds(peak_no,2) = (guess(peak_no,1) + (guess(peak_no,2) * gauss_overlap_threshold)); % Right bound
        %     end
        %
        %     gauss_bounds = [gauss_bounds; gauss_bounds(1,:) + 100]; % Wrap around bin 100 to bin 1 for comparison
        %     guess_ = [guess; guess(1,:)];
        %
        %     keep_peak_ = true(size(guess_,1), 1);
        %     for peak_no = 1:size(guess,1)
        %         if gauss_bounds(peak_no,2) > gauss_bounds(peak_no+1,1) % Check right bound of current peak with left bound of next peak
        %             if guess_(peak_no,3) > guess_(peak_no+1,3) % Check height of current peak with height of next peak
        %                 keep_peak_(i+1) = false; % Drop next peak
        %             elseif guess_(peak_no,3) < guess_(peak_no+1,3)
        %                 keep_peak_(i) = false; % Drop current peak
        %             end
        %         end
        %     end
        %
        %     keep_peak = keep_peak_(2:end-1);
        %     keep_peak = [keep_peak_(1) & keep_peak_(end); keep_peak];
        %     guess = guess(keep_peak, :);
    end % Commented out
    
    fit_data_ = fit_data;
    fit_data_valid_ = fit_data_valid;
    for i = 1:size(guess,1)
        ys = zeros(size(fit_data,1),1);
        ys = ys + guess(i,3) * exp(-(fit_data(:,1)-guess(i,1)).^2 / (2*(guess(i,2).^2)));
        fit_data_(:,2) = fit_data_(:,2) - ys;
        %fit_data_valid_ = fit_data_(251:1250,:);
        fit_data_valid_ = fit_data_(27:126,:);
    end
    
    % Check peak heights compared to baseline
    baseline_mean_height = mean(fit_data_valid_(:,2), 'omitnan');
    peak_height_from_baseline = max(fit_data_valid(:,2), [], 'omitnan') - baseline_mean_height;
    guess = guess(guess(:,3) > (baseline_mean_height + peak_height_from_baseline/2),:);
    
%     figure;
%     plot(fit_data_valid(:,1), fit_data_valid(:,2), 'b');
%     hold on
%     % plot(fit_data_(251:1250, 1), fit_data_(251:1250, 2));
%     plot(fit_data_(27:126, 1), fit_data_(27:126, 2),'r');
    
    figure;
    hold all
    plot(basemapLrw, 'g'); legends{1} = sprintf('Raw');
    
    plot(basemapLsm, 'b'); legends{2} = sprintf('AdSmooth');
    
    plot(mapLsm_FFT_Low_Pass, 'r'); legends{3} = sprintf('AdSm-FFT-Low-Pass');
    legend(legends, 'Location', 'northeastoutside');
    
    mu_lines = [guess(:,1) [1:size(guess,1)]'];
    mu_lines_labels = cellstr("Field " + num2str(mu_lines(:,2)));
    xline(mu_lines(:,1), '--', mu_lines_labels, 'Color', '#D95319');
    ax = gca;
    for peak_no = 1:size(guess,1)
        left = guess(peak_no,1) - guess(peak_no,2);
        right = guess(peak_no,1) + guess(peak_no,2);
        x = [left right right left];
        y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];
        patch(x, y, 'k', 'FaceAlpha', '0.1', 'LineStyle', 'none');
    end
    
    gaus = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;
    
    for i = 1:size(guess,1)
        x = linspace(-25,126,152);
        mu = guess(i,1);
        sig = guess(i,2);
        amp = guess(i,3);
        vo = 0;
        y = gaus(x,mu,sig,amp,vo);
        if mu < 50
            y(101:126) = y(1:26);
        else
            y(27:52) = y(127:152);
        end
        % Plot gaussian
%         hold all
%         plot(x(27:126), y(27:126), 'k-', 'LineWidth',3);
    end
%     axis([1 100 -inf inf]);
%     yline((baseline_mean_height + peak_height_from_baseline/2));
%     hold off
    
    x = linspace(-25,126,152);
    if size(guess,1) == 1
        y = gaus(x,guess(1,1),guess(1,2),guess(1,3),vo);
    elseif size(guess,1) == 2
        y = gaus(x,guess(1,1),guess(1,2),guess(1,3),vo) + gaus(x,guess(2,1),guess(2,2),guess(2,3),vo);
    else
        y = gaus(x,guess(1,1),guess(1,2),guess(1,3),vo) + gaus(x,guess(2,1),guess(2,2),guess(2,3),vo) + gaus(x,guess(3,1),guess(3,2),guess(3,3),vo);
    end
    
    plot(x(27:126), y(27:126), 'k-', 'LineWidth',3);
    
    % Calculate RMSE
    rmse = sqrt(mean((mapLsm_FFT_Low_Pass_padded - y).^2,'omitnan'));
    nrmse = rmse / mean(mapLsm_FFT_Low_Pass_padded, 'omitnan');
    auc = mean(abs(y ./ mapLsm_FFT_Low_Pass_padded), 'omitnan');
    
    data.basemapLrw = basemapLrw;
    data.basemapLsm = basemapLsm;
    data.mapLsm_FFT_Low_Pass = mapLsm_FFT_Low_Pass;
    data.GMM = guess;
    data.rmse = rmse;
    data.nrmse = nrmse;
    data.auc = auc;
    
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
