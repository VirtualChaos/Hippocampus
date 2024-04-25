function [obj, varargout] = plot(obj,varargin)
%@dirfiles/plot Plot function for dirfiles object.
%   OBJ = plot(OBJ) creates plots for the neural response

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		  'ReturnVars',{''}, 'ArgsOnly',0, 'Cmds','', ...
          'FiringRateMapRaw',0, 'FiringRateMapAdSm',0, 'BinFiringRate',0, ...
          'TrialBinFiringRate', 0, 'Baseline',0, 'GMM',0, 'TrialFits',0, 'AlphaAdSm',0, ...
          'FluorescenceTrace',0, 'TrialFluorescenceTrace',0, 'PlaceFieldPositionEntropy',0);
Args.flags = {'LabelsOff','ArgsOnly','FiringRateMapRaw','FiringRateMapAdSm','BinFiringRate', ...
              'TrialBinFiringRate','Baseline','GMM', 'TrialFits', 'AlphaAdSm', ...
              'FluorescenceTrace','TrialFluorescenceTrace','PlaceFieldPositionEntropy'};
[Args,varargin2] = getOptArgs(varargin,Args);

% if user select 'ArgsOnly', return only Args structure for an empty object
if Args.ArgsOnly
    Args = rmfield (Args, 'ArgsOnly');
    varargout{1} = {'Args',Args};
    return;
end

if(~isempty(Args.NumericArguments))
	n = Args.NumericArguments{1};
	if(Args.FiringRateMapRaw)
		% Raw Firing Rate Map
        cell_no = n;
        fns = fieldnames(obj.data.placefieldData);

        subplot(3,1,1);
        imagesc(obj.data.cellData.(fns{cell_no}).maps_raw);
        axis off
        colorbar;
        subplot(3,1,2);
        imagesc(obj.data.cellData.(fns{cell_no}).maps_raw1);
        axis off
        colorbar;
        subplot(3,1,3);
        imagesc(obj.data.cellData.(fns{cell_no}).maps_raw2);
        axis off
        colorbar;
    elseif(Args.FiringRateMapAdSm)
        % Adaptive Smoothed Firing Rate Map
        cell_no = n;
        fns = fieldnames(obj.data.placefieldData);

        subplot(3,1,1);
        imagesc(obj.data.cellData.(fns{cell_no}).maps_adsm);
        axis off
        colorbar;
        subplot(3,1,2);
        imagesc(obj.data.cellData.(fns{cell_no}).maps_adsm1);
        axis off
        colorbar;
        subplot(3,1,3);
        imagesc(obj.data.cellData.(fns{cell_no}).maps_adsm2);
        axis off
        colorbar;
    elseif(Args.BinFiringRate)
        
        subplot(3,1,1);
        % Adaptive Smoothed Firing Rate Map
        cell_no = n;
        cellFiringRate = squeeze(obj.data.binFiringRate(cell_no,:,:));
        plot(mean(cellFiringRate),'k','LineWidth',2);
        hold on
        plot(cellFiringRate','Color',[0.5 0.5 0.5 0.2])
        ylabel('Firing rate','FontSize',12);
%         ylim([0 2.5]);
        
%         text(50,2.15,num2str(mean(cellFiringRate,'all')))
%         text(50,2,num2str(max(mean(cellFiringRate))))
%         text(50,1.85,num2str(mean(obj.data.spiketrain(cell_no,:))))
        
        subplot(3,1,2);
        imagesc(cellFiringRate)
        xlabel('Position Bin','FontSize',12); ylabel('Trial','FontSize',12);
        
        subplot(3,1,3);
        plot(mean(cellFiringRate),'k','LineWidth',2);
        ylabel('Firing rate','FontSize',12);

%         subplot(3,1,2);
%         % Adaptive Smoothed Firing Rate Map
%         cell_no = n;
%         fns = fieldnames(obj.data.placefieldData);
%         
%         binFiringRate = obj.data.cellData.(fns{cell_no}).binFiringRate;
%         binFiringRate(binFiringRate < obj.data.placefieldData.(fns{cell_no}).mean_threshold + obj.data.placefieldData.(fns{cell_no}).std_threshold) = 0;
%         imagesc(binFiringRate); title('Trial Firing Rates'); xlabel('Position Bins'); ylabel('Trials');
%         
%         %imagesc(obj.data.cellData.(fns{cell_no}).binFiringRate); title('Trial Firing Rates'); xlabel('Position Bins'); ylabel('Trials');
%         hold on
%         %colorbar;
%         hold on
        
%         fns = fieldnames(obj.data.placefieldData);
%         
%         if ~isempty(obj.data.placefieldData.(fns{cell_no}).GMM) & obj.data.isplacecell(cell_no,4)
%             mu_lines = [obj.data.placefieldData.(fns{cell_no}).GMM(:,1) [1:size(obj.data.placefieldData.(fns{cell_no}).GMM,1)]'];
%             mu_lines_labels = cellstr("Field " + num2str(mu_lines(:,2)));
%             xline(mu_lines(:,1), '--', mu_lines_labels, 'Color', '#D95319'); hold on
%             ax = gca;
%             for peak_no = 1:size(obj.data.placefieldData.(fns{cell_no}).GMM,1)
%                 left = obj.data.placefieldData.(fns{cell_no}).GMM(peak_no,1) - obj.data.placefieldData.(fns{cell_no}).GMM(peak_no,2);
%                 right = obj.data.placefieldData.(fns{cell_no}).GMM(peak_no,1) + obj.data.placefieldData.(fns{cell_no}).GMM(peak_no,2);
%                 x = [left right right left];
%                 y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];
%                 patch(x, y, 'w', 'FaceAlpha', '0.6', 'LineStyle', 'none'); hold on
%             end
%         end
%         hold off
%         
%         % Gaussian Mixed Model
%         cell_no = n;
%         fns = fieldnames(obj.data.placefieldData);
%         
%         subplot(2,1,2);
%         title('Firing Rate Map'); xlabel('Position (Bin)'); ylabel('Firing Rate');
%         plot(obj.data.placefieldData.(fns{cell_no}).basemapLrw, 'g'); legends{1} = sprintf('Raw'); hold on
%         
%         plot(obj.data.placefieldData.(fns{cell_no}).basemapLsm, 'b'); legends{2} = sprintf('AdSmooth'); hold on
%         
%         %plot(imgaussfilt(obj.data.placefieldData.(fns{cell_no}).basemapLsm, 3), 'c'); legends{2} = sprintf('AdSmooth'); hold on
%                 
%         plot(obj.data.placefieldData.(fns{cell_no}).mapLsm_FFT_Low_Pass, 'r'); legends{3} = sprintf('AdSm-FFT-Low-Pass'); hold on
%         
%         axis([1 100 -inf inf]);
%         yline(obj.data.placefieldData.(fns{cell_no}).mean_threshold,'b-'); hold on
%         yline(obj.data.placefieldData.(fns{cell_no}).mean_threshold + obj.data.placefieldData.(fns{cell_no}).std_threshold,'b--'); hold on
%         yline(obj.data.placefieldData.(fns{cell_no}).mean_threshold - obj.data.placefieldData.(fns{cell_no}).std_threshold,'b--'); hold on
% %         yline((obj.data.placefieldData.(fns{cell_no}).baseline_mean_height + obj.data.placefieldData.(fns{cell_no}).peak_height_from_baseline/2),'r'); hold on
% %         yline(obj.data.placefieldData.(fns{cell_no}).std_threshold,'b'); hold on
% %         yline(obj.data.placefieldData.(fns{cell_no}).std_threshold / obj.data.placefieldData.(fns{cell_no}).Args.PeakThreshold,'b--'); hold on
%                 
%         text(3, 1.3*max(obj.data.placefieldData.(fns{cell_no}).basemapLrw), {"Mean: " + obj.data.placefieldData.(fns{cell_no}).mean_threshold, "StDev: " + obj.data.placefieldData.(fns{cell_no}).std_threshold, ...
%             "VMR: " + obj.data.placefieldData.(fns{cell_no}).vmr_threshold},'FontSize',12); %"CV: " + obj.data.placefieldData.(fns{cell_no}).cv_threshold,
%         ylim([0 1.7*max(obj.data.placefieldData.(fns{cell_no}).basemapLrw)])
%         
%         if ~isempty(obj.data.placefieldData.(fns{cell_no}).GMM) & obj.data.isplacecell(cell_no,4)
%             mu_lines = [obj.data.placefieldData.(fns{cell_no}).GMM(:,1) [1:size(obj.data.placefieldData.(fns{cell_no}).GMM,1)]'];
%             mu_lines_labels = cellstr("Field " + num2str(mu_lines(:,2)));
%             xline(mu_lines(:,1), '--', mu_lines_labels, 'Color', '#D95319'); hold on
%             ax = gca;
%             for peak_no = 1:size(obj.data.placefieldData.(fns{cell_no}).GMM,1)
%                 left = obj.data.placefieldData.(fns{cell_no}).GMM(peak_no,1) - obj.data.placefieldData.(fns{cell_no}).GMM(peak_no,2);
%                 right = obj.data.placefieldData.(fns{cell_no}).GMM(peak_no,1) + obj.data.placefieldData.(fns{cell_no}).GMM(peak_no,2);
%                 x = [left right right left];
%                 y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];
%                 patch(x, y, 'k', 'FaceAlpha', '0.1', 'LineStyle', 'none'); hold on
%             end
%             
%             gaus = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;
%             
%             for i = 1:size(obj.data.placefieldData.(fns{cell_no}).GMM,1)
%                 x = linspace(-24,125,150);
%                 mu = obj.data.placefieldData.(fns{cell_no}).GMM(i,1);
%                 sig = obj.data.placefieldData.(fns{cell_no}).GMM(i,2);
%                 amp = obj.data.placefieldData.(fns{cell_no}).GMM(i,3);
%                 vo = 0;
%                 y = gaus(x,mu,sig,amp,vo);
%                 if mu < 50
%                     y(101:125) = y(1:25);
%                 else
%                     y(26:50) = y(126:150);
%                 end
%                 % Plot gaussian
%                 plot(x(26:125), y(26:125), 'k-', 'LineWidth',3);
%             end
%             
%             x = linspace(-24,125,150);         
%             y = 0;
%             for i = 1:size(obj.data.placefieldData.(fns{cell_no}).GMM,1)
%                 y = gaus(x,obj.data.placefieldData.(fns{cell_no}).GMM(i,1),obj.data.placefieldData.(fns{cell_no}).GMM(i,2),obj.data.placefieldData.(fns{cell_no}).GMM(i,3),vo) + y;
%             end
%             
%             legends{4} = ''; legends{5} = ''; legends{6} = ''; legends{7} = ''; legends{8} = ''; legends{9} = '';
%             %plot(x(26:125), y(26:125), 'k-', 'LineWidth',3); legends{10} = sprintf('GMM'); hold on
%             %legend(legends, 'Location', 'northeastoutside');
%         end
%         hold off
        
    elseif(Args.TrialBinFiringRate)
        delete(findall(gcf,'type','annotation'));
        cell_no = floor((n-1)/obj.data.nTrials) + 1;
        trial = n - ((cell_no-1)*obj.data.nTrials);
        fns = fieldnames(obj.data.placefieldData);
        % fprintf("Cell: %d; Trial: %d\n",cell_no,trial);
        
        plot(obj.data.cellData.(fns{cell_no}).binFiringRate(trial,:), 'b'); title('Trial Firing Rates'); xlabel('Position Bins'); ylabel('Firing Rate (spikes/sec)');
        annotation('textbox', [0.005 0.85 0.095 0.07], 'String', sprintf("Cell: %d\nTrial: %d\n",cell_no,trial));
        
    elseif(Args.Baseline)
        % Baseline Firing
        cell_no = n;
        fns = fieldnames(obj.data.placefieldData);
        
        plot(1:100, obj.data.placefieldData.(fns{cell_no}).mapLsm_FFT_Low_Pass, 'b'); legends{1} = sprintf('Firing Rate Map');
        hold on
        plot(obj.data.placefieldData.(fns{cell_no}).baseline(26:125, 1), obj.data.placefieldData.(fns{cell_no}).baseline(26:125, 2),'r'); legends{2} = sprintf('Estimated Baseline');
        legend(legends);
        hold off
        
    elseif(Args.GMM)
        
        % Gaussian Mixed Model
        cell_no = n;
        fns = fieldnames(obj.data.placefieldData);
        
        %subplot(2,1,1);
        title('Firing Rate Map'); xlabel('Position Bins'); ylabel('Firing Rate (spike/sec)');
        plot(obj.data.placefieldData.(fns{cell_no}).basemapLrw, 'g'); legends{1} = sprintf('Raw'); hold on
        
        plot(obj.data.placefieldData.(fns{cell_no}).basemapLsm, 'b'); legends{2} = sprintf('AdSmooth'); hold on
        
        plot(obj.data.placefieldData.(fns{cell_no}).mapLsm_FFT_Low_Pass, 'r'); legends{3} = sprintf('AdSm-FFT-Low-Pass'); hold on
        
        if ~isempty(obj.data.placefieldData.(fns{cell_no}).GMM) & obj.data.isplacecell(cell_no,5)
            mu_lines = [obj.data.placefieldData.(fns{cell_no}).GMM(:,1) [1:size(obj.data.placefieldData.(fns{cell_no}).GMM,1)]'];
            mu_lines_labels = cellstr("Field " + num2str(mu_lines(:,2)));
            xline(mu_lines(:,1), '--', mu_lines_labels, 'Color', '#D95319'); hold on
            ax = gca;
            for peak_no = 1:size(obj.data.placefieldData.(fns{cell_no}).GMM,1)
                left = obj.data.placefieldData.(fns{cell_no}).GMM(peak_no,1) - obj.data.placefieldData.(fns{cell_no}).GMM(peak_no,2);
                right = obj.data.placefieldData.(fns{cell_no}).GMM(peak_no,1) + obj.data.placefieldData.(fns{cell_no}).GMM(peak_no,2);
                x = [left right right left];
                y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];
                patch(x, y, 'k', 'FaceAlpha', '0.1', 'LineStyle', 'none'); hold on
            end
            
            gaus = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;
            
            for i = 1:size(obj.data.placefieldData.(fns{cell_no}).GMM,1)
                x = linspace(-25,126,152);
                mu = obj.data.placefieldData.(fns{cell_no}).GMM(i,1);
                sig = obj.data.placefieldData.(fns{cell_no}).GMM(i,2);
                amp = obj.data.placefieldData.(fns{cell_no}).GMM(i,3);
                vo = 0;
                y = gaus(x,mu,sig,amp,vo);
                if mu < 50
                    y(101:125) = y(1:25);
                else
                    y(26:50) = y(126:150);
                end
                % Plot gaussian
                % plot(x(26:125), y(26:125), 'k-', 'LineWidth',3);
            end
            axis([1 100 -inf inf]);
            yline((obj.data.placefieldData.(fns{cell_no}).baseline_mean_height + obj.data.placefieldData.(fns{cell_no}).peak_height_from_baseline/2),'r'); hold on
            yline(obj.data.placefieldData.(fns{cell_no}).std_threshold,'b'); hold on
            
            x = linspace(-25,126,152);         
            y = 0;
            for i = 1:size(obj.data.placefieldData.(fns{cell_no}).GMM,1)
                y = gaus(x,obj.data.placefieldData.(fns{cell_no}).GMM(i,1),obj.data.placefieldData.(fns{cell_no}).GMM(i,2),obj.data.placefieldData.(fns{cell_no}).GMM(i,3),vo) + y;
            end
            
            legends{4} = ''; legends{5} = ''; legends{6} = ''; legends{7} = ''; legends{8} = ''; legends{9} = '';
            plot(x(26:125), y(26:125), 'k-', 'LineWidth',3); legends{10} = sprintf('GMM'); hold on
            legend(legends, 'Location', 'northeastoutside');
        end
        hold off
        
%         subplot(2,1,2);
%         title('Firing Rate Map'); xlabel('Position Bins'); ylabel('Firing Rate (spike/sec)');
%         plot(obj.data.placefieldData.(fns{cell_no}).basemapLrw, 'g'); legends{1} = sprintf('Raw'); hold on
%         
%         plot(obj.data.placefieldData.(fns{cell_no}).basemapLsm, 'b'); legends{2} = sprintf('AdSmooth'); hold on
%         
%         plot(obj.data.placefieldData.(fns{cell_no}).mapLsm_FFT_Low_Pass, 'r'); legends{3} = sprintf('AdSm-FFT-Low-Pass'); hold on
%         
%         if ~isempty(obj.data.placefieldData.(fns{cell_no}).GMM_original)
%             mu_lines = [obj.data.placefieldData.(fns{cell_no}).GMM_original(:,1) [1:size(obj.data.placefieldData.(fns{cell_no}).GMM_original,1)]'];
%             mu_lines_labels = cellstr("Field " + num2str(mu_lines(:,2)));
%             xline(mu_lines(:,1), '--', mu_lines_labels, 'Color', '#D95319'); hold on
%             ax = gca;
%             for peak_no = 1:size(obj.data.placefieldData.(fns{cell_no}).GMM_original,1)
%                 left = obj.data.placefieldData.(fns{cell_no}).GMM_original(peak_no,1) - obj.data.placefieldData.(fns{cell_no}).GMM_original(peak_no,2) *1.5;
%                 right = obj.data.placefieldData.(fns{cell_no}).GMM_original(peak_no,1) + obj.data.placefieldData.(fns{cell_no}).GMM_original(peak_no,2) *1.5;
%                 x = [left right right left];
%                 y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];
%                 patch(x, y, 'k', 'FaceAlpha', '0.1', 'LineStyle', 'none'); hold on
%             end
%             
%             gaus = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;
%             
%             for i = 1:size(obj.data.placefieldData.(fns{cell_no}).GMM_original,1)
%                 x = linspace(-25,126,152);
%                 mu = obj.data.placefieldData.(fns{cell_no}).GMM_original(i,1);
%                 sig = obj.data.placefieldData.(fns{cell_no}).GMM_original(i,2);
%                 amp = obj.data.placefieldData.(fns{cell_no}).GMM_original(i,3);
%                 vo = 0;
%                 y = gaus(x,mu,sig,amp,vo);
%                 if mu < 50
%                   y(101:125) = y(1:25);
%                 else
%                   y(26:50) = y(126:150);
%                 end
%                 % Plot gaussian
%                 % plot(x(26:125), y(26:125), 'k-', 'LineWidth',3);
%             end
%             axis([1 100 -inf inf]);
%             yline((obj.data.placefieldData.(fns{cell_no}).baseline_mean_height + obj.data.placefieldData.(fns{cell_no}).peak_height_from_baseline/2),'r'); hold on
%             yline(obj.data.placefieldData.(fns{cell_no}).std_threshold,'b'); hold on
%             
%             x = linspace(-25,126,152);
%            
%             y = 0;
%             for i = 1:size(obj.data.placefieldData.(fns{cell_no}).GMM_original,1)
%                 y = gaus(x,obj.data.placefieldData.(fns{cell_no}).GMM_original(i,1),obj.data.placefieldData.(fns{cell_no}).GMM_original(i,2),obj.data.placefieldData.(fns{cell_no}).GMM_original(i,3),vo) + y;
%             end
%             
%             legends{4} = ''; legends{5} = ''; legends{6} = ''; legends{7} = ''; legends{8} = ''; legends{9} = '';
%             plot(x(26:125), y(26:125), 'k-', 'LineWidth',3); legends{10} = sprintf('GMM'); hold on
%             legend(legends, 'Location', 'northeastoutside');
%         end
%         hold off
        
    elseif(Args.TrialFits)
        delete(findall(gcf,'type','annotation'));
        cell_no = floor((n-1)/obj.data.nTrials) + 1;
        trial = n - ((cell_no-1)*obj.data.nTrials);
        fns = fieldnames(obj.data.placefieldData);
        % fprintf("Cell: %d; Trial: %d\n",cell_no,trial);
        
        plot(obj.data.cellData.(fns{cell_no}).binFiringRate(trial,:), 'b'); title('Trial Firing Rates'); xlabel('Position Bins'); ylabel('Firing Rate (spikes/sec)');
        yline(obj.data.placefieldData.(fns{cell_no}).trial_std_threshold.("t" + trial),'r');
        yline(obj.data.placefieldData.(fns{cell_no}).trial_std_threshold.("t" + trial) / 3,'r--');
        yline(obj.data.placefieldData.(fns{cell_no}).trial_std_threshold.("t" + trial) * 2 / 3,'k--');
        yline(obj.data.placefieldData.(fns{cell_no}).trial_std_threshold.("t" + trial) * 4 / 3,'k--');
        yline(obj.data.placefieldData.(fns{cell_no}).trial_std_threshold.("t" + trial) * 5 / 3,'k--');
        yline(obj.data.placefieldData.(fns{cell_no}).trial_std_threshold.("t" + trial) * 6 / 3,'k--');
        annotation('textbox', [0.005 0.85 0.095 0.07], 'String', sprintf("Cell: %d\nTrial: %d\n",cell_no,trial));
        
    elseif(Args.AlphaAdSm)
        % Adaptive Smoothed Firing Rate Map
        cell_no = n;
        fns = fieldnames(obj.data.placefieldData);
        delete(findall(gcf,'type','annotation'));

        subplot(4,1,1);
        imagesc(obj.data.cellData.(fns{cell_no}).maps_raw);
        axis off
        colorbar;
        subplot(4,1,2);
        imagesc(obj.data.cellData.(fns{cell_no}).maps_adsm);
        axis off
        colorbar;
        subplot(4,1,3);
        plot(obj.data.placefieldData.(fns{cell_no}).basemapLrw, 'g'); legends{1} = sprintf('Raw');
        subplot(4,1,4);
        plot(obj.data.placefieldData.(fns{cell_no}).basemapLsm, 'b'); legends{2} = sprintf('AdSmooth');
        annotation('textbox', [0.005 0.85 0.095 0.07], 'String', sprintf("Alpha: %d\n", obj.data.cellData.(fns{cell_no}).alpha));    
        
    elseif(Args.FluorescenceTrace)
        
        % 95% SIC and GMM criteria
%         place_cell_sort = sortrows(obj.data.placefieldStats,[1 4]);
%         [c, ia, ic] = unique(place_cell_sort(:,1),'last');
%         place_cell_sort = sortrows(place_cell_sort(ia,:),[2 1]);
        
        % Only 95% SIC criteria
        place_cell = obj.data.isplacecell(obj.data.isplacecell(:,4) == 1,:);
        binFiringRate_trialAveraged = squeeze(mean(obj.data.binFiringRate(place_cell(:,1),:,:),2));
        [~,maxFiringRateBin] = max(binFiringRate_trialAveraged,[],2);
        place_cell = [place_cell maxFiringRateBin];
        place_cell_sort = sortrows(place_cell, 6);
        
        subplot(1,2,1)
        maps_all_combined = [];
        fns = fieldnames(obj.data.cellData);
        for cell_no = 1:length(fns)
            maps_all_combined = [maps_all_combined; obj.data.cellData.(fns{cell_no}).maps_adsm];
        end
        maps_placecells_combined = maps_all_combined(place_cell_sort(:,1),:);
        imagesc(maps_placecells_combined); title('AdSm'); xlabel('Bins'); ylabel('Neuron (sorted)');
        
        subplot(1,2,2)
        fluorescence_placecells_combined = obj.data.sessData.dF_F0_corrected(place_cell_sort(:,1),:);
        imagesc(fluorescence_placecells_combined); title('dF\_F0'); xlabel('Time (unit)'); ylabel('Neuron (sorted)');

    elseif(Args.TrialFluorescenceTrace)
        
        trial_no = n;
        trial_start_time = obj.data.sessData.data_trial(trial_no,2);
        trial_end_time = obj.data.sessData.data_trial(trial_no,3);
        [val,trial_start_idx] = min(abs(obj.data.sessData.tsF-trial_start_time));
        [val,trial_end_idx] = min(abs(obj.data.sessData.tsF-trial_end_time));
        
        place_cell_sort = sortrows(obj.data.placefieldStats,[1 4]);
        [c, ia, ic] = unique(place_cell_sort(:,1),'last');
        place_cell_sort = sortrows(place_cell_sort(ia,:),[2 1]);
        
        non_place_cell = setdiff(obj.data.placecellStats(:,1),place_cell_sort(:,1));
        
        maps_all_combined = [];
        fns = fieldnames(obj.data.cellData);
        for cell_no = 1:length(fns)
            maps_all_combined = [maps_all_combined; obj.data.cellData.(fns{cell_no}).maps_adsm];
        end
        maps_placecells_combined = maps_all_combined(place_cell_sort(:,1),:);
        maps_nonplacecells_combined = maps_all_combined(non_place_cell,:);
        
        subplot(1,4,1)
        imagesc(maps_placecells_combined); title('AdSm'); xlabel('Bins'); ylabel('Neuron (sorted)');
        
        subplot(1,4,2)
        imagesc(maps_nonplacecells_combined); title('AdSm'); xlabel('Bins'); ylabel('Neuron (sorted)');
        
        subplot(1,4,3)
        imagesc([maps_placecells_combined; maps_nonplacecells_combined]); title('AdSm'); xlabel('Bins'); ylabel('Neuron (sorted)');
        
        subplot(1,4,4)
        fluorescence_placecells_combined = obj.data.sessData.dF_F0_corrected(place_cell_sort(:,1),trial_start_idx:trial_end_idx);
        imagesc(fluorescence_placecells_combined); title('dF\_F0'); xlabel('Time (unit)'); ylabel('Neuron (sorted)');

    elseif(Args.PlaceFieldPositionEntropy)
        
        subplot(2,1,1)      
        y_data = histcounts(obj.data.placefieldStats(:,2),100);
        y_data = y_data ./ sum(y_data);
        bar([-49:50], [y_data(51:end) y_data(1:50)]);
                
        subplot(2,1,2)
        bounds = [floor(obj.data.placefieldStats(:,2) - obj.data.placefieldStats(:,3)) ceil(obj.data.placefieldStats(:,2) + obj.data.placefieldStats(:,3))];
        combined_bins = [];
        for i = 1:size(bounds,1)
            combined_bins = [combined_bins bounds(i,1):bounds(i,2)];
        end
        combined_bins(combined_bins < 1) = combined_bins(combined_bins < 1) + 100;
        combined_bins(combined_bins > 100) = combined_bins(combined_bins > 100) - 100;
        
        y_data = histcounts(combined_bins,100);
        y_data = y_data ./ sum(y_data);
        bar([-49:50], [y_data(51:end) y_data(1:50)]);
        
% 	else
% 		% code to plot yet another kind of plot
% 		
% 		% label the axis
% 		xlabel('X Axis')
% 		ylabel('Y Axis')

	end	

	% add an appropriate title
% 	sdstr = get(obj,'SessionDirs');

% 	title(getDataOrder('ShortName','DirString',sdstr{1}))
else
	% plot all data
    n = 1;
end

% The following code allows any commands to be executed as part of each plot
if(~isempty(Args.Cmds))
    % save the current figure in case Args.Cmds switches to another figure
    h = gcf;
    eval(Args.Cmds)
    % switch back to previous figure
    figure(h);
end

RR = eval('Args.ReturnVars');
lRR = length(RR);
if(lRR>0)
    for i=1:lRR
        RR1{i}=eval(RR{i});
    end 
    varargout = getReturnVal(Args.ReturnVars, RR1);
else
    varargout = {};
end
