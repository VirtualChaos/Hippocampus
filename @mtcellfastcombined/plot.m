function [obj, varargout] = plot(obj,varargin)
%@dirfiles/plot Plot function for dirfiles object.
%   OBJ = plot(OBJ) creates plots for the neural response

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		  'ReturnVars',{''}, 'ArgsOnly',0, 'Cmds','', ...
          'FiringRateMapRaw',0, 'FiringRateMapAdSm',0, 'BinFiringRate',0, ...
          'TrialBinFiringRate', 0, 'Baseline',0, 'GMM',0, 'TrialFits',0, 'AlphaAdSm',0, ...
          'FluorescenceTrace',0, 'TrialFluorescenceTrace',0, 'PlaceFieldPositionEntropy',0, 'PopVec',0, ...
          'RestRunCorrDist',0,'Links',0,'Corr',0);
Args.flags = {'LabelsOff','ArgsOnly','FiringRateMapRaw','FiringRateMapAdSm','BinFiringRate', ...
              'TrialBinFiringRate','Baseline','GMM', 'TrialFits', 'AlphaAdSm', ...
              'FluorescenceTrace','TrialFluorescenceTrace','PlaceFieldPositionEntropy','PopVec', ...
              'RestRunCorrDist','Links','Corr'};
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
        
    elseif(Args.PopVec)
        
        bin_no = n;
        FiringRate = squeeze(obj.data.binFiringRate(:,:,bin_no));
        
    elseif(Args.RestRunCorrDist)
        
        Cor_now     = obj.data.Cor_Mspks.correlations;
        k = 2; % Green (Rest)
        MeanCor(k)    = mean(Cor_now(k,:));
        MedianCor(k)  = median(Cor_now(k,:));
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
        Leftfit_sym  = a1.*exp(-((xLeft_sym-b1)./c1L).^2);
        Rightfit_sym = a1.*exp(-((xRight_sym-b1)./c1R).^2);
        figure(1);
%         plot(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'Color',[0.4660 0.6740 0.1880]);hold on;
%         plot(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'Color',[0.4660 0.6740 0.1880]);hold on;
        area(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'EdgeColor',[0.4660 0.6740 0.1880],'FaceColor',[0.4660 0.6740 0.1880], ...
            'FaceAlpha',0.5,'LineWidth',2);
        hold on
        area(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'EdgeColor',[0.4660 0.6740 0.1880],'FaceColor',[0.4660 0.6740 0.1880], ...
            'FaceAlpha',0.5,'LineWidth',2);
        %     plot(X,N,'k');xlim([-0.1,0.1])
        hold on
        k = 3; % Yellow (Run)
        MeanCor(k)    = mean(Cor_now(k,:));
        MedianCor(k)  = median(Cor_now(k,:));
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
        Leftfit_sym  = a1.*exp(-((xLeft_sym-b1)./c1L).^2);
        Rightfit_sym = a1.*exp(-((xRight_sym-b1)./c1R).^2);
        figure(1);
%         plot(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'Color',[0.9290 0.6940 0.1250]);hold on;
%         plot(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'Color',[0.9290 0.6940 0.1250]);hold on;
        area(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak), ...
            'EdgeColor',[0.9290 0.6940 0.1250],'FaceColor',[0.9290 0.6940 0.1250],'FaceAlpha',0.5,'LineWidth',2);
        hold on
        area(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end), ...
            'EdgeColor',[0.9290 0.6940 0.1250],'FaceColor',[0.9290 0.6940 0.1250],'FaceAlpha',0.5,'LineWidth',2);
        %     plot(X,N,'k');xlim([-0.1,0.1])
        xlim([-0.1,0.1])
        xlabel("Correlation Coefficient",'FontSize',18); ylabel("Normalised density",'FontSize',18);
        hold off
        
    elseif(Args.Links)
        
        loop_pc_non_pc = n;
        
        switch loop_pc_non_pc
            case 1 % All cells
                rest_proportion_pos_neg = obj.data.Cor_Mspks.links.all.rest_proportion_pos_neg;
                run_proportion_pos_neg = obj.data.Cor_Mspks.links.all.run_proportion_pos_neg;
                k_values = obj.data.Cor_Mspks.links.all.k_values;
            case 2 % PC-PC
                rest_proportion_pos_neg = obj.data.Cor_Mspks.links.pc_pc.rest_proportion_pos_neg;
                run_proportion_pos_neg = obj.data.Cor_Mspks.links.pc_pc.run_proportion_pos_neg;
                k_values = obj.data.Cor_Mspks.links.pc_pc.k_values;
            case 3 % non PC-non PC
                rest_proportion_pos_neg = obj.data.Cor_Mspks.links.nonpc_nonpc.rest_proportion_pos_neg;
                run_proportion_pos_neg = obj.data.Cor_Mspks.links.nonpc_nonpc.run_proportion_pos_neg;
                k_values = obj.data.Cor_Mspks.links.nonpc_nonpc.k_values;
            case 4 % PC-non PC
                rest_proportion_pos_neg = obj.data.Cor_Mspks.links.pc_nonpc.rest_proportion_pos_neg;
                run_proportion_pos_neg = obj.data.Cor_Mspks.links.pc_nonpc.run_proportion_pos_neg;
                k_values = obj.data.Cor_Mspks.links.pc_nonpc.k_values;
            case 5 % non PC-PC
                rest_proportion_pos_neg = obj.data.Cor_Mspks.links.nonpc_pc.rest_proportion_pos_neg;
                run_proportion_pos_neg = obj.data.Cor_Mspks.links.nonpc_pc.run_proportion_pos_neg;
                k_values = obj.data.Cor_Mspks.links.nonpc_pc.k_values;
        end
                
        for loop_links = 1:6
            subplot(3,2,loop_links)
            if loop_links == 1 | loop_links == 3 | loop_links == 5
                proportion_pos_neg = rest_proportion_pos_neg;
            elseif loop_links == 2 | loop_links == 4 | loop_links == 6
                proportion_pos_neg = run_proportion_pos_neg;
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
            % [N,edges] = histcounts(link_dist(:,2),0:5:max(link_dist(:,2)+5),'Normalization','probability'); %1:max(link_dist(:,2))+1
            degree_dist = [[0:0.2:max(link_dist(:,2))]' N'];
            
            bar(edges(1:end-1),degree_dist(:,2));
            xline(k_mean,'k--');
            xline(k_median,'r--');
%             loglog(degree_dist(:,2),'x');
            if loop_links == 1
                title("Resting (All links)")
            elseif loop_links == 2
                title("Running (All links)")
            elseif loop_links == 3
                title("Resting (Positive links)")
            elseif loop_links == 4
                title("Running (Positive links)")
            elseif loop_links == 5
                title("Resting (Negative links)")
            elseif loop_links == 6
                title("Running (Negative links)")
            end
        end
                   
    elseif(Args.Corr)
        
        % Place-Non Place cell connectivity
        temp_cell_list = 1:obj.data.nNeuron;
        pc_list = obj.data.PC_ID;
        non_pc_list = temp_cell_list(~ismember(temp_cell_list,pc_list));
        
        % SIC list
        SIC_list = [[1:size(obj.data.SIC.SICVec,2)]' obj.data.SIC.SICVec'];
        SIC_list = sortrows(SIC_list,2);
        
        cr_res = obj.data.Cor_Mspks.cr_res;
        
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
        
        subplot(3,2,1)
        Cor_now = reshape(pc_cr_res.',1,[]); % Red
        [MeanCor,MedianCor,FWHM,IDpeak,xLeft_sym,Leftfit_sym,xRight_sym,Rightfit_sym] = FWHM_compute(Cor_now);
        %         plot(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'Color',[0.4660 0.6740 0.1880]);hold on;
        %         plot(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'Color',[0.4660 0.6740 0.1880]);hold on;
        area(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'EdgeColor',[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840], ...
            'FaceAlpha',0.5,'LineWidth',2);
        hold on
        area(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'EdgeColor',[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840], ...
            'FaceAlpha',0.5,'LineWidth',2);
        hold on
        xline(MedianCor,'Color',[0.6350 0.0780 0.1840],'LineStyle','--');
        %     plot(X,N,'k');xlim([-0.1,0.1])
        hold on
        Cor_now = reshape(nonpc_cr_res.',1,[]); % Blue
        [MeanCor,MedianCor,FWHM,IDpeak,xLeft_sym,Leftfit_sym,xRight_sym,Rightfit_sym] = FWHM_compute(Cor_now); % Blue
        %         plot(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'Color',[0.9290 0.6940 0.1250]);hold on;
        %         plot(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'Color',[0.9290 0.6940 0.1250]);hold on;
        area(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak), ...
            'EdgeColor',[0.3010 0.7450 0.9330],'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.5,'LineWidth',2);
        hold on
        area(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end), ...
            'EdgeColor',[0.3010 0.7450 0.9330],'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.5,'LineWidth',2);
        %     plot(X,N,'k');xlim([-0.1,0.1])
        hold on
        xline(MedianCor,'Color',[0.3010 0.7450 0.9330],'LineStyle','--');
        xlim([-0.1,0.1])
        xlabel("Correlation Coefficient",'FontSize',18); ylabel("Normalised density",'FontSize',18);
        title("PC-Non PC")
        hold off
        
        subplot(3,2,2)
        Cor_now = reshape(high_SIC_cr_res.',1,[]); % Red
        [MeanCor,MedianCor,FWHM,IDpeak,xLeft_sym,Leftfit_sym,xRight_sym,Rightfit_sym] = FWHM_compute(Cor_now);
        %         plot(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'Color',[0.4660 0.6740 0.1880]);hold on;
        %         plot(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'Color',[0.4660 0.6740 0.1880]);hold on;
        area(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'EdgeColor',[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840], ...
            'FaceAlpha',0.5,'LineWidth',2);
        hold on
        area(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'EdgeColor',[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840], ...
            'FaceAlpha',0.5,'LineWidth',2);
        hold on
        xline(MedianCor,'Color',[0.6350 0.0780 0.1840],'LineStyle','--');
        %     plot(X,N,'k');xlim([-0.1,0.1])
        hold on
        Cor_now = reshape(low_SIC_cr_res.',1,[]); % Blue
        [MeanCor,MedianCor,FWHM,IDpeak,xLeft_sym,Leftfit_sym,xRight_sym,Rightfit_sym] = FWHM_compute(Cor_now); % Blue
        %         plot(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'Color',[0.9290 0.6940 0.1250]);hold on;
        %         plot(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'Color',[0.9290 0.6940 0.1250]);hold on;
        area(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak), ...
            'EdgeColor',[0.3010 0.7450 0.9330],'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.5,'LineWidth',2);
        hold on
        area(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end), ...
            'EdgeColor',[0.3010 0.7450 0.9330],'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.5,'LineWidth',2);
        %     plot(X,N,'k');xlim([-0.1,0.1])
        hold on
        xline(MedianCor,'Color',[0.3010 0.7450 0.9330],'LineStyle','--');
        xlim([-0.1,0.1])
        xlabel("Correlation Coefficient",'FontSize',18); ylabel("Normalised density",'FontSize',18);
        title("SIC")
        hold off
        
        subplot(3,2,3)
        Cor_now = reshape(pc_pc_cr_res.',1,[]); % Red
        [MeanCor,MedianCor,FWHM,IDpeak,xLeft_sym,Leftfit_sym,xRight_sym,Rightfit_sym] = FWHM_compute(Cor_now);
        %         plot(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'Color',[0.4660 0.6740 0.1880]);hold on;
        %         plot(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'Color',[0.4660 0.6740 0.1880]);hold on;
        area(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'EdgeColor',[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840], ...
            'FaceAlpha',0.5,'LineWidth',2);
        hold on
        area(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'EdgeColor',[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840], ...
            'FaceAlpha',0.5,'LineWidth',2);
        hold on
        xline(MedianCor,'Color',[0.6350 0.0780 0.1840],'LineStyle','--');
        %     plot(X,N,'k');xlim([-0.1,0.1])
        hold on
        Cor_now = reshape(pc_nonpc_cr_res.',1,[]); % Blue
        [MeanCor,MedianCor,FWHM,IDpeak,xLeft_sym,Leftfit_sym,xRight_sym,Rightfit_sym] = FWHM_compute(Cor_now); % Blue
        %         plot(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'Color',[0.9290 0.6940 0.1250]);hold on;
        %         plot(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'Color',[0.9290 0.6940 0.1250]);hold on;
        area(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak), ...
            'EdgeColor',[0.3010 0.7450 0.9330],'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.5,'LineWidth',2);
        hold on
        area(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end), ...
            'EdgeColor',[0.3010 0.7450 0.9330],'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.5,'LineWidth',2);
        %     plot(X,N,'k');xlim([-0.1,0.1])
        hold on
        xline(MedianCor,'Color',[0.3010 0.7450 0.9330],'LineStyle','--');
        xlim([-0.1,0.1])
        xlabel("Correlation Coefficient",'FontSize',18); ylabel("Normalised density",'FontSize',18);
        title("PC")
        hold off
        
        subplot(3,2,4)
        Cor_now = reshape(high_high_SIC_cr_res.',1,[]); % Red
        [MeanCor,MedianCor,FWHM,IDpeak,xLeft_sym,Leftfit_sym,xRight_sym,Rightfit_sym] = FWHM_compute(Cor_now);
%         plot(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'Color',[0.4660 0.6740 0.1880]);hold on;
%         plot(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'Color',[0.4660 0.6740 0.1880]);hold on;
        area(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'EdgeColor',[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840], ...
            'FaceAlpha',0.5,'LineWidth',2);
        hold on
        area(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'EdgeColor',[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840], ...
            'FaceAlpha',0.5,'LineWidth',2);
        hold on
        xline(MedianCor,'Color',[0.6350 0.0780 0.1840],'LineStyle','--');
        %     plot(X,N,'k');xlim([-0.1,0.1])
        hold on
        Cor_now = reshape(high_low_SIC_cr_res.',1,[]); % Blue
        [MeanCor,MedianCor,FWHM,IDpeak,xLeft_sym,Leftfit_sym,xRight_sym,Rightfit_sym] = FWHM_compute(Cor_now); % Blue
        %         plot(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'Color',[0.9290 0.6940 0.1250]);hold on;
        %         plot(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'Color',[0.9290 0.6940 0.1250]);hold on;
        area(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak), ...
            'EdgeColor',[0.3010 0.7450 0.9330],'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.5,'LineWidth',2);
        hold on
        area(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end), ...
            'EdgeColor',[0.3010 0.7450 0.9330],'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.5,'LineWidth',2);
        %     plot(X,N,'k');xlim([-0.1,0.1])
        hold on
        xline(MedianCor,'Color',[0.3010 0.7450 0.9330],'LineStyle','--');
        xlim([-0.1,0.1])
        xlabel("Correlation Coefficient",'FontSize',18); ylabel("Normalised density",'FontSize',18);
        title("High SIC")
        hold off
        
        subplot(3,2,5)
        Cor_now = reshape(nonpc_nonpc_cr_res.',1,[]); % Red
        [MeanCor,MedianCor,FWHM,IDpeak,xLeft_sym,Leftfit_sym,xRight_sym,Rightfit_sym] = FWHM_compute(Cor_now);
        %         plot(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'Color',[0.4660 0.6740 0.1880]);hold on;
        %         plot(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'Color',[0.4660 0.6740 0.1880]);hold on;
        area(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'EdgeColor',[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840], ...
            'FaceAlpha',0.5,'LineWidth',2);
        hold on
        area(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'EdgeColor',[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840], ...
            'FaceAlpha',0.5,'LineWidth',2);
        hold on
        xline(MedianCor,'Color',[0.6350 0.0780 0.1840],'LineStyle','--');
        %     plot(X,N,'k');xlim([-0.1,0.1])
        hold on
        Cor_now = reshape(nonpc_pc_cr_res.',1,[]); % Blue
        [MeanCor,MedianCor,FWHM,IDpeak,xLeft_sym,Leftfit_sym,xRight_sym,Rightfit_sym] = FWHM_compute(Cor_now); % Blue
        %         plot(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'Color',[0.9290 0.6940 0.1250]);hold on;
        %         plot(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'Color',[0.9290 0.6940 0.1250]);hold on;
        area(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak), ...
            'EdgeColor',[0.3010 0.7450 0.9330],'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.5,'LineWidth',2);
        hold on
        area(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end), ...
            'EdgeColor',[0.3010 0.7450 0.9330],'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.5,'LineWidth',2);
        %     plot(X,N,'k');xlim([-0.1,0.1])
        hold on
        xline(MedianCor,'Color',[0.3010 0.7450 0.9330],'LineStyle','--');
        xlim([-0.1,0.1])
        xlabel("Correlation Coefficient",'FontSize',18); ylabel("Normalised density",'FontSize',18);
        title("Non PC")
        hold off
        
        subplot(3,2,6)
        Cor_now = reshape(low_low_SIC_cr_res.',1,[]); % Red
        [MeanCor,MedianCor,FWHM,IDpeak,xLeft_sym,Leftfit_sym,xRight_sym,Rightfit_sym] = FWHM_compute(Cor_now);
        %         plot(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'Color',[0.4660 0.6740 0.1880]);hold on;
        %         plot(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'Color',[0.4660 0.6740 0.1880]);hold on;
        area(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'EdgeColor',[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840], ...
            'FaceAlpha',0.5,'LineWidth',2);
        hold on
        area(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'EdgeColor',[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840], ...
            'FaceAlpha',0.5,'LineWidth',2);
        hold on
        xline(MedianCor,'Color',[0.6350 0.0780 0.1840],'LineStyle','--');
%         plot(X,N,'k');xlim([-0.1,0.1])
        hold on
        Cor_now = reshape(low_high_SIC_cr_res.',1,[]); % Blue
        [MeanCor,MedianCor,FWHM,IDpeak,xLeft_sym,Leftfit_sym,xRight_sym,Rightfit_sym] = FWHM_compute(Cor_now); % Blue
        %         plot(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak),'Color',[0.9290 0.6940 0.1250]);hold on;
        %         plot(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end),'Color',[0.9290 0.6940 0.1250]);hold on;
        area(xLeft_sym(1:IDpeak),Leftfit_sym(1:IDpeak), ...
            'EdgeColor',[0.3010 0.7450 0.9330],'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.5,'LineWidth',2);
        hold on
        area(xRight_sym(round(length(xRight_sym)/2):end),Rightfit_sym(round(length(xRight_sym)/2):end), ...
            'EdgeColor',[0.3010 0.7450 0.9330],'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.5,'LineWidth',2);
        %     plot(X,N,'k');xlim([-0.1,0.1])
        hold on
        xline(MedianCor,'Color',[0.3010 0.7450 0.9330],'LineStyle','--');
        xlim([-0.1,0.1])
        xlabel("Correlation Coefficient",'FontSize',18); ylabel("Normalised density",'FontSize',18);
        title("Low SIC")
        hold off

        
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
