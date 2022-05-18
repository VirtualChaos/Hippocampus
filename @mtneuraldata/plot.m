function [obj, varargout] = plot(obj,varargin)
%@dirfiles/plot Plot function for dirfiles object.
%   OBJ = plot(OBJ) creates plots for the neural response

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		  'ReturnVars',{''}, 'ArgsOnly',0, 'Cmds','', ...
          'FiringRateMapRaw',0, 'FiringRateMapAdSm',0, 'BinFiringRate',0, ...
          'TrialBinFiringRate', 0, 'Baseline',0, 'GMM',0, 'TrialGMM',0, 'TrialFits',0, 'AlphaAdSm',0);
Args.flags = {'LabelsOff','ArgsOnly','FiringRateMapRaw','FiringRateMapAdSm','BinFiringRate', ...
              'TrialBinFiringRate','Baseline','GMM', 'TrialGMM', 'TrialFits', 'AlphaAdSm'};
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
        
        subplot(2,1,1);
        % Adaptive Smoothed Firing Rate Map
        cell_no = n;
        fns = fieldnames(obj.data.placefieldData);
        
        imagesc(obj.data.cellData.(fns{cell_no}).binFiringRate); title('Trial Firing Rates'); xlabel('Position Bins'); ylabel('Trials');
        hold on
        %colorbar;
        hold on
        
        fns = fieldnames(obj.data.placefieldData);
        
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
                patch(x, y, 'w', 'FaceAlpha', '0.6', 'LineStyle', 'none'); hold on
            end
        end
        hold off
        
        % Gaussian Mixed Model
        cell_no = n;
        fns = fieldnames(obj.data.placefieldData);
        
        subplot(2,1,2);
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
                    y(101:126) = y(1:26);
                else
                    y(27:52) = y(127:152);
                end
                % Plot gaussian
                % plot(x(27:126), y(27:126), 'k-', 'LineWidth',3);
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
            plot(x(27:126), y(27:126), 'k-', 'LineWidth',3); legends{10} = sprintf('GMM'); hold on
            %legend(legends, 'Location', 'northeastoutside');
        end
        hold off
        
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
        plot(obj.data.placefieldData.(fns{cell_no}).baseline(27:126, 1), obj.data.placefieldData.(fns{cell_no}).baseline(27:126, 2),'r'); legends{2} = sprintf('Estimated Baseline');
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
                    y(101:126) = y(1:26);
                else
                    y(27:52) = y(127:152);
                end
                % Plot gaussian
                % plot(x(27:126), y(27:126), 'k-', 'LineWidth',3);
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
            plot(x(27:126), y(27:126), 'k-', 'LineWidth',3); legends{10} = sprintf('GMM'); hold on
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
%                     y(101:126) = y(1:26);
%                 else
%                     y(27:52) = y(127:152);
%                 end
%                 % Plot gaussian
%                 % plot(x(27:126), y(27:126), 'k-', 'LineWidth',3);
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
%             plot(x(27:126), y(27:126), 'k-', 'LineWidth',3); legends{10} = sprintf('GMM'); hold on
%             legend(legends, 'Location', 'northeastoutside');
%         end
%         hold off
        
    elseif(Args.TrialGMM)
        delete(findall(gcf,'type','annotation'));
        cell_no = floor((n-1)/obj.data.nTrials) + 1;
        trial = n - ((cell_no-1)*obj.data.nTrials);
        fns = fieldnames(obj.data.placefieldData);
        % fprintf("Cell: %d; Trial: %d\n",cell_no,trial);
        
        plot(obj.data.cellData.(fns{cell_no}).binFiringRate(trial,:), 'b'); title('Trial Firing Rates'); xlabel('Position Bins'); ylabel('Firing Rate (spikes/sec)');
        annotation('textbox', [0.005 0.85 0.095 0.07], 'String', sprintf("Cell: %d\nTrial: %d\n",cell_no,trial));
        
        if ~isempty(obj.data.placefieldData.(fns{cell_no}).trial_gmm_fits.("t" + num2str(trial)))
            mu_lines = [obj.data.placefieldData.(fns{cell_no}).trial_gmm_fits.("t" + num2str(trial))(:,1) [1:size(obj.data.placefieldData.(fns{cell_no}).trial_gmm_fits.("t" + num2str(trial)),1)]'];
            mu_lines_labels = cellstr("Field " + num2str(mu_lines(:,2)));
            xline(mu_lines(:,1), '--', mu_lines_labels, 'Color', '#D95319');
            ax = gca;
            for peak_no = 1:size(obj.data.placefieldData.(fns{cell_no}).trial_gmm_fits.("t" + num2str(trial)),1)
                left = obj.data.placefieldData.(fns{cell_no}).trial_gmm_fits.("t" + num2str(trial))(peak_no,1) - obj.data.placefieldData.(fns{cell_no}).trial_gmm_fits.("t" + num2str(trial))(peak_no,2);
                right = obj.data.placefieldData.(fns{cell_no}).trial_gmm_fits.("t" + num2str(trial))(peak_no,1) + obj.data.placefieldData.(fns{cell_no}).trial_gmm_fits.("t" + num2str(trial))(peak_no,2);
                x = [left right right left];
                y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];
                patch(x, y, 'k', 'FaceAlpha', '0.1', 'LineStyle', 'none');
            end
            
            gaus = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;
            
            for i = 1:size(obj.data.placefieldData.(fns{cell_no}).trial_gmm_fits.("t" + num2str(trial)),1)
                x = linspace(-25,126,152);
                mu = obj.data.placefieldData.(fns{cell_no}).trial_gmm_fits.("t" + num2str(trial))(i,1);
                sig = obj.data.placefieldData.(fns{cell_no}).trial_gmm_fits.("t" + num2str(trial))(i,2);
                amp = obj.data.placefieldData.(fns{cell_no}).trial_gmm_fits.("t" + num2str(trial))(i,3);
                vo = 0;
                y = gaus(x,mu,sig,amp,vo);
                if mu < 50
                    y(101:126) = y(1:26);
                else
                    y(27:52) = y(127:152);
                end
                % Plot gaussian
                % plot(x(27:126), y(27:126), 'k-', 'LineWidth',3);
            end
            axis([1 100 -inf inf]);
%             yline((baseline_mean_height + peak_height_from_baseline/2),'r');
%             yline(std_threshold,'b');
            
            x = linspace(-25,126,152);           
            y = 0;
            for i = 1:size(obj.data.placefieldData.(fns{cell_no}).trial_gmm_fits.("t" + num2str(trial)),1)
                y = gaus(x,obj.data.placefieldData.(fns{cell_no}).trial_gmm_fits.("t" + num2str(trial))(i,1),obj.data.placefieldData.(fns{cell_no}).trial_gmm_fits.("t" + num2str(trial))(i,2),obj.data.placefieldData.(fns{cell_no}).trial_gmm_fits.("t" + num2str(trial))(i,3),vo) + y;
            end
            
%             legends{4} = ''; legends{5} = ''; legends{6} = ''; legends{7} = ''; legends{8} = ''; legends{9} = '';
%             plot(x(27:126), y(27:126), 'k-', 'LineWidth',3); legends{10} = sprintf('GMM');
%             legend(legends, 'Location', 'northeastoutside');
        end
        hold off
        
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
