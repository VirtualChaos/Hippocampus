function [obj, varargout] = plot(obj,varargin)
%@dirfiles/plot Plot function for dirfiles object.
%   OBJ = plot(OBJ) creates plots for the individual mice

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		  'ReturnVars',{''}, 'ArgsOnly',0, 'Cmds','', ...
          'MatchedFiringRateMap',0, 'MatchedCorrMatrix',0,'PlaceFieldPositionEntropy',0,...
          'AdSm',0,'Vel',0,'VelDist',0,'Lick',0,'LickBurstWidth',0,'LickRZEstimate',0);
Args.flags = {'LabelsOff','ArgsOnly','MatchedFiringRateMap','MatchedCorrMatrix','PlaceFieldPositionEntropy',...
                'AdSm','Vel','VelDist','Lick','LickBurstWidth','LickRZEstimate'};
[Args,varargin2] = getOptArgs(varargin,Args);

% if user select 'ArgsOnly', return only Args structure for an empty object
if Args.ArgsOnly
    Args = rmfield (Args, 'ArgsOnly');
    varargout{1} = {'Args',Args};
    return;
end

if(~isempty(Args.NumericArguments))
	n = Args.NumericArguments{1};
	if(Args.MatchedFiringRateMap)
        
        sess_list = fieldnames(obj.data.neuralCombined);
        for sess_no = 1:size(obj.data.roiMatchData.allSessionMapping,2)
            
            subplot(2,size(obj.data.roiMatchData.allSessionMapping,2),sess_no);
            cell_no = obj.data.roiMatchData.allSessionMapping(n,sess_no);
            neuralData = obj.data.neuralCombined.(sess_list{sess_no});
            fns = fieldnames(neuralData.placefieldData);
            
            % Adaptive Smoothed Firing Rate Map
            imagesc(neuralData.cellData.(fns{cell_no}).binFiringRate); title("Trial Firing Rates: Cell " + num2str(cell_no)); xlabel('Position Bins'); ylabel('Trials');
            hold on
            %colorbar;
            hold on
            
            if ~isempty(neuralData.placefieldData.(fns{cell_no}).GMM) & neuralData.isplacecell(cell_no,4)
                mu_lines = [neuralData.placefieldData.(fns{cell_no}).GMM(:,1) [1:size(neuralData.placefieldData.(fns{cell_no}).GMM,1)]'];
                mu_lines_labels = cellstr("Field " + num2str(mu_lines(:,2)));
                xline(mu_lines(:,1), '--', mu_lines_labels, 'Color', '#D95319'); hold on
                ax = gca;
                for peak_no = 1:size(neuralData.placefieldData.(fns{cell_no}).GMM,1)
                    left = neuralData.placefieldData.(fns{cell_no}).GMM(peak_no,1) - neuralData.placefieldData.(fns{cell_no}).GMM(peak_no,2);
                    right = neuralData.placefieldData.(fns{cell_no}).GMM(peak_no,1) + neuralData.placefieldData.(fns{cell_no}).GMM(peak_no,2);
                    x = [left right right left];
                    y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];
                    patch(x, y, 'w', 'FaceAlpha', '0.6', 'LineStyle', 'none'); hold on
                end
            end
            hold off
            
            % Gaussian Mixed Model
            subplot(2,size(obj.data.roiMatchData.allSessionMapping,2),sess_no+size(obj.data.roiMatchData.allSessionMapping,2));
            title('Firing Rate Map'); xlabel('Position Bins'); ylabel('Firing Rate (spike/sec)');
            
            plot(neuralData.placefieldData.(fns{cell_no}).basemapLrw, 'g'); legends{1} = sprintf('Raw'); hold on
            
            plot(neuralData.placefieldData.(fns{cell_no}).basemapLsm, 'b'); legends{2} = sprintf('AdSmooth'); hold on
            
            plot(imgaussfilt(neuralData.placefieldData.(fns{cell_no}).basemapLsm, 3), 'c'); legends{2} = sprintf('AdSmooth'); hold on
            
            plot(neuralData.placefieldData.(fns{cell_no}).mapLsm_FFT_Low_Pass, 'r'); legends{3} = sprintf('AdSm-FFT-Low-Pass'); hold on
            
            axis([1 100 -inf inf]);
            yline((neuralData.placefieldData.(fns{cell_no}).baseline_mean_height + neuralData.placefieldData.(fns{cell_no}).peak_height_from_baseline/2),'r'); hold on
            yline(neuralData.placefieldData.(fns{cell_no}).std_threshold,'b'); hold on
            yline(neuralData.placefieldData.(fns{cell_no}).std_threshold / neuralData.placefieldData.(fns{cell_no}).Args.PeakThreshold,'b--'); hold on
            
            if ~isempty(neuralData.placefieldData.(fns{cell_no}).GMM) & neuralData.isplacecell(cell_no,4)
                mu_lines = [neuralData.placefieldData.(fns{cell_no}).GMM(:,1) [1:size(neuralData.placefieldData.(fns{cell_no}).GMM,1)]'];
                mu_lines_labels = cellstr("Field " + num2str(mu_lines(:,2)));
                xline(mu_lines(:,1), '--', mu_lines_labels, 'Color', '#D95319'); hold on
                ax = gca;
                for peak_no = 1:size(neuralData.placefieldData.(fns{cell_no}).GMM,1)
                    left = neuralData.placefieldData.(fns{cell_no}).GMM(peak_no,1) - neuralData.placefieldData.(fns{cell_no}).GMM(peak_no,2);
                    right = neuralData.placefieldData.(fns{cell_no}).GMM(peak_no,1) + neuralData.placefieldData.(fns{cell_no}).GMM(peak_no,2);
                    x = [left right right left];
                    y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];
                    patch(x, y, 'k', 'FaceAlpha', '0.1', 'LineStyle', 'none'); hold on
                end
                
                gaus = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;
                
                for i = 1:size(neuralData.placefieldData.(fns{cell_no}).GMM,1)
                    x = linspace(-25,126,152);
                    mu = neuralData.placefieldData.(fns{cell_no}).GMM(i,1);
                    sig = neuralData.placefieldData.(fns{cell_no}).GMM(i,2);
                    amp = neuralData.placefieldData.(fns{cell_no}).GMM(i,3);
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
                
                x = linspace(-25,126,152);
                y = 0;
                for i = 1:size(neuralData.placefieldData.(fns{cell_no}).GMM,1)
                    y = gaus(x,neuralData.placefieldData.(fns{cell_no}).GMM(i,1),neuralData.placefieldData.(fns{cell_no}).GMM(i,2),neuralData.placefieldData.(fns{cell_no}).GMM(i,3),vo) + y;
                end
                
                legends{4} = ''; legends{5} = ''; legends{6} = ''; legends{7} = ''; legends{8} = ''; legends{9} = '';
                plot(x(27:126), y(27:126), 'k-', 'LineWidth',3); legends{10} = sprintf('GMM'); hold on
                %legend(legends, 'Location', 'northeastoutside');
            end
            hold off
        end
        
    elseif(Args.MatchedCorrMatrix)
        sess_list = fieldnames(obj.data.neuralCombined);
        rawFiringMap = [];
        adsmFiringMap = [];
        lowpassFiringMap = [];
        for sess_no = 1:size(obj.data.roiMatchData.allSessionMapping,2)
            cell_no = obj.data.roiMatchData.allSessionMapping(n,sess_no);
            neuralData = obj.data.neuralCombined.(sess_list{sess_no});
            fns = fieldnames(neuralData.placefieldData);
            
            rawFiringMap = [rawFiringMap; neuralData.placefieldData.(fns{cell_no}).basemapLrw];
            adsmFiringMap = [adsmFiringMap; neuralData.placefieldData.(fns{cell_no}).basemapLsm];
            lowpassFiringMap = [lowpassFiringMap; neuralData.placefieldData.(fns{cell_no}).mapLsm_FFT_Low_Pass];
        end
        
        rawCorrMatrix = zeros(size(obj.data.roiMatchData.allSessionMapping,2));
        adsmCorrMatrix = zeros(size(obj.data.roiMatchData.allSessionMapping,2));
        lowpassCorrMatrix = zeros(size(obj.data.roiMatchData.allSessionMapping,2));
        for i = 1:size(obj.data.roiMatchData.allSessionMapping,2)
            for j = 1:size(obj.data.roiMatchData.allSessionMapping,2)
                a = rawFiringMap(i,:);
                b = rawFiringMap(j,:);
                rawCorrMatrix(i,j) = sum(a.*b,'omitnan') / sqrt(sum(a.^2,'omitnan')*sum(b.^2,'omitnan'));
                a = adsmFiringMap(i,:);
                b = adsmFiringMap(j,:);
                adsmCorrMatrix(i,j) = sum(a.*b,'omitnan') / sqrt(sum(a.^2,'omitnan')*sum(b.^2,'omitnan'));
                a = lowpassFiringMap(i,:);
                b = lowpassFiringMap(j,:);
                lowpassCorrMatrix(i,j) = sum(a.*b,'omitnan') / sqrt(sum(a.^2,'omitnan')*sum(b.^2,'omitnan'));
            end
        end
        
        colormap hot
        subplot(3,1,1)
        imagesc(rawCorrMatrix);; title("Raw Firing Map Correlation Matrix"); colorbar;
        subplot(3,1,2)
        imagesc(adsmCorrMatrix);; title("Adaptive Smoothed Firing Map Correlation Matrix"); colorbar;
        subplot(3,1,3)
        imagesc(lowpassCorrMatrix);; title("Low Pass Firing Map Correlation Matrix"); colorbar;
        
    elseif(Args.PlaceFieldPositionEntropy)
        sess_list = fieldnames(obj.data.neuralCombined);
        for sess_no = 1:size(sess_list,1)
            subplot(2,ceil(size(sess_list,1)/2),sess_no);
            title("Place field positions"); xlabel('Position Bins'); ylabel('Count');
            neuralData = obj.data.neuralCombined.(sess_list{sess_no});
            
            bounds = [floor(neuralData.placefieldStats(:,2) - neuralData.placefieldStats(:,3)) ceil(neuralData.placefieldStats(:,2) + neuralData.placefieldStats(:,3))];
            combined_bins = [];
            for i = 1:size(bounds,1)
               combined_bins = [combined_bins bounds(i,1):bounds(i,2)];
            end
            combined_bins(combined_bins < 1) = combined_bins(combined_bins < 1) + 100;
            combined_bins(combined_bins > 100) = combined_bins(combined_bins > 100) - 100;
            
            y_data = histcounts(combined_bins,100);
            y_data = y_data ./ sum(y_data);
            bar([-49:50], [y_data(51:end) y_data(1:50)]);            
        end
        
    elseif(Args.AdSm)
        sess_list = fieldnames(obj.data.neuralCombined);
        for sess_no = 1:size(sess_list,1)
            subplot(1,size(sess_list,1),sess_no);
            neuralData = obj.data.neuralCombined.(sess_list{sess_no});
            
            place_cell_sort = sortrows(neuralData.placefieldStats,[1 4]);
            [c, ia, ic] = unique(place_cell_sort(:,1),'last');
            place_cell_sort = sortrows(place_cell_sort(ia,:),[2 1]);
            
            non_place_cell = setdiff(neuralData.placecellStats(:,1),place_cell_sort(:,1));
            
            maps_all_combined = [];
            fns = fieldnames(neuralData.cellData);
            for cell_no = 1:length(fns)
                maps_all_combined = [maps_all_combined; neuralData.cellData.(fns{cell_no}).maps_adsm];
            end
            maps_placecells_combined = maps_all_combined(place_cell_sort(:,1),:);
            maps_nonplacecells_combined = maps_all_combined(non_place_cell,:);
            %imagesc([maps_placecells_combined; maps_nonplacecells_combined]); title('AdSm'); xlabel('Bins'); ylabel('Neuron (sorted)');
            %imagesc(maps_all_combined); title('AdSm'); xlabel('Bins'); ylabel('Neuron (sorted)');
            
            [M,I] = max(maps_all_combined,[],2,'omitnan');
            temp_idx = sortrows([[1:size(I,1)]' I],[2 1]);
            imagesc(maps_all_combined(temp_idx(:,1),:)); title('AdSm'); xlabel('Bins'); ylabel('Neuron (sorted)');
        end
        
    elseif(Args.Vel)
        sess_list = fieldnames(obj.data.sessionCombined);
                
        mean_velocity = zeros(length(fieldnames(obj.data.sessionCombined)),1);
        mean_track_velocity = zeros(length(fieldnames(obj.data.sessionCombined)),1);
        for sess_no = 1:length(fieldnames(obj.data.sessionCombined))
            sessionData = obj.data.sessionCombined.(sess_list{sess_no});
            mean_velocity(sess_no) = mean(sessionData.velocity_averaged(:,4), 'omitnan');
            mean_track_velocity(sess_no) = mean(sessionData.velocity_binned, 'omitnan');
        end
        
        subplot(2,1,1)
        plot(obj.data.sessionDays, mean_velocity); title('Mean Velocity'); xlabel('Sessions (Day)'); ylabel('Velocity (cm/s)');
        
        subplot(2,1,2)
        plot(obj.data.sessionDays, mean_track_velocity); title('Mean Velocity (Binned)'); xlabel('Sessions (Day)'); ylabel('Velocity (cm/s)');
        
    elseif(Args.VelDist)
        sess_list = fieldnames(obj.data.sessionCombined);
                
        for sess_no = 1:length(fieldnames(obj.data.sessionCombined))
            subplot(ceil(length(fieldnames(obj.data.sessionCombined))/4),4,sess_no)
            sessionData = obj.data.sessionCombined.(sess_list{sess_no});
%             edges = 0:0.5:ceil(max(sessionData.velocity_averaged(:,4)));
%             histogram(sessionData.velocity_averaged(sessionData.velocity_averaged(:,4) ~= 0,4),edges); title("D" + obj.data.sessionDays(sess_no)); ylabel('Velocity (cm/s)');
            
            edges = 0:0.5:ceil(max(sessionData.velocity_binned));
            histogram(sessionData.velocity_binned,edges); title("D" + obj.data.sessionDays(sess_no)); ylabel('Velocity (cm/s)');
        end
        
    elseif(Args.Lick)
        sess_no = n;
        sess_list = fieldnames(obj.data.sessionCombined);
        sessionData = obj.data.sessionCombined.(sess_list{sess_no});
        
        subplot(2,1,1);
        plot(sessionData.session_data_exclude_zero_trials(:,1), sessionData.session_data_exclude_zero_trials(:,5)*2, 'b'); xlabel('Time (s)');
        hold on
        plot(sessionData.session_data_exclude_zero_trials(:,1), sessionData.session_data_exclude_zero_trials(:,6), 'r'); title('Licks'); xlabel('Time (s)');
        hold off
        
        subplot(2,1,2);
        yyaxis left
        plot(sessionData.session_data_exclude_zero_trials(:,1), sessionData.session_data_exclude_zero_trials(:,7)*100, 'b'); title('Lick Bursts'); xlabel('Time (s)'); ylabel('Distance (Bin)'); xline(sessionData.data_trial(:,3), 'r');
        hold on
        plot(sessionData.session_data_exclude_zero_trials(:,1), sessionData.lick_burst(:)*80, 'k-');
        yyaxis right
        plot(sessionData.session_data_exclude_zero_trials(:,1),sessionData.velocity_averaged(:,4)); title('Velocity'); xlabel('Time (s)'); ylabel('Velocity (cm/s)');
        hold off
         
    elseif(Args.LickBurstWidth)
        sess_list = fieldnames(obj.data.sessionCombined);
                
        for sess_no = 1:length(fieldnames(obj.data.sessionCombined))
            subplot(ceil(length(fieldnames(obj.data.sessionCombined))/4),4,sess_no)
            sessionData = obj.data.sessionCombined.(sess_list{sess_no});
%             edges = 0:0.5:ceil(max(sessionData.lick_burst_idx(:,3)));
%             histogram(sessionData.lick_burst_idx(:,3),edges); title("D" + obj.data.sessionDays(sess_no)); ylabel('Lick Burst Widths (s)');
            
            sessionData.lick_burst_idx(:,6) = sessionData.lick_burst_idx(:,5) - sessionData.lick_burst_idx(:,4);
            edges = 0:ceil(max(sessionData.lick_burst_idx(:,6)));
            histogram(sessionData.lick_burst_idx(:,6),edges); title("D" + obj.data.sessionDays(sess_no)); ylabel('Lick Burst Widths (Bins)');
        end
        
    elseif(Args.LickRZEstimate)
        sess_list = fieldnames(obj.data.sessionCombined);
        
        for sess_no = 1:length(fieldnames(obj.data.sessionCombined))
            ha = subplot(ceil(length(fieldnames(obj.data.sessionCombined))/10),10,sess_no);
            sessionData = obj.data.sessionCombined.(sess_list{sess_no});
            sessionData.lick_burst_idx(:,6) = sessionData.lick_burst_idx(:,5) - sessionData.lick_burst_idx(:,4);
           
            bounds = [sessionData.lick_burst_idx(:,4) (sessionData.lick_burst_idx(:,4)+sessionData.lick_burst_idx(:,6))];
            bounds(bounds(:,1) > bounds(:,2),:) = [bounds(bounds(:,1) > bounds(:,2),2) bounds(bounds(:,1) > bounds(:,2),1)];
            combined_bins = [];
            for i = 1:size(bounds,1)
                combined_bins = [combined_bins bounds(i,1):bounds(i,2)];
            end
            combined_bins(combined_bins < 1) = combined_bins(combined_bins < 1) + 100;
            combined_bins(combined_bins > 100) = combined_bins(combined_bins > 100) - 100;
            
            y_data = histcounts(combined_bins,100);
            %y_data = histcounts(mice_data(:,5),100);
            y_data = y_data ./ sum(y_data);
            bar([-49:50], [y_data(51:end) y_data(1:50)],0.5,'r');
            title("D" + obj.data.sessionDays(sess_no));
            ylabel("Density (Lick Burst Widths)"); xlabel("Distance (Bin)");
            
            trial_no = sessionData.nTrials;
            text(5, 1.2*max(y_data), {"Trials: " + trial_no},'FontSize',9)
            ylim(ha,[0 1.4*max(y_data)])
        end
        
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
