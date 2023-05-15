function [obj, varargout] = plot(obj,varargin)
%@dirfiles/plot Plot function for dirfiles object.
%   OBJ = plot(OBJ) creates plots for the individual mice

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		  'ReturnVars',{''}, 'ArgsOnly',0, 'Cmds','', ...
          'MatchedFiringRateMap',0, 'MatchedCorrMatrix',0,'PlaceFieldPositionEntropy',0,...
          'AdSm',0,'Vel',0,'VelBinned',0,'VelDist',0,'VelPosition',0,'Lick',0,'LickPosition',0,'LickPrc',0,...
          'LickBurstWidth',0,'LickRZEstimate',0,'VelStateEstimate',0,'MatchedStats',0);
Args.flags = {'LabelsOff','ArgsOnly','MatchedFiringRateMap','MatchedCorrMatrix','PlaceFieldPositionEntropy',...
                'AdSm','Vel','VelBinned','VelDist','VelPosition','Lick','LickPosition','LickPrc','LickBurstWidth','LickRZEstimate',...
                'VelStateEstimate','MatchedStats'};
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
        vmr_combined = zeros(size(obj.data.roiMatchData.allSessionMapping,2),1);
        for sess_no = 1:size(obj.data.roiMatchData.allSessionMapping,2)
            
            subplot(3,size(obj.data.roiMatchData.allSessionMapping,2),sess_no);
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
            subplot(3,size(obj.data.roiMatchData.allSessionMapping,2),sess_no+size(obj.data.roiMatchData.allSessionMapping,2));
            
            plot(neuralData.placefieldData.(fns{cell_no}).basemapLrw, 'g'); legends{1} = sprintf('Raw'); hold on
        
            plot(neuralData.placefieldData.(fns{cell_no}).basemapLsm, 'b'); legends{2} = sprintf('AdSmooth'); hold on
            
            %plot(imgaussfilt(neuralData.placefieldData.(fns{cell_no}).basemapLsm, 3), 'c'); legends{2} = sprintf('AdSmooth'); hold on
            
            plot(neuralData.placefieldData.(fns{cell_no}).mapLsm_FFT_Low_Pass, 'r'); legends{3} = sprintf('AdSm-FFT-Low-Pass'); hold on
            
            axis([1 100 -inf inf]);
            yline(neuralData.placefieldData.(fns{cell_no}).mean_threshold,'b-'); hold on
            yline(neuralData.placefieldData.(fns{cell_no}).mean_threshold + neuralData.placefieldData.(fns{cell_no}).std_threshold,'b--'); hold on
            yline(neuralData.placefieldData.(fns{cell_no}).mean_threshold - neuralData.placefieldData.(fns{cell_no}).std_threshold,'b--'); hold on
            %         yline((neuralData.placefieldData.(fns{cell_no}).baseline_mean_height + neuralData.placefieldData.(fns{cell_no}).peak_height_from_baseline/2),'r'); hold on
            %         yline(neuralData.placefieldData.(fns{cell_no}).std_threshold,'b'); hold on
            %         yline(neuralData.placefieldData.(fns{cell_no}).std_threshold / neuralData.placefieldData.(fns{cell_no}).Args.PeakThreshold,'b--'); hold on
            
            vmr_combined(sess_no) = neuralData.placefieldData.(fns{cell_no}).vmr_threshold;
            
            text(3, 1.3*max(neuralData.placefieldData.(fns{cell_no}).basemapLrw), {"Mean: " + neuralData.placefieldData.(fns{cell_no}).mean_threshold, "StDev: " + neuralData.placefieldData.(fns{cell_no}).std_threshold, ...
                "VMR: " + neuralData.placefieldData.(fns{cell_no}).vmr_threshold},'FontSize',12); %"CV: " + neuralData.placefieldData.(fns{cell_no}).cv_threshold,
            ylim([0 1.7*max(neuralData.placefieldData.(fns{cell_no}).basemapLrw)])
            
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
                    x = linspace(-24,125,150);
                    mu = neuralData.placefieldData.(fns{cell_no}).GMM(i,1);
                    sig = neuralData.placefieldData.(fns{cell_no}).GMM(i,2);
                    amp = neuralData.placefieldData.(fns{cell_no}).GMM(i,3);
                    vo = 0;
                    y = gaus(x,mu,sig,amp,vo);
                    if mu < 50
                        y(101:125) = y(1:25);
                    else
                        y(26:50) = y(126:150);
                    end
                    % Plot gaussian
                    plot(x(26:125), y(26:125), 'k-', 'LineWidth',3);
                end
                
                x = linspace(-24,125,150);
                y = 0;
                for i = 1:size(neuralData.placefieldData.(fns{cell_no}).GMM,1)
                    y = gaus(x,neuralData.placefieldData.(fns{cell_no}).GMM(i,1),neuralData.placefieldData.(fns{cell_no}).GMM(i,2),neuralData.placefieldData.(fns{cell_no}).GMM(i,3),vo) + y;
                end
                
                legends{4} = ''; legends{5} = ''; legends{6} = ''; legends{7} = ''; legends{8} = ''; legends{9} = '';
                %plot(x(26:125), y(26:125), 'k-', 'LineWidth',3); legends{10} = sprintf('GMM'); hold on
                %legend(legends, 'Location', 'northeastoutside');
            end
            title(obj.data.sessionCombined.(sess_list{sess_no}).date); xlabel('Position (Bin)'); ylabel('Firing Rate');
            hold off
        end
        
        subplot(3,size(obj.data.roiMatchData.allSessionMapping,2),[1+2*size(obj.data.roiMatchData.allSessionMapping,2) 3*size(obj.data.roiMatchData.allSessionMapping,2)]);
        plot(vmr_combined);
        
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
                rawCorrMatrix = rawCorrMatrix.*~eye(size(rawCorrMatrix));
                a = adsmFiringMap(i,:);
                b = adsmFiringMap(j,:);
                adsmCorrMatrix(i,j) = sum(a.*b,'omitnan') / sqrt(sum(a.^2,'omitnan')*sum(b.^2,'omitnan'));
                adsmCorrMatrix = adsmCorrMatrix.*~eye(size(adsmCorrMatrix));
                a = lowpassFiringMap(i,:);
                b = lowpassFiringMap(j,:);
                lowpassCorrMatrix(i,j) = sum(a.*b,'omitnan') / sqrt(sum(a.^2,'omitnan')*sum(b.^2,'omitnan'));
                lowpassCorrMatrix = lowpassCorrMatrix.*~eye(size(lowpassCorrMatrix));
            end
        end
        
        colormap hot
        subplot(3,1,1)
        imagesc(rawCorrMatrix);; title("Raw Firing Map Correlation Matrix"); colorbar;
        subplot(3,1,2)
        imagesc(adsmCorrMatrix);; title("Adaptive Smoothed Firing Map Correlation Matrix"); colorbar;
        subplot(3,1,3)
        imagesc(lowpassCorrMatrix);; title("Low Pass Firing Map Correlation Matrix"); colorbar;
        
    elseif(Args.MatchedStats)
        
        sess_list = fieldnames(obj.data.neuralCombined);
        mean_all_combined = zeros(size(obj.data.roiMatchData.allSessionMapping,1),length(sess_list));
        std_all_combined = zeros(size(obj.data.roiMatchData.allSessionMapping,1),length(sess_list));
        vmr_all_combined = zeros(size(obj.data.roiMatchData.allSessionMapping,1),length(sess_list));
        for cell_list_idx = 1:size(obj.data.roiMatchData.allSessionMapping,1)
            for sess_no = 1:size(obj.data.roiMatchData.allSessionMapping,2)
                cell_no = obj.data.roiMatchData.allSessionMapping(cell_list_idx,sess_no);
                neuralData = obj.data.neuralCombined.(sess_list{sess_no});
                fns = fieldnames(neuralData.placefieldData);
                if ~isempty(neuralData.placefieldData.(fns{cell_no}))
                    mean_all_combined(cell_list_idx,sess_no) = neuralData.placefieldData.(fns{cell_no}).mean_threshold;
                    std_all_combined(cell_list_idx,sess_no) = neuralData.placefieldData.(fns{cell_no}).std_threshold;
                    vmr_all_combined(cell_list_idx,sess_no) = neuralData.placefieldData.(fns{cell_no}).vmr_threshold;
                else
                    mean_all_combined(cell_list_idx,sess_no) = NaN;
                    std_all_combined(cell_list_idx,sess_no) = NaN;
                    vmr_all_combined(cell_list_idx,sess_no) = NaN;
                end
            end
        end
        
        subplot(3,1,1)
        plot(mean_all_combined','color',[0,0,0,0.05]);
        hold on
        plot(mean(mean_all_combined,'omitnan'),'r','LineWidth',4);
        hold off
        title('Mean'); xlabel('Sessions');
        subplot(3,1,2)
        plot(std_all_combined','color',[0,0,0,0.05]);
        hold on
        plot(mean(std_all_combined,'omitnan'),'r','LineWidth',4);
        hold off
        title('StDev'); xlabel('Sessions');
        subplot(3,1,3)
        plot(vmr_all_combined','color',[0,0,0,0.05]);
        hold on
        plot(mean(vmr_all_combined,'omitnan'),'r','LineWidth',4);
        hold off
        title('VMR'); xlabel('Sessions');
               
    elseif(Args.PlaceFieldPositionEntropy)
        sess_list = fieldnames(obj.data.neuralCombined);
        
        if n == 1
            for sess_no = 1:size(sess_list,1)
                subplot(2,ceil(size(sess_list,1)/2),sess_no);
                title("Place field positions"); xlabel('Position Bins'); ylabel('Count');
                neuralData = obj.data.neuralCombined.(sess_list{sess_no});
                sessionData = obj.data.sessionCombined.(sess_list{sess_no});
                                
                y_data = histcounts(neuralData.placefieldStats(:,2),100);
                y_data = y_data ./ sum(y_data);
                bar([-49:50], [y_data(51:end) y_data(1:50)]);
                title("D" + obj.data.sessionDays(sess_no) + "; " + sessionData.date);
            end
        elseif n == 2
            for sess_no = 1:size(sess_list,1)
                subplot(2,ceil(size(sess_list,1)/2),sess_no);
                title("Place field positions"); xlabel('Position Bins'); ylabel('Count');
                neuralData = obj.data.neuralCombined.(sess_list{sess_no});
                sessionData = obj.data.sessionCombined.(sess_list{sess_no});
                
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
                title("D" + obj.data.sessionDays(sess_no) + "; " + sessionData.date);
            end
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
        %prc_vel = zeros(length(fieldnames(obj.data.sessionCombined)),1);
        mean_trial_velocity = zeros(length(fieldnames(obj.data.sessionCombined)),1);
        for sess_no = 1:length(fieldnames(obj.data.sessionCombined))
            sessionData = obj.data.sessionCombined.(sess_list{sess_no});
%             mid_nTrial = floor(sessionData.nTrials/2);
%             velocity_subset = sessionData.velocity_averaged(find(sessionData.velocity_averaged(:,2) > mid_nTrial),:);
%             mean_velocity(sess_no) = mean(velocity_subset(velocity_subset(:,4) ~= 0,4), 'omitnan');
            
            mean_velocity(sess_no) = mean(sessionData.velocity_averaged(sessionData.velocity_averaged(:,4) ~= 0,4), 'omitnan');
            mean_track_velocity(sess_no) = mean(sessionData.velocity_binned, 'all', 'omitnan');
            %prc_vel(sess_no) = nnz(sessionData.velocity_averaged(sessionData.velocity_averaged(:,4) ~= 0,4) < mean_velocity(sess_no)) / nnz(sessionData.velocity_averaged(sessionData.velocity_averaged(:,4) ~= 0,4));
            mean_trial_velocity(sess_no) = mean(220./sessionData.data_trial(:,4), 'omitnan');

        end
        
        subplot(5,1,1)
        plot(mean_velocity); title('Mean Velocity'); xlabel('Sessions'); ylabel('Velocity (cm/s)');
        yline(median(mean_velocity,'omitnan'),'r')
        yline(mean(mean_velocity,'omitnan'),'g')
        
        subplot(5,1,2)
        plot(mean_track_velocity); title('Mean Velocity (Binned)'); xlabel('Sessions'); ylabel('Velocity (cm/s)');
        yline(median(mean_track_velocity,'omitnan'),'r')
        yline(mean(mean_track_velocity,'omitnan'),'g')
        
        subplot(5,1,3)
        plot(obj.data.sessionDays, mean_velocity); title('Mean Velocity'); xlabel('Sessions (Day)'); ylabel('Velocity (cm/s)');
        yline(median(mean_velocity,'omitnan'),'r')
        yline(mean(mean_velocity,'omitnan'),'g')
        
        subplot(5,1,4)
        plot(obj.data.sessionDays, mean_track_velocity); title('Mean Velocity (Binned)'); xlabel('Sessions (Day)'); ylabel('Velocity (cm/s)');
        yline(median(mean_track_velocity,'omitnan'),'r')
        yline(mean(mean_track_velocity,'omitnan'),'g')
        
        subplot(5,1,5)
        plot(mean_trial_velocity); title('Mean Trial Velocity'); xlabel('Sessions'); ylabel('Velocity (cm/s)');
        yline(median(mean_trial_velocity,'omitnan'),'r')
        yline(mean(mean_trial_velocity,'omitnan'),'g')
        yline(1,'k--')
        
    elseif(Args.VelBinned)
        
%         sess_list = fieldnames(obj.data.sessionCombined);
%         sessionData = obj.data.sessionCombined.(sess_list{n});
%                
%         velocity_binned = zeros(sessionData.nTrials, sessionData.Args.BinSize);
%         for i = 1:sessionData.nTrials
%             for j = 1:sessionData.Args.BinSize
%                 velocity_binned = sessionData.velocity_binned;
%             end
%         end
%         
%         % Smoothing
%         for row = 1:sessionData.nTrials
%             velocity_binned(row,:) = smooth(velocity_binned(row,:));
%         end
%                 
%         imagesc(velocity_binned); colormap hot
%         title(sessionData.date);
        
%         for block_no = 1:12
%             subplot(3,4,block_no)
%             sess_list = fieldnames(obj.data.sessionCombined);
%             session_no = 12 * (n - 1) + block_no;
%             
%             if session_no <= length(fieldnames(obj.data.sessionCombined))
%                 sessionData = obj.data.sessionCombined.(sess_list{session_no});
%                 
% %                 velocity_binned = zeros(sessionData.nTrials, sessionData.Args.BinSize);
% %                 for i = 1:sessionData.nTrials
% %                     for j = 1:sessionData.Args.BinSize
% %                         velocity_binned = sessionData.velocity_binned;
% %                     end
% %                 end
%                 velocity_binned = sessionData.velocity_binned;
%                 velocity_binned = [velocity_binned(51:end); velocity_binned(1:50)];
%                 
%                 % Smoothing
%                 for row = 1:sessionData.nTrials
%                     velocity_binned(row,:) = smooth(velocity_binned(row,:));
%                 end
%                 
%                 imagesc(velocity_binned'); colormap hot; colorbar
%                 title(sessionData.date);
%             end
%         end
        
        sess_list = fieldnames(obj.data.sessionCombined);
        index = reshape(1:ceil(length(fieldnames(obj.data.sessionCombined))/5)*5,5,[]).';
        for sess_no = 1:length(fieldnames(obj.data.sessionCombined))
            subplot(ceil(length(fieldnames(obj.data.sessionCombined))/5),5,index(sess_no))
            sessionData = obj.data.sessionCombined.(sess_list{sess_no});
            
            velocity_binned = sessionData.velocity_binned;
            velocity_binned = [velocity_binned(51:end); velocity_binned(1:50)];
                        
            imagesc(velocity_binned'); colormap hot; colorbar
            title("D" + obj.data.sessionDays(sess_no) + "; " + sessionData.date);
        end
            
    elseif(Args.VelDist)
        sess_list = fieldnames(obj.data.sessionCombined);
        
        index = reshape(1:ceil(length(fieldnames(obj.data.sessionCombined))/5)*5,5,[]).';
        for sess_no = 1:length(fieldnames(obj.data.sessionCombined))
            subplot(ceil(length(fieldnames(obj.data.sessionCombined))/5),5,index(sess_no))
            sessionData = obj.data.sessionCombined.(sess_list{sess_no});
            
%             edges = 0:0.5:ceil(max(sessionData.velocity_averaged(:,4)));
%             histogram(sessionData.velocity_averaged(sessionData.velocity_averaged(:,4) ~= 0,4),edges); title("D" + obj.data.sessionDays(sess_no));% ylabel('Velocity (cm/s)');
            
%             edges = 0:0.5:ceil(max(sessionData.velocity_binned));
%             histogram(sessionData.velocity_binned,edges); title("D" + obj.data.sessionDays(sess_no));% ylabel('Velocity (cm/s)');
            
%             [y_data,edges] = histcounts(sessionData.velocity_averaged(sessionData.velocity_averaged(:,4) ~= 0,4),edges,'Normalization','pdf'); title("D" + obj.data.sessionDays(sess_no));
%             edges = edges(2:end) - (edges(2)-edges(1))/2;
%             plot(edges,y_data);
%             xlim([0 30]);
%             xline(mean(sessionData.velocity_averaged(sessionData.velocity_averaged(:,4) ~= 0,4), 'omitnan'));

            edges = 0:0.5:8;
            histogram(220./sessionData.data_trial(:,4),edges); title("D" + obj.data.sessionDays(sess_no));% ylabel('Velocity (cm/s)');
            xline(mean(220./sessionData.data_trial(:,4), 'omitnan'));
        end
      
    elseif(Args.VelPosition)
        
        sess_list = fieldnames(obj.data.sessionCombined);
        index = reshape(1:ceil(length(fieldnames(obj.data.sessionCombined))/5)*5,5,[]).';
        for sess_no = 1:length(fieldnames(obj.data.sessionCombined))
            subplot(ceil(length(fieldnames(obj.data.sessionCombined))/5),5,index(sess_no))
            sessionData = obj.data.sessionCombined.(sess_list{sess_no});
            
            vel_binned = reshape(sessionData.data_bin(:,8),100,sessionData.nTrials)';
            vel_binned = [zeros(1,100); vel_binned;]; % Pad one trial (100 bins) before
            vel_binned = reshape(vel_binned.',1,[]); % Flatten array
            vel_binned = circshift(vel_binned,-50); % Shift array 50 bins behind
            vel_binned = reshape(vel_binned,100,sessionData.nTrials+1)'; % Restore array
            imagesc(vel_binned); colorbar;
            title("D" + obj.data.sessionDays(sess_no) + "; " + sessionData.date);
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
        
    elseif(Args.LickPosition)
%         sess_no = n;
%         sess_list = fieldnames(obj.data.sessionCombined);
%         sessionData = obj.data.sessionCombined.(sess_list{sess_no});
        
        sess_list = fieldnames(obj.data.sessionCombined);
        index = reshape(1:ceil(length(fieldnames(obj.data.sessionCombined))/5)*5,5,[]).';
        for sess_no = 1:length(fieldnames(obj.data.sessionCombined))
            subplot(ceil(length(fieldnames(obj.data.sessionCombined))/5),5,index(sess_no))
            sessionData = obj.data.sessionCombined.(sess_list{sess_no});
        
            lick_count_binned = sessionData.lick_count_binned';
            lick_count_binned = [zeros(1,100); lick_count_binned;]; % Pad one trial (100 bins) before
            lick_count_binned = reshape(lick_count_binned.',1,[]); % Flatten array
            lick_count_binned = circshift(lick_count_binned,-50); % Shift array 50 bins behind
            lick_count_binned = reshape(lick_count_binned,100,sessionData.nTrials+1)'; % Restore array
            imagesc(lick_count_binned); colorbar;
            title("D" + obj.data.sessionDays(sess_no) + "; " + sessionData.date);
        end
        
    elseif(Args.LickPrc)
        sess_list = fieldnames(obj.data.sessionCombined);
        
        lick_prc_RZ = zeros(length(fieldnames(obj.data.sessionCombined)),3);
        for sess_no = 1:length(fieldnames(obj.data.sessionCombined))
            sessionData = obj.data.sessionCombined.(sess_list{sess_no});
            lick_prc = zeros(100,1);
            for i = 1:100
                lick_prc(i) = nnz(sessionData.lick_timestamps_spliced(find(sessionData.lick_timestamps_spliced(:,3) == i), [1])) / nnz(sessionData.lick_timestamps_adjusted(:,1));
            end
            lick_prc_RZ(sess_no,1) = sum(lick_prc(96:100));
            lick_prc_RZ(sess_no,2) = sum(lick_prc(1:5));
            lick_prc_RZ(sess_no,3) = lick_prc_RZ(sess_no,1) + lick_prc_RZ(sess_no,2);
        end
        
        subplot(2,1,1)
        plot(lick_prc_RZ(:,3),'x'); title('Lick % RZ'); xlabel('Sessions'); ylabel('%');
        
        subplot(2,1,2)
        plot(obj.data.sessionDays, lick_prc_RZ(:,3),'x'); title('Lick % RZ'); xlabel('Sessions (Day)'); ylabel('%');

    elseif(Args.LickBurstWidth)
        sess_list = fieldnames(obj.data.sessionCombined);
        
        lick_burst_width_RZ = [];
        entropy_combined = [];
        centre_score_combined = [];
        for sess_no = 1:length(fieldnames(obj.data.sessionCombined))
              if n == 1
                ha = subplot(ceil(length(fieldnames(obj.data.sessionCombined))/5),5,sess_no);
              end
%             
%             sessionData = obj.data.sessionCombined.(sess_list{sess_no});
% %             edges = 0:0.5:ceil(max(sessionData.lick_burst_idx(:,3)));
% %             histogram(sessionData.lick_burst_idx(:,3),edges); title("D" + obj.data.sessionDays(sess_no)); ylabel('Lick Burst Widths (s)');
%             
%             sessionData.lick_burst_idx(:,6) = sessionData.lick_burst_idx(:,5) - sessionData.lick_burst_idx(:,4);
%             edges = 0:ceil(max(sessionData.lick_burst_idx(:,6)));
%             histogram(sessionData.lick_burst_idx(:,6),edges); title("D" + obj.data.sessionDays(sess_no)); ylabel('Lick Burst Widths (Bins)');
            
             sessionData = obj.data.sessionCombined.(sess_list{sess_no});
             sess_data = [sessionData.lick_burst_idx(:,3) mod(sessionData.lick_burst_idx(:,4),100) mod(sessionData.lick_burst_idx(:,5),100)];
             sess_data(:,4) = sess_data(:,3) - sess_data(:,2);
             sess_data(sess_data(:,2) == 0,2) = 100;
             sess_data(sess_data(:,3) == 0,3) = 100;
             
             bounds = [sess_data(:,2) (sess_data(:,2)+sess_data(:,4))];
             bounds(bounds(:,1) > bounds(:,2),:) = [bounds(bounds(:,1) > bounds(:,2),2) bounds(bounds(:,1) > bounds(:,2),1)];
             combined_bins = [];
             for i = 1:size(bounds,1)
                 combined_bins = [combined_bins bounds(i,1):bounds(i,2)];
             end
             combined_bins(combined_bins < 1) = combined_bins(combined_bins < 1) + 100;
             combined_bins(combined_bins > 100) = combined_bins(combined_bins > 100) - 100;
             
             y_data = histcounts(combined_bins,100);
             y_data = y_data ./ sum(y_data);
             adjusted_y_data = [y_data(51:end) y_data(1:50)];
             
             trial_no = sessionData.nTrials;
             entropy_score = round(exp(-sum(y_data .* log(y_data),2,'omitnan'))/100,3,'significant');
             centre_score = round(sum(adjusted_y_data .* [linspace(-1,1,50) linspace(1,-1,50)],'omitnan'),3,'significant');
             
             entropy_combined = [entropy_combined entropy_score];
             centre_score_combined = [centre_score_combined centre_score];
             
             if n == 1
                 plot([-49:50], [y_data(51:end) y_data(1:50)],'r');
                 title(sessionData.date);
                 ylabel("Density"); xlabel("Distance (Bin)");
                 
                 text(5, 1.3*max(y_data), {"Trials: " + trial_no, "Entropy: " + entropy_score,"Centre Score: " + centre_score},'FontSize',9);
                 ylim(ha,[0 1.7*max(y_data)])
             end
                                    
             lick_burst_width_RZ = [lick_burst_width_RZ; [sum(y_data(96:100)) sum(y_data(1:5))]];
             %y_data_combined = [y_data_combined; y_data];
        end
        lick_burst_width_RZ(:,3) = lick_burst_width_RZ(:,1) + lick_burst_width_RZ(:,2);
        
        if n == 2
            subplot(2,1,1)
            plot(lick_burst_width_RZ(:,3),'x'); title('Lick Burst RZ'); xlabel('Sessions'); ylabel('%');
            
            subplot(2,1,2)
            plot(obj.data.sessionDays, lick_burst_width_RZ(:,3),'x'); title('Lick Burst RZ'); xlabel('Sessions (Day)'); ylabel('%');
        end
        
        if n == 3
            subplot(4,1,1)
            plot(entropy_combined,'x'); title('Entropy'); xlabel('Sessions');
            
            subplot(4,1,2)
            plot(centre_score_combined,'x'); title('Centre Score'); xlabel('Sessions');
            
            subplot(4,1,3)
            plot(obj.data.sessionDays, entropy_combined,'x'); title('Entropy'); xlabel('Sessions (Days');
            
            subplot(4,1,4)
            plot(obj.data.sessionDays, centre_score_combined,'x'); title('Centre Score'); xlabel('Sessions (Days)');
        end
        
        if n == 4
            y_data_combined = sum(adjusted_y_data,1,'omitnan');
            bar([-49:50], y_data_combined,0.5,'r');
            ylabel("Density (Lick Burst Widths)"); xlabel("Distance (Bin)");
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
        
        
    elseif(Args.VelStateEstimate)
        
        hold off
        sess_list = fieldnames(obj.data.sessionCombined);
                
        mean_velocity = zeros(length(fieldnames(obj.data.sessionCombined)),1);
        mean_track_velocity = zeros(length(fieldnames(obj.data.sessionCombined)),1);
        %prc_vel = zeros(length(fieldnames(obj.data.sessionCombined)),1);
        mean_trial_velocity = zeros(length(fieldnames(obj.data.sessionCombined)),1);
        lick_prc_RZ = zeros(length(fieldnames(obj.data.sessionCombined)),3);
        for sess_no = 1:length(sess_list)
            sessionData = obj.data.sessionCombined.(sess_list{sess_no});
%             mid_nTrial = floor(sessionData.nTrials/2);
%             velocity_subset = sessionData.velocity_averaged(find(sessionData.velocity_averaged(:,2) > mid_nTrial),:);
%             mean_velocity(sess_no) = mean(velocity_subset(velocity_subset(:,4) ~= 0,4), 'omitnan');
            
            mean_velocity(sess_no) = mean(sessionData.velocity_averaged(sessionData.velocity_averaged(:,4) ~= 0,4), 'omitnan');
            mean_track_velocity(sess_no) = mean(sessionData.velocity_binned, 'all', 'omitnan');
            %prc_vel(sess_no) = nnz(sessionData.velocity_averaged(sessionData.velocity_averaged(:,4) ~= 0,4) < mean_velocity(sess_no)) / nnz(sessionData.velocity_averaged(sessionData.velocity_averaged(:,4) ~= 0,4));
            mean_trial_velocity(sess_no) = mean(220./sessionData.data_trial(:,4), 'omitnan');
            
            lick_prc = zeros(length(sess_list),100);
            for bin_no = 1:100
                lick_prc(sess_no,bin_no) = (nnz(sessionData.lick_timestamps_spliced(find(sessionData.lick_timestamps_spliced(:,3) == bin_no), [1])) / nnz(sessionData.lick_timestamps_adjusted(:,1))) * 100;
                if isnan(lick_prc(sess_no,bin_no))
                    lick_prc(sess_no,bin_no) = 0;
                end
            end
            lick_prc_RZ(sess_no,1) = sum(lick_prc(sess_no,96:100));
            lick_prc_RZ(sess_no,2) = sum(lick_prc(sess_no,1:5));
            lick_prc_RZ(sess_no,3) = lick_prc_RZ(sess_no,1) + lick_prc_RZ(sess_no,2);
        end
        
        mean_vel = mean_velocity;
        padSize = floor(5/2); % 5 days sliding window
        mean_vel_padded = padarray(mean_vel, padSize, 'replicate');
        mean_vel_window_smooth = zeros(size(mean_vel,1), 1);
        for i = 1+padSize:size(mean_vel, 1)+padSize
            mean_vel_window_smooth(i - padSize) = mean(mean_vel_padded((i - padSize):(i + padSize)));
        end
        
        upper_vel_thres = mean(mean_vel(ceil(length(mean_vel)/2):end),'omitnan');
        %lower_vel_thres = mean_vel(1);
        lower_vel_thres = min(mean_vel(1:ceil(length(mean_vel)/2)),[],'omitnan');
        
        % Assumption 1:
        % Initial state (low velocity regime) can be estimated using mean velocity of first session
        % Assumption 2:
        % Final state (high velocity regime) can be estimated using mean velocity of second half of all sessions
        
        mu = zeros(length(sess_list),2);
        mu(1,:) = [lower_vel_thres upper_vel_thres];
        sigma = zeros(length(sess_list),1,1,2);
        sigma(1,:) = cat(3,[1],[1]);
        mix_prop = zeros(length(sess_list),2);
        mix_prop(1,:) = [0.99 0.01];%[0.99 0.01];
        
        decision_state = zeros(length(sess_list),5);
        decision_state2 = zeros(length(sess_list),7);
        
        %index = reshape(1:ceil(length(fieldnames(miceData.data_trial))/5)*5,5,[]).';
        for sess_no = 1:length(sess_list)
            %subplot(ceil(length(fieldnames(miceData.data_trial))/5),5,index(sess_no))
                       
            if sess_no > 1
                mu(sess_no,:) = [lower_vel_thres upper_vel_thres];
                sigma(sess_no,:) = cat(3,[1],[1]);
                mix_prop(sess_no,:) = decision_state(sess_no-1,1:2);
                if mix_prop(sess_no,1) == 0
                    mix_prop(sess_no,1) = mix_prop(sess_no,1) + eps;
                end
            end
    
            % Obtain prior
            gm = gmdistribution(mu(sess_no,:)',sigma(sess_no,:,:,:),mix_prop(sess_no,:));
            gmPDF = pdf(gm,[0:0.1:40]');
            gmPDF = gmPDF / sum(gmPDF);
            %plot(0:0.1:40,gmPDF,'b-');
            %title(fns{mice_id(n,2)} + " - D" + miceData.sessionDays(sess_no));
            %%xlabel('Velocity (cm/s)'); xlim([0 40]);
            %hold on
            
            % Calculate observation (likelihood)
%             bootstrap_vel = bootstrp(100,@median,sessionData.velocity_binned);
%             [fi,xi] = ksdensity(bootstrap_vel);
%             fi = fi ./ sum(fi);% / 20
%             %plot(xi,fi,'kx'); xlim([0 40]);
%             %hold on
%             %f = fit(xi',fi','gauss1');
%             %plot(f,'k-');
%             GMModel = fitgmdist([xi]',1);
%             gmPDF2 = pdf(GMModel,[0:0.1:40]');
%             gmPDF2 = gmPDF2 ./ sum(gmPDF2);
%             %plot(0:0.1:40,gmPDF2,'r-');
            
            %observation = mean(sessionData.velocity_binned(51:end),'omitnan');
            %observation = mean(bootstrap_vel,'omitnan');
            %observation = mean_vel_window_smooth(sess_no);
            observation = mean_vel(sess_no);
            [gm_posterior,nlogL] = posterior(gm,observation);
            %[gm_posterior,nlogL] = posterior(gm,gmPDF2);
            
            covNames = {'diagonal','full'};
            CovType = find(strncmpi(gm.CovType,covNames,length(gm.CovType)));
            nlog_lh = wdensity_dupe(observation,gm.mu, gm.Sigma, gm.PComponents, gm.SharedCov, CovType) * -1; %nLogL for individual components
            
            decision_state(sess_no,:) = [gm_posterior nlogL nlog_lh];
                        
        end

        upper_vel_thres2 = mean(mean_vel_window_smooth(ceil(length(mean_vel_window_smooth)/2):end),'omitnan');
        lower_vel_thres2 = min(mean_vel_window_smooth(1:ceil(length(mean_vel_window_smooth)/2)),[],'omitnan');
        
        %index = reshape(1:ceil(length(fieldnames(miceData.data_trial))/5)*5,5,[]).';
        for sess_no = 1:length(sess_list)
            %subplot(ceil(length(fieldnames(miceData.data_trial))/5),5,index(sess_no))
            
            if sess_no > 1
                mu(sess_no,:) = [lower_vel_thres2 upper_vel_thres2];
                sigma(sess_no,:) = cat(3,[1],[1]);
                mix_prop(sess_no,:) = decision_state2(sess_no-1,1:2);
                %mix_prop(sess_no,1) = mix_prop(sess_no,1) / (2^decision_state(sess_no-1,3)); % Scaling factor - Pressure to switch to high velocity state
                if mix_prop(sess_no,1) == 0
                    mix_prop(sess_no,1) = mix_prop(sess_no,1) + eps;
                end
            end
            
            % Obtain prior
            gm = gmdistribution(mu(sess_no,:)',sigma(sess_no,:,:,:),mix_prop(sess_no,:));
            gmPDF = pdf(gm,[0:0.1:40]');
            gmPDF = gmPDF / sum(gmPDF);
            %plot(0:0.1:40,gmPDF,'b-');
            %title(fns{mice_id(n,2)} + " - D" + miceData.sessionDays(sess_no));
            %%xlabel('Velocity (cm/s)'); xlim([0 40]);
            %hold on
            
            % Calculate observation (likelihood)
%             bootstrap_vel = bootstrp(100,@median,sessionData.velocity_binned);
%             [fi,xi] = ksdensity(bootstrap_vel);
%             fi = fi ./ sum(fi);% / 20
%             %plot(xi,fi,'kx'); xlim([0 40]);
%             %hold on
%             %f = fit(xi',fi','gauss1');
%             %plot(f,'k-');
%             GMModel = fitgmdist([xi]',1);
%             gmPDF2 = pdf(GMModel,[0:0.1:40]');
%             gmPDF2 = gmPDF2 ./ sum(gmPDF2);
%             %plot(0:0.1:40,gmPDF2,'r-');
            
            %observation = mean(sessionData.velocity_binned,'omitnan');
            %observation = mean(bootstrap_vel,'omitnan');
            observation = mean_vel_window_smooth(sess_no);
            [gm_posterior,nlogL] = posterior(gm,observation);
            %[gm_posterior,nlogL] = posterior(gm,gmPDF2);
            
            covNames = { 'diagonal','full'};
            CovType = find(strncmpi(gm.CovType,covNames,length(gm.CovType)));
            [nlog_lh nlog_lh_update] = wdensity_dupe(observation,gm.mu, gm.Sigma, gm.PComponents, gm.SharedCov, CovType); %nLogL for individual components
            nlog_lh = nlog_lh * -1; nlog_lh_update = nlog_lh_update * -1; 
            decision_state2(sess_no,:) = [gm_posterior nlogL nlog_lh nlog_lh_update];

        end
        
        delete(findall(gcf,'Type','textbox'))
        subplot(5,1,1)
        yyaxis left
        plot(mean_vel,'b.'); xlabel('Sessions'); ylabel('Velocity (cm/s)');
        hold on
        plot(mean_vel_window_smooth,'rx'); xlabel('Sessions'); ylabel('Velocity (cm/s)');
        hold on
        %plot(mean_vel_smooth,'b-');
        yline(lower_vel_thres,'g--')
        yline(upper_vel_thres,'r--')
        hold on
        yyaxis right
        plot(decision_state(:,2),'k-'); ylim([0 1.2]); ylabel('State');
        hold on
        plot(decision_state2(:,2),'c-');
        hold on
        
        subplot(5,1,2)
        plot(decision_state(:,3),'k-'); xlabel('Sessions'); ylabel('nLogL'); % Overall nLogL
        hold on
        plot(decision_state2(:,3),'c-'); % Overall nLogL
        
        subplot(5,1,3)
        hold all
        plot(decision_state(:,4),'k--','LineWidth',2); xlabel('Sessions'); ylabel('nLogL');  % nLogL of component 1
        plot(decision_state(:,5),'k:','LineWidth',2); % nLogL of component 2
        
        plot(decision_state2(:,4),'c--','LineWidth',2); %nLogL of component 1
        plot(decision_state2(:,5),'c:','LineWidth',2); %nLogL of component 2
        
        TextLocation(sprintf("""- -"" - Component 1 (Low vel)\n"". ."" - Component 2 (High vel)"),'Location','northeast');
        
        subplot(5,1,4)
        plot(decision_state(:,4) - decision_state(:,5),'k--','LineWidth',2); xlabel('Sessions'); ylabel('nLogL');
        hold on
        plot(decision_state2(:,4) - decision_state2(:,5),'c--','LineWidth',2); xlabel('Sessions'); ylabel('nLogL');
        hold on
        yline(0,'r:','LineWidth',2)
        ylim([-10 50]); xlim([0 size(mean_vel,1)]); 
        
        subplot(5,1,5)
        plot(lick_prc_RZ(:,3),'rx'); xlabel('Sessions'); ylabel('RZ Lick (%)');
        hold on
        
        % fit to a*(1-exp(-b.*x)) model;
        expc1 = fittype('a*(1-exp(-b.*x))','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b'});% a*(1-exp(-b.*x))
        [cf,gof] = fit([1:size(lick_prc_RZ(:,3),1)]',lick_prc_RZ(:,3),expc1,'Lower',[0,0],'Upper',[Inf,Inf]);
        plot(cf,'r-'); ylim([0,110]); %p1.LineWidth=1; %xlim([0,22]);
        hold off
        delete(findobj('type','legend'))
        yline(0.75*cf(length(sess_list)),'r--');
        
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
