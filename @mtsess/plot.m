function [obj, varargout] = plot(obj,varargin)
%@dirfiles/plot Plot function for dirfiles object.
%   OBJ = plot(OBJ) creates plots for the session

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		  'ReturnVars',{''}, 'ArgsOnly',0, 'Cmds','', ...
          'VelRaw',0, 'TrialVelRaw',0, 'Vel',0, 'TrialVel',0, 'TrialVelFilt',0, 'VelBinned',0, ...
          'VelCount',0, 'WaterLick',0, 'TrialWaterLick',0, 'LickDistribution',0, 'TrialLickDistribution',0, ...
          'LickBinned',0, 'TrialLickBinned',0, 'LickRate',0, 'TrialLickRate',0, 'TrialLick',0, 'LickRZ',0, ...
          'LickBurst',0,'TrialLickBurst',0,'AnticipatoryLick',0,'SummaryPlot',0,'TrialdFF0',0,'PopVec',0,'FRVel',0);
Args.flags = {'LabelsOff','ArgsOnly','VelRaw','TrialVelRaw','Vel', 'TrialVel', 'TrialVelFilt', ...
              'VelBinned', 'VelCount', 'WaterLick','TrialWaterLick', ...
              'TrialLickDistribution', 'LickDistribution', 'LickBinned', 'TrialLickBinned', ...
              'LickRate', 'TrialLickRate', 'TrialLick', 'LickRZ', 'LickBurst', 'TrialLickBurst', ...
              'AnticipatoryLick','SummaryPlot','TrialdFF0','PopVec','FRVel'};
[Args,varargin2] = getOptArgs(varargin,Args);

% if user select 'ArgsOnly', return only Args structure for an empty object
if Args.ArgsOnly
    Args = rmfield (Args, 'ArgsOnly');
    varargout{1} = {'Args',Args};
    return;
end

if(~isempty(Args.NumericArguments))
	n = Args.NumericArguments{1};
	if(Args.VelRaw)
		% Trial Distance-Raw Velocity
        subplot(2,1,1);
        plot(obj.data.session_data_exclude_zero_trials(:,1), obj.data.session_data_exclude_zero_trials(:,7)*220); title('Distance'); xlabel('Time (s)'); ylabel('Distance (cm)'); xline(obj.data.data_trial(:,3));
        subplot(2,1,2);
        plot(obj.data.session_data_exclude_zero_trials(:,1), obj.data.session_data_exclude_zero_trials(:,8)); title('Velocity'); xlabel('Time (s)'); ylabel('Velocity (cm/s)'); xline(obj.data.data_trial(:,3));

    elseif(Args.TrialVelRaw)
		% Distance-Raw Velocity
        trial = n;
        plot_start_idx = obj.data.TrialTime_idx(trial,1);
        plot_end_idx = obj.data.TrialTime_idx(trial,2);
        subplot(2,1,1);
        plot(obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,7)*220); title('Distance'); xlabel('Time (s)'); ylabel('Distance (cm)'); xline(obj.data.data_trial(trial,3));
        subplot(2,1,2);
        plot(obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,8)); title('Velocity'); xlabel('Time (s)'); ylabel('Velocity (cm/s)'); xline(obj.data.data_trial(trial,3));

    elseif(Args.Vel)
        % Distance-Velocity
        subplot(3,1,1);
        plot(obj.data.session_data_exclude_zero_trials(:,1), obj.data.session_data_exclude_zero_trials(:,7)*220); title('Distance'); xlabel('Time (s)'); ylabel('Distance (cm)'); xline(obj.data.data_trial(:,3));
        subplot(3,1,2);
        plot(obj.data.session_data_exclude_zero_trials(:,1), obj.data.velocity_averaged(:,4)); title('Velocity'); xlabel('Time (s)'); ylabel('Velocity (cm/s)'); xline(obj.data.data_trial(:,3));
        subplot(3,1,3);
        plot(obj.data.session_data_exclude_zero_trials(:,1), obj.data.velocity_averaged_filt(:,4)); title('Velocity'); xlabel('Time (s)'); ylabel('Velocity (Thresholded) (cm/s)'); xline(obj.data.data_trial(:,3));

    elseif(Args.TrialVel)
        % Trial Distance-Velocity
        trial = n;
        plot_start_idx = obj.data.TrialTime_idx(trial,1);
        plot_end_idx = obj.data.TrialTime_idx(trial,2);
        subplot(3,1,1);
        plot(obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,7)*220); title('Distance'); xlabel('Time (s)'); ylabel('Distance (cm)'); xline(obj.data.data_trial(trial,3));
        subplot(3,1,2);
        plot(obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), obj.data.velocity_averaged(plot_start_idx:plot_end_idx,4)); title('Velocity'); xlabel('Time (s)'); ylabel('Velocity (cm/s)'); xline(obj.data.data_trial(trial,3));
        subplot(3,1,3);
        plot(obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), obj.data.velocity_averaged_filt(plot_start_idx:plot_end_idx,4)); title('Velocity (Thresholded)'); xlabel('Time (s)'); ylabel('Velocity (cm/s)'); xline(obj.data.data_trial(trial,3));

        
    elseif(Args.TrialVelFilt)
        % Distance-Filtered Velocity
        trial = n;
        plot_start_idx = obj.data.TrialTime_idx(trial,1);
        plot_end_idx = obj.data.TrialTime_idx(trial,2);
        subplot(2,1,1);
        plot(obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,7)*220); title('Distance'); xlabel('Time (s)'); ylabel('Distance (cm)'); xline(obj.data.data_trial(trial,3));
        subplot(2,1,2);
        plot(obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), obj.data.velocity_averaged_filt(plot_start_idx:plot_end_idx,4)); title('Velocity'); xlabel('Time (s)'); ylabel('Velocity (cm/s)'); xline(obj.data.data_trial(trial,3));

 	elseif(Args.VelBinned)
        % Velocity Binned
        %plot(1:size(obj.data.velocity_binned(:,1),1), obj.data.velocity_binned(:,1), 'b'); title('Mean Velocity along Track Distance'); xlabel('Distance (Bin)'); ylabel('Mean Velocity (cm/s)');
%         for i = 1:size(obj.data.velocity_binned,1)
%             obj.data.velocity_binned(i,:) = smooth(obj.data.velocity_binned(i,:),3);
%         end
        imagesc(obj.data.velocity_binned); colorbar;
        title('Mean Velocity (cm/s) along Track Distance','FontSize',20); xlabel('Distance (Bin)','FontSize',18); ylabel('Lap','FontSize',18);
    
    elseif(Args.VelCount)
        % Velocity Binned
        histogram(obj.data.velocity_averaged(obj.data.velocity_averaged(:,4) ~= 0,4)); title('Velocity Distribution'); xlabel('Velocity (cm/s)'); ylabel('Count');

    elseif(Args.WaterLick)   
        % Water-Lick
        plot(obj.data.session_data_exclude_zero_trials(:,1), obj.data.session_data_exclude_zero_trials(:,5)*2, 'b'); xlabel('Time (s)');
        hold on
        plot(obj.data.session_data_exclude_zero_trials(:,1), obj.data.session_data_exclude_zero_trials(:,10), 'r'); title('Licking'); xlabel('Time (s)');% xline(obj.data.data_trial(:,3));
        hold off
        
     elseif(Args.TrialWaterLick)   
        % Trial Water-Lick
        trial = n;
        plot_start_idx = obj.data.TrialTime_idx(trial,1);
        plot_end_idx = obj.data.TrialTime_idx(trial,2);
        plot(obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,5)*2, 'b'); xlabel('Time (s)');
        hold on
        plot(obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,6), 'r'); title('Licking'); xlabel('Time (s)');% xline(obj.data.data_trial(:,3));
        hold off

    elseif(Args.LickDistribution)
        % Lick Distribution
        edges = -1000:0.1:1000;
        histogram(obj.data.lick_timestamps_adjusted, edges); title('Time distribution of licks'); xlabel('Time (s)');

    elseif(Args.TrialLickDistribution)
        % Trial Lick Distribution
        trial = n;
        edges = -1000:0.01:1000;
        histogram(obj.data.lick_timestamps_adjusted(obj.data.lick_timestamps_spliced(:,2) == trial), edges); title('Time distribution of licks'); xlabel('Time (s)');

    elseif(Args.LickBinned)
        % Lick Binned
        histogram(obj.data.lick_timestamps_spliced(:,3), 1:100); title('Bin distribution of licks'); xlabel('Bin No.');
        xlim([1,100]);
        
    elseif(Args.TrialLickBinned)
        % Trial Lick Binned
        trial = n;
        histogram(obj.data.lick_timestamps_spliced(obj.data.lick_timestamps_spliced(:,2) == trial,3), 1:100); title('Bin distribution of licks'); xlabel('Bin No.');
        xlim([1,100]);

    elseif(Args.LickRate)
        % Lick Rate
        plot(1:size(obj.data.lick_binned,1), obj.data.lick_binned(:,3), 'b'); title('Lick Rate along Track Distance'); xlabel('Distance (Bin)'); ylabel('Lick Rate (s^-1)');

    elseif(Args.TrialLickRate)
        % Trial Lick Rate
        trial = n;
        lick_rate = obj.data.lick_count_binned(:,trial) / obj.data.data_bin(obj.data.data_bin(:,2) == trial,6);
        plot(1:size(obj.data.lick_binned,1), lick_rate, 'b'); title('Lick Rate along Track Distance'); xlabel('Distance (Bin)'); ylabel('Lick Rate (s^-1)');

    elseif(Args.TrialLick)   
        % Trial Lick
        trial = n;
        plot_start_idx = obj.data.TrialTime_idx(trial,1);
        plot_end_idx = obj.data.TrialTime_idx(trial,2);
        plot(obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,7)*220, 'b'); title('Distance'); xlabel('Time (s)'); ylabel('Distance (cm)'); xline(obj.data.data_trial(trial,3), 'r');
        hold on
        plot(obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), obj.data.lick_count_vel_filt(plot_start_idx:plot_end_idx,4)*250, 'g');
        hold off

    elseif(Args.LickRZ)   
        % Lick Reward Zone
        lick_reward_zone_prc = zeros(100,1);
        rz_timing = 5000; % Change reward zone timing
        for i = 1:100
            lick_reward_zone_prc(i) = nnz(obj.data.lick_timestamps_spliced(find(obj.data.lick_timestamps_spliced(:,3) == i & abs(obj.data.lick_timestamps_adjusted(:,1)) <= rz_timing), [1])) / nnz(obj.data.lick_timestamps_adjusted(:,1));
        end

        plot(lick_reward_zone_prc); title('Lick % (Reward Zone) < 5 sec'); xlabel('Bins'); ylabel('No. of Licks (%)');

    elseif(Args.LickBurst)
        % Lick Bursts
        subplot(2,1,1);
        plot(obj.data.session_data_exclude_zero_trials(:,1), obj.data.session_data_exclude_zero_trials(:,5)*2, 'b'); xlabel('Time (s)');
        hold on
        plot(obj.data.session_data_exclude_zero_trials(:,1), obj.data.session_data_exclude_zero_trials(:,10), 'r'); title('Licks'); xlabel('Time (s)');
        ax = gca;
        ax.YAxis.Visible = 'off';
        hold off
        
        subplot(2,1,2);
        %yyaxis left
        plot(obj.data.session_data_exclude_zero_trials(:,1), obj.data.session_data_exclude_zero_trials(:,7)*100, 'b'); title('Lick Bursts'); xlabel('Time (s)'); ylabel('Distance (Bin)'); xline(obj.data.data_trial(:,3), 'r');
        hold on
        plot(obj.data.session_data_exclude_zero_trials(:,1), obj.data.lick_burst(:)*80, 'k-');
%         yyaxis right
%         plot(obj.data.session_data_exclude_zero_trials(:,1), obj.data.velocity_averaged(:,4)); title('Velocity'); xlabel('Time (s)'); ylabel('Velocity (cm/s)');
        hold off
        
    elseif(Args.TrialLickBurst)
        % Trial Lick Frequency
        trial = n;
        plot_start_idx = obj.data.TrialTime_idx(trial,1);
        plot_end_idx = obj.data.TrialTime_idx(trial,2);
                
        subplot(5,1,1);
        plot(obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,5)*2, 'r'); xlabel('Time (s)');
        hold on
        plot(obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,6), 'k'); title('Licks'); xlabel('Time (s)');
        ax = gca;
        ax.YAxis.Visible = 'off';
        hold off
        
        subplot(5,1,2);
        plot(obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,5)*30, 'r'); xlabel('Time (s)');
        hold on
        plot(obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), obj.data.lick_freq(plot_start_idx:plot_end_idx), 'k'); title('Lick Frequency'); xlabel('Time (s)');ylabel('Hertz (Hz)');
        hold off
        
        subplot(5,1,3);
        yyaxis left
        plot(obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,7)*100, 'b'); title('Distance'); xlabel('Time (s)'); ylabel('Distance (Bin)'); xline(obj.data.data_trial(trial,3), 'r');
        hold on
        plot(obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), obj.data.lick_count(plot_start_idx:plot_end_idx,4)*100, 'k');
        yyaxis right
        plot(obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), obj.data.velocity_averaged(plot_start_idx:plot_end_idx,4)); title('Velocity'); xlabel('Time (s)'); ylabel('Velocity (cm/s)');
        hold off
        
        subplot(5,1,4);
        yyaxis left
        plot(obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,7)*100, 'b'); title('Distance'); xlabel('Time (s)'); ylabel('Distance (Bin)'); xline(obj.data.data_trial(trial,3), 'r');
        hold on
        yyaxis right
        plot(obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), obj.data.lick_freq(plot_start_idx:plot_end_idx), 'k'); title('Lick Frequency'); xlabel('Time (s)');ylabel('Hertz (Hz)');
        hold off
        
        subplot(5,1,5);
        yyaxis left
        plot(obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,7)*100, 'b'); title('Distance'); xlabel('Time (s)'); ylabel('Distance (Bin)'); xline(obj.data.data_trial(trial,3), 'r');
        hold on
        plot(obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), obj.data.lick_burst(plot_start_idx:plot_end_idx)*80, 'k-');
        yyaxis right
        plot(obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), obj.data.velocity_averaged(plot_start_idx:plot_end_idx,4)); title('Lick Burst'); xlabel('Time (s)'); ylabel('Velocity (cm/s)');
        hold off
        
    elseif(Args.AnticipatoryLick)
        
        trial = n;
        lick_timestamps_adjusted = obj.data.lick_timestamps_adjusted;
        
        edges = -10:0.5:10.5;
        [N,edges] = histcounts(lick_timestamps_adjusted, edges, 'Normalization', 'probability');
        [M,I] = max(N);
        plot(edges(1:end-1),N,'k-');
        xline(0,'k--');
        ylabel('Lick Density'); xlabel('Time (s)');
        hold on
        x_space = edges(1:end-1);
        b0 = [0,edges(I),1];
        lb = [0,-5,0];
        ub = [1,5,5];
        
%         if ~isnan(N)
%             [fitobject,gof] = fit(x_space', N', 'gauss1', ...
%                 'StartPoint', b0, ...
%                 'Lower', lb, ...
%                 'Upper', ub, ...
%                 'Robust','off');
%             y_pred = fitobject.a1*exp(-((x_space-fitobject.b1)/fitobject.c1).^2);
%             plot(x_space,y_pred,'r--');
%             text(0.9*min(x_space), 0.9*max(N), {"Mean: " + fitobject.b1 + ", StDev: " + fitobject.c1, "rmse: " + gof.rmse},'FontSize',6);
%         end
        
    elseif(Args.SummaryPlot)
        
        trial = n;
        plot_start_idx = obj.data.TrialTime_idx(trial,1);
        plot_end_idx = obj.data.TrialTime_idx(trial,2);
        subplot(3,1,1);
        plot(obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,7)*220);
        title('Distance','FontSize',20); xlabel('Time (s)','FontSize',18); ylabel('Distance (cm)','FontSize',18); xline(obj.data.data_trial(trial,3));
        subplot(3,1,2);
        plot(obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), obj.data.velocity_averaged(plot_start_idx:plot_end_idx,4));
        title('Velocity','FontSize',20); xlabel('Time (s)','FontSize',18); ylabel('Velocity (cm/s)','FontSize',18); xline(obj.data.data_trial(trial,3));
        subplot(3,1,3);
        plot(obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,1), obj.data.session_data_exclude_zero_trials(plot_start_idx:plot_end_idx,5)*2, 'b'); xlabel('Time (s)');
        
    elseif(Args.TrialdFF0)
        trial = n;
        trial = 9;
        plot_start_idx = find(obj.data.tsF > obj.data.data_trial(trial,2) + obj.data.actual_start_time,1,'first') + 350;
        plot_end_idx = find(obj.data.tsF > obj.data.data_trial(trial,3) + obj.data.actual_start_time,1,'first') - 1 - 60;
        
        plot_idx_sess = obj.data.tsFindex(plot_start_idx:plot_end_idx);
        
        tbl = timetable(seconds(obj.data.tsF(plot_start_idx:plot_end_idx)));
        
        varname = {'Distance (cm)'};
        tbl = [tbl table(obj.data.session_data_raw(plot_idx_sess,7)*220,'VariableNames',varname)];
        
        varname = {'Velocity (cm/s)'};
        tbl = [tbl table(obj.data.session_data_raw(plot_idx_sess,10),'VariableNames',varname)];
        
%         neuron_list = [1:20];
        neuron_list = [163 175 180 309 121 466 19 2 1 286 400];
%         neuron_list = [((n-1)*20)+1:((n-1)*20)+20];
        for neuron_no = 1:size(neuron_list,2)
            varname = {convertStringsToChars('n' + string(neuron_list(neuron_no)))};
            temp_neural_response = obj.data.dF_F0_corrected(neuron_list(neuron_no),plot_start_idx:plot_end_idx)';
            temp_neural_response = normalize(temp_neural_response,"range");
            tbl = [tbl table(temp_neural_response,'VariableNames',varname)];
        end
        
        s = stackedplot(tbl,"FontSize",14);
        s.LineProperties(1).Color = [0 0 0];
        s.LineProperties(2).Color = [0.8500 0.3250 0.0980];
        for neuron_no = 1:size(neuron_list,2)
            s.DisplayLabels(neuron_no+2) = {''};
        end
    
        %for neuron_no = 1:5
%             subplot(5,1,neuron_no);
%             plot(obj.data.tsF(plot_start_idx:plot_end_idx), obj.data.dF_F0_corrected(neuron_no,plot_start_idx:plot_end_idx));
%             if neuron_no == 1
%                 title('Fluorescence Trace','FontSize',20);
%             end
%             if neuron_no == 3
%                 ylabel('Fluorescence Intensity (a.u.)','FontSize',18);
%             end
%             if neuron_no == 5
%                 xlabel('Time (s)','FontSize',18); 
%             end
%         end
        
    elseif(Args.PopVec)
        
        bin_no = n;
        delete(findall(gcf,'type','annotation'));
        
        session_data_raw_dsp = obj.data.session_data_raw(obj.data.tsFindex,:);
        bin_dsp_idx = find(session_data_raw_dsp(:,3) == bin_no);
        velocity_bin_dsp = session_data_raw_dsp(bin_dsp_idx,10);
        dF_bin = obj.data.dF_F0_corrected(:,bin_dsp_idx);
        dF_bin = dF_bin(:,velocity_bin_dsp >= 5); % Velocity filter
%         cellfastcombinedData = mtcellfastcombined('auto').data;
        cellfastcombinedData = obj.data.cellfastcombinedData;
        
        decoded_bins = cellfastcombinedData.BayeDecoding_shuffle.pos_decoded_bin_shuffle(cellfastcombinedData.BayeDecoding_shuffle.bin_actual == bin_no);
        
% %         KLDIV_matrix = zeros(cellfastcombinedData.data.nNeuron,cellfastcombinedData.data.nTrials);
%         max_magnitude = zeros(cellfastcombinedData.data.nNeuron,cellfastcombinedData.data.nTrials);
%         for trial_no = 1:cellfastcombinedData.data.nTrials
%             plot_start_idx = find(obj.data.tsF > obj.data.data_trial(trial_no,2) + obj.data.actual_start_time,1,'first');
%             plot_end_idx = find(obj.data.tsF > obj.data.data_trial(trial_no,3) + obj.data.actual_start_time,1,'first') - 1;
%             if isempty(plot_end_idx)
%                 plot_end_idx = size(obj.data.tsF,1);
%             end
%             
%             for neuron_no = 1:cellfastcombinedData.data.nNeuron
%                 temp_dFF0 = obj.data.dF_F0_corrected(neuron_no,plot_start_idx:plot_end_idx);
%                 plot_idx_sess = obj.data.tsFindex(plot_start_idx:plot_end_idx);
%                 dFF0_bin_idx = obj.data.session_data_raw(plot_idx_sess,3) == bin_no;
%                 dFF0_bin = temp_dFF0(dFF0_bin_idx);
%                 if ~isempty(dFF0_bin)
%                     max_magnitude(neuron_no,trial_no) = max(dFF0_bin);
%                 else
%                     max_magnitude(neuron_no,trial_no) = NaN;
%                 end
%                 
% %                 ideal_square_impulse = temp*max(temp_dFF0);
% %                 KLDIV_matrix(neuron_no,trial_no) = KLDiv((temp_dFF0+eps),(ideal_square_impulse+eps)');
% %                 KLDiv((temp_dFF0+eps))')
%                 
% %                 figure; plot(temp_dFF0); hold on; plot(ideal_square_impulse,'red');
%             end
%         end
%         
% %         KLDIV_matrix_trialAveraged = mean(KLDIV_matrix,2,'omitnan');
%         max_magnitude_trialAveraged = mean(max_magnitude,2,'omitnan');
        
        popvector = zeros(cellfastcombinedData.nNeuron,3);
        binToAngle = linspace(3.6,360,100) - 1.8;
        for cell_no = 1:cellfastcombinedData.nNeuron
            cellFiringRate = squeeze(cellfastcombinedData.binFiringRate(cell_no,:,:));
            [~,maxidx] = max(mean(cellFiringRate));
            popvector(cell_no,1) = maxidx; % Preferred angle based on mean binned firing rate
            popvector(cell_no,2) = binToAngle(maxidx); % Preferred angle in polar coordinates
            popvector(cell_no,3) = mean(dF_bin(cell_no,:),'omitnan'); % dFF0 magnitude
%             popvector(cell_no,3) = mean(cellFiringRate(:,bin_no),'omitnan');
        end
        
%         popvector(:,3) = max_magnitude_trialAveraged; % dFF0 magnitude
%         popvector = popvector(popvector(:,3) > mean(popvector(:,3)),:);
%         popvector = popvector(popvector(:,3) > prctile(popvector(:,3),75),:);
%         popvector = popvector(cellfastcombinedData.data.FieldNo' ~= 0,:);
        

%         popvector(:,4) = temp2.FieldNo';
%         test = sortrows(popvector,3);
%         figure; plot(test(:,3),smooth(test(:,2)));
%         xline(prctile(test(:,3),25),'r--'); xline(prctile(test(:,3),50),'r');
%         xline(prctile(test(:,3),75),'r--'); xline(prctile(test(:,3),90),'r--');
%         yline(binToAngle(bin_no),'b--'); yline(sum(popvector(:,2).*popvector(:,3)) / sum(popvector(:,3)),'b-'); 
        
        polarplot(deg2rad(popvector(:,2)),popvector(:,3),'o')
        actual_theta = deg2rad(binToAngle(bin_no));
        centre = 0.518;
        radius = 0.2;
        x = [centre centre+radius*sin(actual_theta)];
        y = [centre centre+radius*cos(actual_theta)];
        ta1 = annotation('textarrow',x,y,'String',' Actual ','FontSize',13,'Linewidth',2,'Color','g')
%         sum(popvector(:,2).*popvector(:,3)) / sum(popvector(:,3))
%         decoded_theta = deg2rad(sum(popvector(:,2).*popvector(:,3)) / sum(popvector(:,3)));
        decoded_theta = deg2rad(binToAngle(round(circular_mean(decoded_bins))));
        x = [centre centre+radius*sin(decoded_theta)];
        y = [centre centre+radius*cos(decoded_theta)];
        ta2 = annotation('textarrow',x,y,'String',' Decoded ','FontSize',13,'Linewidth',2,'Color','r')
        title('Position (cm)','FontSize',18);
        
        pax = gca;
        pax.ThetaDir = 'clockwise';
        pax.ThetaZeroLocation = 'top';
        pax.ThetaTickLabel = {'0'; '18.3'; '36.7'; '55'; '73.3'; '91.7'; '110'; '128.3'; '146.7'; '165'; '183.3'; '201.7'};
        
    elseif(Args.FRVel)
        
        cell_no = n;
        
        velocity_averaged_dsp = obj.data.session_data_raw(obj.data.tsFindex,10);
        cell_dFF0 = obj.data.dF_F0_corrected(cell_no,:)';
        
        hold on
        [R,Pvalue] = corr(velocity_averaged_dsp,cell_dFF0,'Type','Pearson')
        plot(velocity_averaged_dsp,cell_dFF0,'k.','MarkerSize',2);
%         xlim([0,2.5]); ylim([0 55]);
        f = fit(velocity_averaged_dsp,cell_dFF0,'poly1');
        plot(f,'k');
        xlabel("Velocity (cm/s)",'FontSize',18); ylabel("dFF0",'FontSize',18);
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
