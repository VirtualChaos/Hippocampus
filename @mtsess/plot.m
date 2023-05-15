function [obj, varargout] = plot(obj,varargin)
%@dirfiles/plot Plot function for dirfiles object.
%   OBJ = plot(OBJ) creates plots for the session

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		  'ReturnVars',{''}, 'ArgsOnly',0, 'Cmds','', ...
          'VelRaw',0, 'TrialVelRaw',0, 'Vel',0, 'TrialVel',0, 'TrialVelFilt',0, 'VelBinned',0, ...
          'VelCount',0, 'WaterLick',0, 'TrialWaterLick',0, 'LickDistribution',0, 'TrialLickDistribution',0, ...
          'LickBinned',0, 'TrialLickBinned',0, 'LickRate',0, 'TrialLickRate',0, 'TrialLick',0, 'LickRZ',0, ...
          'LickBurst',0,'TrialLickBurst',0);
Args.flags = {'LabelsOff','ArgsOnly','VelRaw','TrialVelRaw','Vel', 'TrialVel', 'TrialVelFilt', ...
              'VelBinned', 'VelCount', 'WaterLick','TrialWaterLick', ...
              'TrialLickDistribution', 'LickDistribution', 'LickBinned', 'TrialLickBinned', ...
              'LickRate', 'TrialLickRate', 'TrialLick', 'LickRZ', 'LickBurst', 'TrialLickBurst'};
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
        plot(1:size(obj.data.velocity_binned(:,1),1), obj.data.velocity_binned(:,1), 'b'); title('Mean Velocity along Track Distance'); xlabel('Distance (Bin)'); ylabel('Mean Velocity (cm/s)');
 	
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
