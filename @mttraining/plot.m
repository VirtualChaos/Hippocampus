function [obj, varargout] = plot(obj,varargin)
%@dirfiles/plot Plot function for dirfiles object.
%   OBJ = plot(OBJ) creates plots for the mice group

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		  'ReturnVars',{''}, 'ArgsOnly',0, 'Cmds','', ...
          'MeanVel',0, 'MedianVel',0, 'TrialVel',0, 'MiceVelDist',0, 'Vel',0, 'MiceVelBinned',0, 'LickPrc',0, 'LickBurstWidth',0, 'MiceLickBurstWidth',0, ...
          'MiceLickPosition',0,'LickTime',0,'MiceLickTime',0,'VelLickPrc',0, 'MiceVelLickPrc',0,'VelStateEstimate',0);
Args.flags = {'LabelsOff','ArgsOnly','MeanVel','MedianVel','TrialVel','MiceVelDist','Vel','MiceVelBinned','LickPrc', 'LickBurstWidth', ...
          'MiceLickBurstWidth','MiceLickPosition','LickTime','MiceLickTime','VelLickPrc','MiceVelLickPrc','VelStateEstimate'};
[Args,varargin2] = getOptArgs(varargin,Args);

% if user select 'ArgsOnly', return only Args structure for an empty object
if Args.ArgsOnly
    Args = rmfield (Args, 'ArgsOnly');
    varargout{1} = {'Args',Args};
    return;
end

if(~isempty(Args.NumericArguments))
	n = Args.NumericArguments{1};
	if(Args.MeanVel)
        fns = fieldnames(obj.data.sessionCombined);
        group_id = fns{n};
        
        groupData = obj.data.sessionCombined.(group_id);
        fns = fieldnames(groupData);
        for mice_no = 1:length(fieldnames(groupData))
            subplot(length(fieldnames(groupData)),1,mice_no)
            mean_vel = groupData.(fns{mice_no}).average_track_velocity(:,1);
            mean_vel_smooth = smooth(groupData.(fns{mice_no}).average_track_velocity(:,1),0.1);
            
            plot(mean_vel,'k.'); xlabel('Sessions'); ylabel('Velocity (cm/s)');
            hold on
            plot(mean_vel_smooth,'b-');
            title(fns{mice_no});
            yline(mean(mean_vel,'omitnan'),'g')
            yline(mean(mean_vel(ceil(length(mean_vel)/2):end),'omitnan'),'r') 
        end
        
    elseif(Args.MedianVel)
        fns = fieldnames(obj.data.sessionCombined);
        group_id = fns{n};
        
        groupData = obj.data.sessionCombined.(group_id);
        fns = fieldnames(groupData);
        for mice_no = 1:length(fieldnames(groupData))
            subplot(length(fieldnames(groupData)),1,mice_no)
            median_vel = groupData.(fns{mice_no}).average_velocity(:,2);
            plot(median_vel); xlabel('Sessions'); ylabel('Velocity (cm/s)');
            title(fns{mice_no});
            yline(median(median_vel,'omitnan'),'g')
            yline(mean(median_vel,'omitnan'),'r')      
        end
        
    elseif(Args.TrialVel)
        fns = fieldnames(obj.data.sessionCombined);
        group_id = fns{n};
        
        groupData = obj.data.sessionCombined.(group_id);
        fns = fieldnames(groupData);
        for mice_no = 1:length(fieldnames(groupData))
            subplot(length(fieldnames(groupData)),1,mice_no)
            trial_vel = groupData.(fns{mice_no}).average_trial_velocity(:,2);
            plot(trial_vel); xlabel('Sessions'); ylabel('Velocity (cm/s)');
            title(fns{mice_no});
            yline(median(trial_vel,'omitnan'),'g')
            yline(mean(trial_vel,'omitnan'),'r')      
        end
        
    elseif(Args.MiceVelDist)
        fns = fieldnames(obj.data.sessionCombined);
        mice_id = [];
        for group_no = 1:4
            temp1 = repmat(group_no,length(fieldnames(obj.data.sessionCombined.(fns{group_no}))),1);
            temp2 = [1:length(fieldnames(obj.data.sessionCombined.(fns{group_no})))]';
            mice_id = [mice_id; [temp1 temp2]];
        end
        
        group_id = (fns{mice_id(n,1)});
        groupData = obj.data.sessionCombined.(group_id);
        
        fns = fieldnames(groupData);
        miceData = groupData.(fns{mice_id(n,2)});
        sess_list = fieldnames(miceData.data_trial);
        
        index = reshape(1:ceil(length(fieldnames(miceData.data_trial))/5)*5,5,[]).';
        for sess_no = 1:length(fieldnames(miceData.data_trial))
            subplot(ceil(length(fieldnames(miceData.data_trial))/5),5,index(sess_no))
            sessionData.velocity_binned = miceData.velocity_binned.(sess_list{sess_no});
            
%             edges = 0:0.5:ceil(max(sessionData.velocity_averaged(:,4)));
%             histogram(sessionData.velocity_averaged(sessionData.velocity_averaged(:,4) ~= 0,4),edges); title(fns{mice_id(n,2)} + " - D" + miceData.sessionDays(sess_no));% ylabel('Velocity (cm/s)');
%             xline(mean(sessionData.velocity_averaged(sessionData.velocity_averaged(:,4) ~= 0,4), 'omitnan'),'LineWidth',2);
            
            edges = 0:0.5:ceil(max(sessionData.velocity_binned));
            histogram(sessionData.velocity_binned,edges); title(fns{mice_id(n,2)} + " - D" + miceData.sessionDays(sess_no));
            xline(mean(sessionData.velocity_binned, 'omitnan'),'LineWidth',2); xlim([0 30]);
            
            %             [y_data,edges] = histcounts(sessionData.velocity_averaged(sessionData.velocity_averaged(:,4) ~= 0,4),edges,'Normalization','pdf'); title(fns{mice_id(n,2)} + " - D" + miceData.sessionDays(sess_no));
            %             edges = edges(2:end) - (edges(2)-edges(1))/2;
            %             plot(edges,y_data);
            %             xlim([0 30]);
            %             xline(mean(sessionData.velocity_averaged(sessionData.velocity_averaged(:,4) ~= 0,4), 'omitnan'));
            
%             edges = 0:0.5:8;
%             histogram(220./sessionData.data_trial(:,4),edges); title(fns{mice_id(n,2)} + " - D" + miceData.sessionDays(sess_no));% ylabel('Velocity (cm/s)');
%             xline(mean(220./sessionData.data_trial(:,4), 'omitnan'),'LineWidth',2);
        end
        
    elseif(Args.MiceVelBinned)
        fns = fieldnames(obj.data.sessionCombined);
        mice_id = [];
        for group_no = 1:4
            temp1 = repmat(group_no,length(fieldnames(obj.data.sessionCombined.(fns{group_no}))),1);
            temp2 = [1:length(fieldnames(obj.data.sessionCombined.(fns{group_no})))]';
            mice_id = [mice_id; [temp1 temp2]];
        end
        
        group_id = (fns{mice_id(n,1)});
        groupData = obj.data.sessionCombined.(group_id);
        
        fns = fieldnames(groupData);
        miceData = groupData.(fns{mice_id(n,2)});
        sess_list = fieldnames(miceData.data_trial);
        
        mice_no = char(fns(mice_id(n,2)));
        mice_no = str2num(mice_no(3:4));
        sessionDate = obj.data.sessionDate(find(obj.data.sessionDate(:,1) == mice_no),:);
        
        index = reshape(1:ceil(length(fieldnames(miceData.data_trial))/5)*5,5,[]).';
        for sess_no = 1:length(fieldnames(miceData.data_trial))
            subplot(ceil(length(fieldnames(miceData.data_trial))/5),5,index(sess_no))
            velocity_binned = miceData.velocity_binned.(sess_list{sess_no});
            if size(velocity_binned,2) > 1 % Average across trials
                velocity_binned = mean(velocity_binned,1,'omitnan')';
            end
            velocity_binned = [velocity_binned(51:end); velocity_binned(1:50)];
                        
            %imagesc(velocity_binned'); colormap hot; colorbar
            plot(velocity_binned');
            xline(51,'r--'); yline(mean(velocity_binned),'b--')
            title(fns{mice_id(n,2)} + " - D" + miceData.sessionDays(sess_no) + "; " + sessionDate(sess_no,2));   
        end
    
        elseif(Args.Vel)
        
            for group_no = 1:4
                %             ha = subplot(2,2,group_no)
                switch group_no
                    case 1
                        ha = subplot(2,2,2)
                    case 2
                        ha = subplot(2,2,1)
                    case 3
                        ha = subplot(2,2,4)
                    case 4
                        ha = subplot(2,2,3)
                end
                fns = fieldnames(obj.data.sessionCombined);
                group_id = fns{group_no};
                
                groupData = obj.data.sessionCombined.(group_id);
                fns2 = fieldnames(groupData);
                y_data_combined_mice = [];
                nTrials = [];
                for mice_no = 1:length(fieldnames(groupData))
                    miceData = groupData.(fns2{mice_no});
                    sess_list = fieldnames(miceData.data_trial);
                    
                    mean_vel = miceData.average_track_velocity(:,1);
                    vel_binned_combined = zeros(length(sess_list),1);
                    
                    padSize = floor(5/2); % 5 days sliding window
                    mean_vel_padded = padarray(mean_vel, padSize, 'replicate');
                    mean_vel_window_smooth = zeros(size(mean_vel,1), 1);
                    for i = 1+padSize:size(mean_vel, 1)+padSize
                        mean_vel_window_smooth(i - padSize) = mean(mean_vel_padded((i - padSize):(i + padSize)));
                    end
                    
                    for sess_no = 1:length(fieldnames(miceData.data_trial))
                        sessionData.velocity_binned = miceData.velocity_binned.(sess_list{sess_no});
                        vel_binned_combined(sess_no) = mean(sessionData.velocity_binned(51:end),'omitnan');
                    end
                    upper_vel_thres = mean(vel_binned_combined(ceil(length(vel_binned_combined)/2):end),'omitnan');
                    lower_vel_thres = min(vel_binned_combined(1:ceil(length(vel_binned_combined)/2)),[],'omitnan');
                    
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
                    
                    mice_velocity_binned = zeros(length(sess_list),100);
                    for sess_no = 1:length(fieldnames(miceData.data_trial))
                        sessionData.velocity_binned = miceData.velocity_binned.(sess_list{sess_no});
                        
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
                        
                        %observation = mean(sessionData.velocity_binned(51:end),'omitnan');
                        %observation = mean(bootstrap_vel,'omitnan');
                        observation = mean_vel_window_smooth(sess_no);
                        [gm_posterior,nlogL] = posterior(gm,observation);
                        %[gm_posterior,nlogL] = posterior(gm,gmPDF2);
                        
                        covNames = {'diagonal','full'};
                        CovType = find(strncmpi(gm.CovType,covNames,length(gm.CovType)));
                        nlog_lh = wdensity_dupe(observation,gm.mu, gm.Sigma, gm.PComponents, gm.SharedCov, CovType) * -1; %nLogL for individual components
                        
                        decision_state(sess_no,:) = [gm_posterior nlogL nlog_lh];
                        
                        velocity_binned = miceData.velocity_binned.(sess_list{sess_no});
                        if size(velocity_binned,2) > 1 % Average across trials
                            velocity_binned = mean(velocity_binned,1,'omitnan')';
                        end
                        mice_velocity_binned(sess_no,:) = velocity_binned;
                    end
                    
                    trained_sess_idx = find(decision_state(:,2) > 0.95,1,'first');
                    y_data_combined_mice = [y_data_combined_mice; mean(mice_velocity_binned(trained_sess_idx:end,:),'omitnan')];
                    nTrials = [nTrials; sum(groupData.(fns2{mice_no}).nTrials(trained_sess_idx:end),'omitnan')];
                end
                adjusted_y_data_combined_mice = [y_data_combined_mice(:,51:end) y_data_combined_mice(:,1:50)];
                
                adjusted_y_data_combined_group = mean(adjusted_y_data_combined_mice,'omitnan');
                bar([-49:50],adjusted_y_data_combined_group);
                title(group_id);
                hold on
                trial_no = sum(nTrials,'omitnan');
                
                %             f = fit([-49:50]',adjusted_y_data_combined_group','gauss1');
                %             plot(f,'k-');
                %             delete(findobj('type','legend'))
                %
                %             text(5, 1.3*max(adjusted_y_data_combined_group), {"Trials: " + trial_no, "Peak: " + f.b1 + ", StDev: " + f.c1},'FontSize',12);
                %             xlabel('Position (Bins)'); ylabel('Density');
                %             xlim([-49 50]); ylim(ha,[0 1.5*max(adjusted_y_data_combined_group)]);
                
                %yyaxis right
                ft = fittype('piecewiseLine(x,a,b,k)');
                f = fit([-49:50]',adjusted_y_data_combined_group',ft,'StartPoint',[1,1,0],'Lower',[1,1,0],'Upper',[Inf,Inf,Inf])
                mu = f.a * f.b;
                sigma = f.a * (f.b)^2;
                ylim([0 21])
                
            end
            
    elseif(Args.LickPrc)
        fns = fieldnames(obj.data.sessionCombined);
        group_id = fns{n};
        
        groupData = obj.data.sessionCombined.(group_id);
        fns = fieldnames(groupData);
        for mice_no = 1:length(fieldnames(groupData))
            subplot(length(fieldnames(groupData)),1,mice_no)
            lick_prc_RZ = groupData.(fns{mice_no}).lick_prc_RZ;
            plot(lick_prc_RZ(:,3)*100,'rx'); xlabel('Sessions'); ylabel('RZ Lick (%)');
            title(fns{mice_no});    
        end    
        
    elseif(Args.LickBurstWidth)
        
        for group_no = 1:4
%             ha = subplot(2,2,group_no)
            switch group_no
                case 1
                    ha = subplot(2,2,2)
                case 2
                    ha = subplot(2,2,1)
                case 3
                    ha = subplot(2,2,4)
                case 4
                    ha = subplot(2,2,3)
            end
            fns = fieldnames(obj.data.sessionCombined);
            group_id = fns{group_no};
            
            groupData = obj.data.sessionCombined.(group_id);
            fns2 = fieldnames(groupData);
            y_data_combined_mice = [];
            nTrials = [];
            for mice_no = 1:length(fieldnames(groupData))
                miceData = groupData.(fns2{mice_no});
                sess_list = fieldnames(miceData.data_trial);
                
                mean_vel = miceData.average_track_velocity(:,1);
                vel_binned_combined = zeros(length(sess_list),1);
                
                padSize = floor(5/2); % 5 days sliding window
                mean_vel_padded = padarray(mean_vel, padSize, 'replicate');
                mean_vel_window_smooth = zeros(size(mean_vel,1), 1);
                for i = 1+padSize:size(mean_vel, 1)+padSize
                    mean_vel_window_smooth(i - padSize) = mean(mean_vel_padded((i - padSize):(i + padSize)));
                end
                
                for sess_no = 1:length(fieldnames(miceData.data_trial))
                    sessionData.velocity_binned = miceData.velocity_binned.(sess_list{sess_no});
                    vel_binned_combined(sess_no) = mean(sessionData.velocity_binned(51:end),'omitnan');
                end
                upper_vel_thres = mean(vel_binned_combined(ceil(length(vel_binned_combined)/2):end),'omitnan');
                lower_vel_thres = min(vel_binned_combined(1:ceil(length(vel_binned_combined)/2)),[],'omitnan');
                
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
                
                for sess_no = 1:length(fieldnames(miceData.data_trial))
                    sessionData.velocity_binned = miceData.velocity_binned.(sess_list{sess_no});
                    
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
                                        
                    %observation = mean(sessionData.velocity_binned(51:end),'omitnan');
                    %observation = mean(bootstrap_vel,'omitnan');
                    observation = mean_vel_window_smooth(sess_no);
                    [gm_posterior,nlogL] = posterior(gm,observation);
                    %[gm_posterior,nlogL] = posterior(gm,gmPDF2);
                    
                    covNames = {'diagonal','full'};
                    CovType = find(strncmpi(gm.CovType,covNames,length(gm.CovType)));
                    nlog_lh = wdensity_dupe(observation,gm.mu, gm.Sigma, gm.PComponents, gm.SharedCov, CovType) * -1; %nLogL for individual components
                    
                    decision_state(sess_no,:) = [gm_posterior nlogL nlog_lh];
                    
                end
                
                trained_sess_idx = find(decision_state(:,2) > 0.95,1,'first');
                y_data_combined_mice = [y_data_combined_mice; mean(groupData.(fns2{mice_no}).y_data_combined(trained_sess_idx:end,:),'omitnan')];
                nTrials = [nTrials; sum(groupData.(fns2{mice_no}).nTrials(trained_sess_idx:end),'omitnan')];
            end
            adjusted_y_data_combined_mice = [y_data_combined_mice(:,51:end) y_data_combined_mice(:,1:50)];
            
            adjusted_y_data_combined_group = mean(adjusted_y_data_combined_mice,'omitnan');
            bar([-49:50],adjusted_y_data_combined_group);
            title(group_id);
            hold on
            trial_no = sum(nTrials,'omitnan');
            
%             f = fit([-49:50]',adjusted_y_data_combined_group','gauss1');
%             plot(f,'k-');
%             delete(findobj('type','legend'))
% 
%             text(5, 1.3*max(adjusted_y_data_combined_group), {"Trials: " + trial_no, "Peak: " + f.b1 + ", StDev: " + f.c1},'FontSize',12);
%             xlabel('Position (Bins)'); ylabel('Density');
%             xlim([-49 50]); ylim(ha,[0 1.5*max(adjusted_y_data_combined_group)]);
                        
            %yyaxis right
            ft = fittype('piecewiseLine(x,a,b,k)');
            f = fit([-49:50]',adjusted_y_data_combined_group',ft,'StartPoint',[1,1,0],'Lower',[1,1,0],'Upper',[Inf,Inf,Inf])
            mu = f.a * f.b;
            sigma = f.a * (f.b)^2;
            ylim([0 0.25]);
        
        end
        
        
    elseif(Args.MiceLickBurstWidth)
        
        hold off
        fns = fieldnames(obj.data.sessionCombined);
        mice_id = [];
        for group_no = 1:4
            temp1 = repmat(group_no,length(fieldnames(obj.data.sessionCombined.(fns{group_no}))),1);
            temp2 = [1:length(fieldnames(obj.data.sessionCombined.(fns{group_no})))]';
            mice_id = [mice_id; [temp1 temp2]];
        end
        
        group_id = (fns{mice_id(n,1)});
        groupData = obj.data.sessionCombined.(group_id);
        
        fns = fieldnames(groupData);
        miceData = groupData.(fns{mice_id(n,2)});
        sess_list = fieldnames(miceData.data_trial);
        
        mean_vel = miceData.average_track_velocity(:,1);
        
        padSize = floor(5/2); % 5 days sliding window
        mean_vel_padded = padarray(mean_vel, padSize, 'replicate');
        mean_vel_window_smooth = zeros(size(mean_vel,1), 1);
        for i = 1+padSize:size(mean_vel, 1)+padSize
            mean_vel_window_smooth(i - padSize) = mean(mean_vel_padded((i - padSize):(i + padSize)));
        end
        
        upper_vel_thres = mean(mean_vel_window_smooth(ceil(length(mean_vel_window_smooth)/2):end),'omitnan');
        lower_vel_thres = min(mean_vel_window_smooth(1:ceil(length(mean_vel_window_smooth)/2)),[],'omitnan');
        
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
        
        for sess_no = 1:length(fieldnames(miceData.data_trial))
            sessionData.velocity_binned = miceData.velocity_binned.(sess_list{sess_no});
            
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

            %observation = mean(bootstrap_vel,'omitnan');
            %observation = mean_vel(sess_no);
            observation = mean_vel_window_smooth(sess_no);
            [gm_posterior,nlogL] = posterior(gm,observation);
            %[gm_posterior,nlogL] = posterior(gm,gmPDF2);
            
            covNames = {'diagonal','full'};
            CovType = find(strncmpi(gm.CovType,covNames,length(gm.CovType)));
            nlog_lh = wdensity_dupe(observation,gm.mu, gm.Sigma, gm.PComponents, gm.SharedCov, CovType) * -1; %nLogL for individual components
            
            decision_state(sess_no,:) = [gm_posterior nlogL nlog_lh];
        end
        
        trained_sess_idx = find(decision_state(:,2) > 0.95,1,'first');
        if isempty(trained_sess_idx) % VelStateEstimate could not identify the post-trained point
            trained_sess_idx = ceil(0.80 * length(fieldnames(miceData.data_trial))); % Assume mice was trained in the last 20% of recorded sessions
        end
        %y_data_combined = mean(miceData.y_data_combined(trained_sess_idx:end,:),'omitnan');
        y_data_combined = mean(miceData.lick_prc(trained_sess_idx:end,:),'omitnan') ./ 100;
        nTrials = sum(miceData.nTrials(trained_sess_idx:end),'omitnan');
    
        adjusted_y_data_combined = [y_data_combined(51:end) y_data_combined(1:50)];
                       
        %bar([-49:50],adjusted_y_data_combined);
        x_space = linspace(0.05,5,100);
        %x_space = linspace(0.02,2,100);
        
        y_data_combined_baseline_corr = y_data_combined;
%         baseline = 0;
        baseline = prctile(y_data_combined,50);
%         y_data_combined_baseline_corr = y_data_combined - baseline;
%         y_data_combined_baseline_corr(y_data_combined_baseline_corr < 0) = 0;
               
        subplot(1,2,1)
        x_space = linspace(1,100,100);
        %bar(x_space,y_data_combined_baseline_corr);
        bar(x_space,y_data_combined);
        title(fns{mice_id(n,2)});
        hold on
        trial_no = sum(nTrials,'omitnan');
        xlim([0 max(x_space)]); ylim([0 1.3*max(y_data_combined_baseline_corr)]);
        
%         b0 = [5,0.8];
%         lb = [0,0.8];
%         ub = [10,0.8];
        
%         b0 = [5,2];
%         lb = [0,0];
%         ub = [10,5];
        
        b0 = [3,2];
        lb = [3,0];
        ub = [3,5];
        
        fun = @(k,theta,x) ((1 / gamma(k) .* (theta^k)) .* (x.^(k-1)) .* exp(-(x./theta)));
        [fitobject,gof] = fit(x_space', y_data_combined_baseline_corr', fun, ...
                          'StartPoint', b0, ...
                          'Lower', lb, ...
                          'Upper', ub, ...
                          'Robust','LAR');
        
        plot(x_space,gampdf(x_space,fitobject.k,fitobject.theta)+baseline);
        text(0.2*max(x_space), 1.1*max(y_data_combined), {"Trials: " + trial_no, "k: " + fitobject.k + ", theta: " + fitobject.theta, "rmse: " + gof.rmse},'FontSize',12);
        xlabel('Position (Bins)'); ylabel('Density');
        xlim([0 max(x_space)]); ylim([0 1.3*max(y_data_combined)]);
        
        y_data_combined_baseline_corr = adjusted_y_data_combined;
%         y_data_combined_baseline_corr = adjusted_y_data_combined - baseline;
%         y_data_combined_baseline_corr(y_data_combined_baseline_corr < 0) = 0;
        
        subplot(1,2,2)
        x_space = linspace(1,100,100);
        %bar(x_space,y_data_combined_baseline_corr);
        bar(x_space,adjusted_y_data_combined);
        title(fns{mice_id(n,2)});
        hold on
        trial_no = sum(nTrials,'omitnan');
        xlim([0 max(x_space)]); ylim([0 1.3*max(y_data_combined_baseline_corr)]);
        
        b0 = [0,50,10];
        lb = [0,30,0];
        ub = [20,70,50];

        [fitobject,gof] = fit(x_space', y_data_combined_baseline_corr', 'gauss1', ...
                          'StartPoint', b0, ...
                          'Lower', lb, ...
                          'Upper', ub, ...
                          'Robust','off');
        
        y_pred = fitobject.a1*exp(-((x_space-fitobject.b1)/fitobject.c1).^2);
        plot(x_space,y_pred+baseline);
        text(0.2*max(x_space), 1.1*max(y_data_combined), {"Trials: " + trial_no, "Mean: " + fitobject.b1 + ", StDev: " + fitobject.c1, "rmse: " + gof.rmse},'FontSize',12);
        xlabel('Position (Bins)'); ylabel('Density');
        xlim([0 max(x_space)]); ylim([0 1.3*max(y_data_combined)]);
    
    elseif(Args.MiceLickPosition)
        fns = fieldnames(obj.data.sessionCombined);
        mice_id = [];
        for group_no = 1:4
            temp1 = repmat(group_no,length(fieldnames(obj.data.sessionCombined.(fns{group_no}))),1);
            temp2 = [1:length(fieldnames(obj.data.sessionCombined.(fns{group_no})))]';
            mice_id = [mice_id; [temp1 temp2]];
        end
        
        group_id = (fns{mice_id(n,1)});
        groupData = obj.data.sessionCombined.(group_id);
        
        fns = fieldnames(groupData);
        miceData = groupData.(fns{mice_id(n,2)});
        sess_list = fieldnames(miceData.data_trial);
        
        mice_no = char(fns(mice_id(n,2)));
        mice_no = str2num(mice_no(3:4));
        sessionDate = obj.data.sessionDate(find(obj.data.sessionDate(:,1) == mice_no),:);
        
        index = reshape(1:ceil(length(fieldnames(miceData.data_trial))/5)*5,5,[]).';
        for sess_no = 1:length(fieldnames(miceData.data_trial))
            subplot(ceil(length(fieldnames(miceData.data_trial))/5),5,index(sess_no))
            lick_count_binned = miceData.lick_count_binned.(sess_list{sess_no});
            
            lick_count_binned = [zeros(1,100); lick_count_binned;]; % Pad one trial (100 bins) before
            lick_count_binned = reshape(lick_count_binned.',1,[]); % Flatten array
            lick_count_binned = circshift(lick_count_binned,-50); % Shift array 50 bins behind
            lick_count_binned = reshape(lick_count_binned,100,miceData.nTrials(sess_no)+1)'; % Restore array
            lick_binned_combined = sum(lick_count_binned,1);
            % imagesc(lick_count_binned); colorbar;
            plot(lick_binned_combined ./ sum(lick_binned_combined));
            
            title(fns{mice_id(n,2)} + " - D" + miceData.sessionDays(sess_no) + "; " + sessionDate(sess_no,2));
        end
    
    elseif(Args.LickTime)
        trained_sessionDate = load('/Volumes/HippocampusNew/NTU/Training_data/trained_sessionDate.mat').trained_sessionDate;
        
        fns = fieldnames(obj.data.sessionCombined);
        mice_id = [];
        for group_no = 1:4
            temp1 = repmat(group_no,length(fieldnames(obj.data.sessionCombined.(fns{group_no}))),1);
            temp2 = [1:length(fieldnames(obj.data.sessionCombined.(fns{group_no})))]';
            mice_id = [mice_id; [temp1 temp2]];
        end
        
        group_id = (fns{n});
        groupData = obj.data.sessionCombined.(group_id);
        
        fns2 = fieldnames(groupData);
        for mice_idx = 1:length(fns2)
            subplot(ceil(length(fns2)/3),3,mice_idx)
            miceData = groupData.(fns2{mice_id(mice_idx,2)});
            sess_list = fieldnames(miceData.data_trial);
            
            mice_no = char(fns2(mice_id(mice_idx,2)));
            mice_no = str2num(mice_no(3:4));
            sessionDate = obj.data.sessionDate(find(obj.data.sessionDate(:,1) == mice_no),:);
            trained_sess_idx = trained_sessionDate(trained_sessionDate(:,3) == mice_no,4);
            
            edges = -10:0.5:10.5;
            N_combined = [];
            for sess_no = trained_sess_idx:length(fieldnames(miceData.data_trial))
                lick_timestamps_adjusted = miceData.lick_timestamps_adjusted.(sess_list{sess_no});
                
                [N,edges] = histcounts(lick_timestamps_adjusted, edges);
                N_combined = [N_combined; N];
            end
            N_combined = sum(N_combined);
            plot(edges(1:end-1),N_combined,'k-');
            xline(0,'k--');
            ylabel('Lick count'); xlabel('Time (s)');
            male_mice_id = [59 60 61 64 65 55 56 58 62 63 53 54 68 69 83 84 85];
            if ismember(mice_no,male_mice_id) % Male
                title(fns2{mice_id(mice_idx,2)});
            else % Female
                title(fns2{mice_id(mice_idx,2)},'Color','r');
            end
        end
        
    elseif(Args.MiceLickTime)
        trained_sessionDate = load('/Volumes/HippocampusNew/NTU/Training_data/trained_sessionDate.mat').trained_sessionDate;
        
        fns = fieldnames(obj.data.sessionCombined);
        mice_id = [];
        for group_no = 1:4
            temp1 = repmat(group_no,length(fieldnames(obj.data.sessionCombined.(fns{group_no}))),1);
            temp2 = [1:length(fieldnames(obj.data.sessionCombined.(fns{group_no})))]';
            mice_id = [mice_id; [temp1 temp2]];
        end
        
        group_id = (fns{mice_id(n,1)});
        groupData = obj.data.sessionCombined.(group_id);
        
        fns = fieldnames(groupData);
        miceData = groupData.(fns{mice_id(n,2)});
        sess_list = fieldnames(miceData.data_trial);
        
        mice_no = char(fns(mice_id(n,2)));
        mice_no = str2num(mice_no(3:4));
        sessionDate = obj.data.sessionDate(find(obj.data.sessionDate(:,1) == mice_no),:);
        trained_sess_idx = trained_sessionDate(trained_sessionDate(:,3) == mice_no,4);
        index = reshape(1:ceil(length(fieldnames(miceData.data_trial))/5)*5,5,[]).';
        for sess_no = 1:length(fieldnames(miceData.data_trial))
            subplot(ceil(length(fieldnames(miceData.data_trial))/5),5,index(sess_no))
            lick_timestamps_adjusted = miceData.lick_timestamps_adjusted.(sess_list{sess_no});
            
            edges = -10:0.5:10.5;
            [N,edges] = histcounts(lick_timestamps_adjusted, edges, 'Normalization', 'probability');
            [M,I] = max(N);
            plot(edges(1:end-1),N,'k-');
            xline(0,'k--');
            ylabel('Lick count PDF'); xlabel('Time (s)');
            if sess_no < trained_sess_idx % Pre-trained
                title(fns{mice_id(n,2)} + " - D" + miceData.sessionDays(sess_no) + "; " + sessionDate(sess_no,2));
            else % Trained
                title(fns{mice_id(n,2)} + " - D" + miceData.sessionDays(sess_no) + "; " + sessionDate(sess_no,2),'Color','r');
            end
            hold on
            x_space = edges(1:end-1);
            b0 = [0,edges(I),1];
            lb = [0,-5,0];
            ub = [1,5,5];
            
            if ~isnan(N)
                [fitobject,gof] = fit(x_space', N', 'gauss1', ...
                    'StartPoint', b0, ...
                    'Lower', lb, ...
                    'Upper', ub, ...
                    'Robust','off');
                y_pred = fitobject.a1*exp(-((x_space-fitobject.b1)/fitobject.c1).^2);
                plot(x_space,y_pred,'r--');
                text(0.9*min(x_space), 0.9*max(N), {"Mean: " + fitobject.b1 + ", StDev: " + fitobject.c1, "rmse: " + gof.rmse},'FontSize',6);
            end
            
        end
        
%         edges = -10:0.5:10.5;
%         N_combined = [];
%         for sess_no = trained_sess_idx:length(fieldnames(miceData.data_trial))
%             lick_timestamps_adjusted = miceData.lick_timestamps_adjusted.(sess_list{sess_no});
% 
%             [N,edges] = histcounts(lick_timestamps_adjusted, edges);
%             N_combined = [N_combined; N];
%         end
%         N_combined = sum(N_combined);
%         plot(edges(1:end-1),N_combined,'k-');
%         xline(0,'k--');
%         ylabel('Lick count'); xlabel('Time (s)');
%         title(fns{mice_id(n,2)});
            
    elseif(Args.VelLickPrc)
        fns = fieldnames(obj.data.sessionCombined);
        group_id = fns{n};
        
        groupData = obj.data.sessionCombined.(group_id);
        fns = fieldnames(groupData);
        for mice_no = 1:length(fieldnames(groupData))
            subplot(length(fieldnames(groupData)),1,mice_no)
            mean_vel = groupData.(fns{mice_no}).average_velocity(:,1);
            mean_vel_smooth = smooth(groupData.(fns{mice_no}).average_velocity(:,1),0.1);
            
            yyaxis left
            plot(mean_vel,'b.'); xlabel('Sessions'); ylabel('Velocity (cm/s)');
            hold on
            plot(mean_vel_smooth,'b-');
            title(fns{mice_no});
            yline(median(mean_vel,'omitnan'),'g')
            %yline(mean(mean_vel,'omitnan'),'r')     
            yyaxis right
            lick_prc_RZ = groupData.(fns{mice_no}).lick_prc_RZ;
            lick_prc_RZ = lick_prc_RZ * 100;
            plot(lick_prc_RZ(:,3),'rx'); xlabel('Sessions'); ylabel('RZ Lick (%)');
            hold on
            
            % fit to a*(1-exp(-b.*x)) model;
            expc1 = fittype('a*(1-exp(-b.*x))','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b'});% a*(1-exp(-b.*x))
            [cf,gof] = fit([1:size(lick_prc_RZ(:,3),1)]',lick_prc_RZ(:,3),expc1,'Lower',[0,0],'Upper',[Inf,Inf]);
            p1=plot(cf,'r-'); ylim([0,120]); p1.LineWidth=1; %xlim([0,22]);
        end
        
    elseif(Args.MiceVelLickPrc)
        
        hold off
        fns = fieldnames(obj.data.sessionCombined);
        mice_id = [];
        for group_no = 1:4
            temp1 = repmat(group_no,length(fieldnames(obj.data.sessionCombined.(fns{group_no}))),1);
            temp2 = [1:length(fieldnames(obj.data.sessionCombined.(fns{group_no})))]';
            mice_id = [mice_id; [temp1 temp2]];
        end
        
        group_id = (fns{mice_id(n,1)});
        groupData = obj.data.sessionCombined.(group_id);
        
        fns = fieldnames(groupData);
        miceData = groupData.(fns{mice_id(n,2)});
        sess_list = fieldnames(miceData.data_trial);
        
        mean_vel = miceData.average_track_velocity(:,1);
        mean_vel_smooth = smooth(miceData.average_track_velocity(:,1),0.1);
        vel_binned_combined = zeros(length(sess_list),1);
        vel_binned_metric = zeros(length(sess_list),1);
        
        padSize = floor(5/2); % 5 days sliding window
        mean_vel_padded = padarray(mean_vel, padSize, 'replicate');
        mean_vel_window_smooth = zeros(size(mean_vel,1), 1);
        for i = 1+padSize:size(mean_vel, 1)+padSize
            mean_vel_window_smooth(i - padSize) = mean(mean_vel_padded((i - padSize):(i + padSize)));
        end
        
        for sess_no = 1:length(fieldnames(miceData.data_trial))
            sessionData.velocity_binned = miceData.velocity_binned.(sess_list{sess_no});
            vel_binned_combined(sess_no) = mean(sessionData.velocity_binned(51:end),'omitnan');
            if size(sessionData.velocity_binned,2) > 1 % Average across trials
                sessionData.velocity_binned = mean(sessionData.velocity_binned,1,'omitnan')';
            end
            velocity_binned = [sessionData.velocity_binned(51:end); sessionData.velocity_binned(1:50)];
            vel_binned_metric(sess_no) = mean(velocity_binned([[1:50] [54:end]])) - mean(velocity_binned(51:53));
        end
        
        %vel_metric = mean_vel;
        vel_metric = vel_binned_metric;
        
        upper_vel_thres = mean(vel_metric(ceil(length(vel_metric)/2):end),'omitnan');
        %lower_vel_thres = mean_vel(1);
        lower_vel_thres = min(vel_metric(1:ceil(length(vel_metric)/2)),[],'omitnan');
        
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
        for sess_no = 1:length(fieldnames(miceData.data_trial))
            %subplot(ceil(length(fieldnames(miceData.data_trial))/5),5,index(sess_no))
            sessionData.velocity_binned = miceData.velocity_binned.(sess_list{sess_no});
                       
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
            observation = vel_metric(sess_no);
            [gm_posterior,nlogL] = posterior(gm,observation);
            %[gm_posterior,nlogL] = posterior(gm,gmPDF2);
            
            covNames = {'diagonal','full'};
            CovType = find(strncmpi(gm.CovType,covNames,length(gm.CovType)));
            nlog_lh = wdensity_dupe(observation,gm.mu, gm.Sigma, gm.PComponents, gm.SharedCov, CovType) * -1; %nLogL for individual components
            
            decision_state(sess_no,:) = [gm_posterior nlogL nlog_lh];
                        
        end

        vel_metric2 = mean_vel_window_smooth;
        upper_vel_thres2 = mean(vel_metric2(ceil(length(vel_metric2)/2):end),'omitnan');
        lower_vel_thres2 = min(vel_metric2(1:ceil(length(vel_metric2)/2)),[],'omitnan');
        
        %index = reshape(1:ceil(length(fieldnames(miceData.data_trial))/5)*5,5,[]).';
        for sess_no = 1:length(fieldnames(miceData.data_trial))
            %subplot(ceil(length(fieldnames(miceData.data_trial))/5),5,index(sess_no))
            sessionData.velocity_binned = miceData.velocity_binned.(sess_list{sess_no});
            
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
            %observation = mean_vel_window_smooth(sess_no);
            observation = vel_metric2(sess_no);
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
        plot(vel_metric,'b.'); xlabel('Sessions'); ylabel('Velocity (cm/s)');
        hold on
        %plot(mean_vel_window_smooth,'rx'); xlabel('Sessions'); ylabel('Velocity (cm/s)');
        title(fns{mice_id(n,2)});
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
        lick_prc_RZ = miceData.lick_prc_RZ;
        plot(lick_prc_RZ(:,3),'rx'); xlabel('Sessions'); ylabel('RZ Lick (%)');
        hold on
        
        % fit to a*(1-exp(-b.*x)) model;
        expc1 = fittype('a*(1-exp(-b.*x))','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b'});% a*(1-exp(-b.*x))
        [cf,gof] = fit([1:size(lick_prc_RZ(:,3),1)]',lick_prc_RZ(:,3),expc1,'Lower',[0,0],'Upper',[Inf,Inf]);
        plot(cf,'r-'); ylim([0,110]); %p1.LineWidth=1; %xlim([0,22]);
        hold off
        delete(findobj('type','legend'))
        yline(0.75*cf(length(sess_list)),'r--');
        
    elseif(Args.VelStateEstimate)
        
        fns = fieldnames(obj.data.sessionCombined);
        mice_id = [];
        for group_no = 1:4
            temp1 = repmat(group_no,length(fieldnames(obj.data.sessionCombined.(fns{group_no}))),1);
            temp2 = [1:length(fieldnames(obj.data.sessionCombined.(fns{group_no})))]';
            mice_id = [mice_id; [temp1 temp2]];
        end
        
        group_id = (fns{mice_id(n,1)});
        groupData = obj.data.sessionCombined.(group_id);
        
        fns = fieldnames(groupData);
        miceData = groupData.(fns{mice_id(n,2)});
        sess_list = fieldnames(miceData.data_trial);
        
        mean_vel = miceData.average_track_velocity(:,1);
        upper_vel_thres = mean(mean_vel(ceil(length(mean_vel)/2):end),'omitnan');
        lower_vel_thres = mean_vel(1);

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
        
        decision_state = zeros(length(sess_list),3);

        index = reshape(1:ceil(length(fieldnames(miceData.data_trial))/5)*5,5,[]).';
        for sess_no = 1:length(fieldnames(miceData.data_trial))
            subplot(ceil(length(fieldnames(miceData.data_trial))/5),5,index(sess_no))
            sessionData.velocity_binned = miceData.velocity_binned.(sess_list{sess_no});
                       
            if sess_no > 1
                mu(sess_no,:) = [lower_vel_thres upper_vel_thres];
                sigma(sess_no,:) = cat(3,[1],[1]);
                mix_prop(sess_no,:) = decision_state(sess_no-1,1:2);
                mix_prop(sess_no,1) = mix_prop(sess_no,1) / (2^decision_state(sess_no-1,3)); % Scaling factor - Pressure to switch to high velocity state
                if mix_prop(sess_no,1) == 0
                    mix_prop(sess_no,1) = mix_prop(sess_no,1) + eps;
                end
            end
    
            % Obtain prior
            gm = gmdistribution(mu(sess_no,:)',sigma(sess_no,:,:,:),mix_prop(sess_no,:));
            gmPDF = pdf(gm,[0:0.1:40]');
            gmPDF = gmPDF / sum(gmPDF);
            plot(0:0.1:40,gmPDF,'b-');
            xline(gm.mu(1)); xline(gm.mu(2));
            title(fns{mice_id(n,2)} + " - D" + miceData.sessionDays(sess_no));
            xlabel('Velocity (cm/s)'); xlim([0 40]);
            hold on
            
            % Calculate likelihood
            bootstrap_vel = bootstrp(100,@median,sessionData.velocity_binned);
            [fi,xi] = ksdensity(bootstrap_vel);
            fi = fi ./ sum(fi);% / 20
            plot(xi,fi,'kx'); xlim([0 40]);
            hold on
            %f = fit(xi',fi','gauss1');
            %plot(f,'k-');
            GMModel = fitgmdist([xi]',1);
            gmPDF2 = pdf(GMModel,[0:0.1:40]');
            gmPDF2 = gmPDF2 ./ sum(gmPDF2);
            plot(0:0.1:40,gmPDF2,'r-');
            
            %[gm_posterior,nlogL] = posterior(gm,mean(sessionData.velocity_binned,'omitnan'));
            [gm_posterior,nlogL] = posterior(gm,mean(bootstrap_vel,'omitnan'));
            decision_state(sess_no,:) = [gm_posterior nlogL];
            %[gm_posterior,nlogL] = posterior(gm,gmPDF2);
        
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
