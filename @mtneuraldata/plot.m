function [obj, varargout] = plot(obj,varargin)
%@dirfiles/plot Plot function for dirfiles object.
%   OBJ = plot(OBJ) creates a raster plot of the neuronal
%   response.

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		  'ReturnVars',{''}, 'ArgsOnly',0, 'Cmds','', ...
          'FiringRateMapRaw',0, 'FiringRateMapAdSm',0, 'BinFiringRate',0, ...
          'Baseline',0, 'GMM',0);
Args.flags = {'LabelsOff','ArgsOnly','FiringRateMapRaw','FiringRateMapAdSm','BinFiringRate', ...
              'Baseline','GMM'};
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
        fns = fieldnames(obj.data.cellData);

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
        fns = fieldnames(obj.data.cellData);

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
        % Adaptive Smoothed Firing Rate Map
        cell_no = n;
        fns = fieldnames(obj.data.cellData);
        
        imagesc(obj.data.cellData.(fns{cell_no}).binFiringRate); title('Trial Firing Rates'); xlabel('Position Bins'); ylabel('Trials');
        colorbar;
    elseif(Args.Baseline)
        % Baseline Firing
        cell_no = n;
        fns = fieldnames(obj.data.cellData);
        
        plot(1:100, obj.data.placefieldData.(fns{cell_no}).mapLsm_FFT_Low_Pass, 'b'); legends{1} = sprintf('Firing Rate Map');
        hold on
        plot(obj.data.placefieldData.(fns{cell_no}).baseline(27:126, 1), obj.data.placefieldData.(fns{cell_no}).baseline(27:126, 2),'r'); legends{2} = sprintf('Estimated Baseline');
        legend(legends);
        hold off
        
    elseif(Args.GMM)
        % Gaussian Mixed Model
        cell_no = n;
        fns = fieldnames(obj.data.cellData);
        
        hold all
        title('Firing Rate Map'); xlabel('Position Bins'); ylabel('Firing Rate (spike/sec)');
        plot(obj.data.placefieldData.(fns{cell_no}).basemapLrw, 'g'); legends{1} = sprintf('Raw');
        
        plot(obj.data.placefieldData.(fns{cell_no}).basemapLsm, 'b'); legends{2} = sprintf('AdSmooth');
        
        plot(obj.data.placefieldData.(fns{cell_no}).mapLsm_FFT_Low_Pass, 'r'); legends{3} = sprintf('AdSm-FFT-Low-Pass');
        
        if ~isempty(obj.data.placefieldData.(fns{cell_no}).GMM)
            mu_lines = [obj.data.placefieldData.(fns{cell_no}).GMM(:,1) [1:size(obj.data.placefieldData.(fns{cell_no}).GMM,1)]'];
            mu_lines_labels = cellstr("Field " + num2str(mu_lines(:,2)));
            xline(mu_lines(:,1), '--', mu_lines_labels, 'Color', '#D95319');
            ax = gca;
            for peak_no = 1:size(obj.data.placefieldData.(fns{cell_no}).GMM,1)
                left = obj.data.placefieldData.(fns{cell_no}).GMM(peak_no,1) - obj.data.placefieldData.(fns{cell_no}).GMM(peak_no,2);
                right = obj.data.placefieldData.(fns{cell_no}).GMM(peak_no,1) + obj.data.placefieldData.(fns{cell_no}).GMM(peak_no,2);
                x = [left right right left];
                y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];
                patch(x, y, 'k', 'FaceAlpha', '0.1', 'LineStyle', 'none');
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
            yline((baseline_mean_height + peak_height_from_baseline/2),'r');
            yline(std_threshold,'b');
            
            x = linspace(-25,126,152);
            if size(obj.data.placefieldData.(fns{cell_no}).GMM,1) == 1
                y = gaus(x,obj.data.placefieldData.(fns{cell_no}).GMM(1,1),obj.data.placefieldData.(fns{cell_no}).GMM(1,2),obj.data.placefieldData.(fns{cell_no}).GMM(1,3),vo);
            elseif size(obj.data.placefieldData.(fns{cell_no}).GMM,1) == 2
                y = gaus(x,obj.data.placefieldData.(fns{cell_no}).GMM(1,1),obj.data.placefieldData.(fns{cell_no}).GMM(1,2),obj.data.placefieldData.(fns{cell_no}).GMM(1,3),vo) + gaus(x,obj.data.placefieldData.(fns{cell_no}).GMM(2,1),obj.data.placefieldData.(fns{cell_no}).GMM(2,2),obj.data.placefieldData.(fns{cell_no}).GMM(2,3),vo);
            else
                y = gaus(x,obj.data.placefieldData.(fns{cell_no}).GMM(1,1),obj.data.placefieldData.(fns{cell_no}).GMM(1,2),obj.data.placefieldData.(fns{cell_no}).GMM(1,3),vo) + gaus(x,obj.data.placefieldData.(fns{cell_no}).GMM(2,1),obj.data.placefieldData.(fns{cell_no}).GMM(2,2),obj.data.placefieldData.(fns{cell_no}).GMM(2,3),vo) + gaus(x,obj.data.placefieldData.(fns{cell_no}).GMM(3,1),obj.data.placefieldData.(fns{cell_no}).GMM(3,2),obj.data.placefieldData.(fns{cell_no}).GMM(3,3),vo);
            end
            
            y = 0;
            for i = 1:size(obj.data.placefieldData.(fns{cell_no}).GMM,1)
                y = gaus(x,obj.data.placefieldData.(fns{cell_no}).GMM(i,1),obj.data.placefieldData.(fns{cell_no}).GMM(i,2),obj.data.placefieldData.(fns{cell_no}).GMM(i,3),vo) + y;
            end
            
            legends{4} = ''; legends{5} = ''; legends{6} = ''; legends{7} = ''; legends{8} = ''; legends{9} = '';
            plot(x(27:126), y(27:126), 'k-', 'LineWidth',3); legends{10} = sprintf('GMM');
            legend(legends, 'Location', 'northeastoutside');
        end
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
