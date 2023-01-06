function [obj, varargout] = plot(obj,varargin)
%@dirfiles/plot Plot function for dirfiles object.
%   OBJ = plot(OBJ) creates a raster plot of the neuronal
%   response.

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		  'ReturnVars',{''}, 'ArgsOnly',0, 'Cmds','', ...
          'BinFiringRate',0,'FiringMap',0);
Args.flags = {'LabelsOff','ArgsOnly','BinFiringRate','FiringMap'};
[Args,varargin2] = getOptArgs(varargin,Args);

% if user select 'ArgsOnly', return only Args structure for an empty object
if Args.ArgsOnly
    Args = rmfield (Args, 'ArgsOnly');
    varargout{1} = {'Args',Args};
    return;
end

if(~isempty(Args.NumericArguments))
	% plot one data set at a time
	n = Args.NumericArguments{1};
	if(Args.BinFiringRate)
        
        subplot(2,1,1);
        imagesc(obj.data.binFiringRate); title('Firing Rate'); xlabel('Distance (Bin)'); ylabel('Trial');
        
        subplot(2,1,2);
        plot(obj.data.maps_raw, 'g'); legends{1} = sprintf('Raw');
        hold on
        plot(obj.data.maps_adsm, 'b'); legends{2} = sprintf('AdSmooth');
        plot(imgaussfilt(obj.data.maps_raw, 3),'r'); legends{3} = sprintf('Gaussian (sigma = 3)');
        legend(legends);
        
	elseif(Args.FiringMap)

        subplot(3,2,1);
        imagesc(obj.data.maps_raw); title('Raw'); xlabel('Distance (Bin)');
        colorbar;
        subplot(3,2,3);
        imagesc(obj.data.maps_raw1); title('Raw (1st half)'); xlabel('Distance (Bin)');
        colorbar;
        subplot(3,2,5);
        imagesc(obj.data.maps_raw2); title('Raw (2nd half)'); xlabel('Distance (Bin)');
        colorbar;
        
        subplot(3,2,2);
        imagesc(obj.data.maps_adsm);  title('AdSm'); xlabel('Distance (Bin)');
        colorbar;
        subplot(3,2,4);
        imagesc(obj.data.maps_adsm1); title('AdSm (1st half)'); xlabel('Distance (Bin)');
        colorbar;
        subplot(3,2,6);
        imagesc(obj.data.maps_adsm2); title('AdSm (2nd half)'); xlabel('Distance (Bin)');
        colorbar;
        
	else
		% code to plot yet another kind of plot
		
		% label the axis
		xlabel('X Axis')
		ylabel('Y Axis')
	end	

	% add an appropriate title
	sdstr = get(obj,'SessionDirs');
	title(getDataOrder('ShortName','DirString',sdstr{1}))
else
	% plot all data
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
