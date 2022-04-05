function [obj, varargout] = plot(obj,varargin)
%@dirfiles/plot Plot function for dirfiles object.
%   OBJ = plot(OBJ) creates plots for the individual mice

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		  'ReturnVars',{''}, 'ArgsOnly',0, 'Cmds','', ...
          'VelRaw',0, 'Vel',0, 'VelBinned',0, 'VelCount',0, 'WaterLick',0, 'LickDistribution',0, ...
          'LickBinned',0, 'LickRate',0, 'LickRZ', 0);
Args.flags = {'LabelsOff','ArgsOnly','VelRaw','Vel','VelBinned','VelCount','WaterLick','LickDistribution', ...
            'LickBinned','LickRate','LickRZ',};
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
        
    elseif(Args.Vel)
        
    elseif(Args.VelBinned)
        
    elseif(Args.VelCount)
        
    elseif(Args.WaterLick)
        
    elseif(Args.LickDistribution)
        
    elseif(Args.LickBinned)
        
    elseif(Args.LickRate)
        
    elseif(Args.LickRZ)
    
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
