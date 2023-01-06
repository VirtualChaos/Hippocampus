function [r,varargout] = get(obj,varargin)
%dirfiles/get Get function for dirfiles objects
%dirfiles/GET Returns object properties
%   VALUE = GET(OBJ,PROP_NAME) returns an object 
%   property.
%   In dirfiles, PROP_NAME can be one of the following:
%      'ObjectLevel'
%	 'AnalysisLevel'
%
%   Dependencies: 

Args = struct('ObjectLevel',0, 'AnalysisLevel',0, 'MeanVel',0, 'MedianVel',0, 'TrialVel',0, 'MiceVelDist',0, 'LickPrc',0, ...
              'MiceLickBurstWidth',0, 'LickBurstWidth',0, 'VelLickPrc',0, 'MiceVelLickPrc',0,'VelStateEstimate',0);
Args.flags ={'ObjectLevel','AnalysisLevel'};
Args = getOptArgs(varargin,Args);

% set variables to default
r = [];

if(Args.ObjectLevel)
	% specifies that the object should be created in the session directory
	r = levelConvert('levelNo',1);
elseif(Args.AnalysisLevel)
	% specifies that the AnalysisLevel of the object is 'AllIntragroup'
	r = 'Single';
elseif (Args.MeanVel | Args.MedianVel | Args.TrialVel | Args.LickPrc | Args.VelLickPrc)
    r = length(fieldnames(obj.data.sessionCombined));
elseif (Args.MiceVelDist | Args.MiceVelLickPrc | Args.MiceLickBurstWidth | Args.VelStateEstimate)
    r = obj.data.nMice;
else
	% if we don't recognize and of the options, pass the call to parent
	% in case it is to get number of events, which has to go all the way
	% nptdata/get
	r = get(obj.nptdata,varargin{:});
end