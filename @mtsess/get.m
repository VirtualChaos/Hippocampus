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

Args = struct('ObjectLevel',0, 'AnalysisLevel',0, 'VelRaw',0, 'TrialVelRaw', 0, 'TrialVel',0, 'Vel',0, ...
              'TrialVelFilt',0, 'VelBinned',0, 'VelCount',0, 'WaterLick',0, 'TrialWaterLick',0, ...
              'LickDistribution',0, 'TrialLickDistribution',0, 'LickBinned',0, 'TrialLickBinned',0, ...
              'LickRate',0, 'TrialLickRate',0, 'TrialLick',0, 'LickRZ', 0);
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
elseif (Args.TrialVelRaw | Args.TrialVel | Args.TrialVelFilt | Args.TrialWaterLick | ...
        Args.TrialLickDistribution | Args.TrialLickBinned | Args.TrialLickRate | Args.TrialLick)
    r = obj.data.nTrials;
else
	% if we don't recognize and of the options, pass the call to parent
	% in case it is to get number of events, which has to go all the way
	% nptdata/get
	r = get(obj.nptdata,varargin{:});
end