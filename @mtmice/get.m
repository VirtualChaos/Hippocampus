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

Args = struct('ObjectLevel',0, 'AnalysisLevel',0, 'MatchedFiringRateMap',0,'MatchedCorrMatrix',0,'PlaceFieldPositionEntropy',0,...
                'AdSm',0,'Lick',0);
Args.flags = {'ObjectLevel','AnalysisLevel'};
Args = getOptArgs(varargin,Args);

% set variables to default
r = [];

if(Args.ObjectLevel)
	% specifies that the object should be created in the session directory
	r = levelConvert('levelNo',1);
elseif(Args.AnalysisLevel)
	% specifies that the AnalysisLevel of the object is 'AllIntragroup'
	r = 'Single';
elseif (Args.MatchedFiringRateMap | Args.MatchedCorrMatrix)
    r = size(obj.data.roiMatchData.allSessionMapping,1);
elseif(Args.Lick)
    r = size(fieldnames(obj.data.sessionCombined),1);
else
	% if we don't recognize and of the options, pass the call to parent
	% in case it is to get number of events, which has to go all the way
	% nptdata/get
	r = get(obj.nptdata,varargin{:});
end