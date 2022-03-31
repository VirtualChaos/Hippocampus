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

Args = struct('ObjectLevel',0, 'AnalysisLevel',0, 'FiringRateMapRaw',0, 'FiringRateMapAdSm',0, ...
              'BinFiringRate',0, 'TrialBinFiringRate',0, 'Baseline',0, 'GMM',0, 'TrialGMM',0, 'TrialFits',0,  'AlphaAdSm', 0);
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
elseif (Args.FiringRateMapRaw | Args.FiringRateMapAdSm | Args.BinFiringRate | Args.Baseline | Args.GMM | Args.AlphaAdSm)
    % r = length(fieldnames(obj.data.cellData));
    r = obj.data.nNeuron;
elseif (Args.TrialBinFiringRate | Args.TrialGMM | Args.TrialFits)
    r = obj.data.nTrials*obj.data.nNeuron;
else
	% if we don't recognize and of the options, pass the call to parent
	% in case it is to get number of events, which has to go all the way
	% nptdata/get
	r = get(obj.nptdata,varargin{:});
end