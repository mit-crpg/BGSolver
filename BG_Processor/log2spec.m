%% log2spec
% Logging according to specification function.
%
% log2spec(str,strstat,fids,runopts)
%
% This function checks whether or not the processor is sufficiently verbose to
% log the string; if it is, the string is logged according to a log
% specification. If the string is logged, it is either output to the MATLAB
% command line, to a log file, or both.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Master logging function

function log2spec(varargin)
% The purpose of this function is to check whether the processor is sufficiently
% verbose to log the string; if it is, the string is logged according to a log
% specification. If the string is logged, it is either output to the MATLAB
% command line, to a log file, or both.
%
% INPUTS/OUTPUTS
% str - String to log according to specification.
% strstat - String verbosity status. numeric. Possible values:
%  1 - Major report.
%  2 - Minor report.
% fids - File IDs data structure. May include:
%  .logfid - Log file ID.
% runopts - Runtime and Recording Options data structure. Important fields:
%  .verbose - Processor's verbosity setting. int8. REQUIRED. Possible values:
%   0 - The processor is silent. Neither major nor minor steps are logged.
%   1 - The processor is neither silent nor verbose. Major steps are logged.
%   2 - The processor is verbose. Major and minor are logged.
%  .log - Processor's logging setting. int8. OPTIONAL. Possible values:
%   0 - Outputs are logged in the MATLAB command line.
%   1 - Outputs are logged in a log file.
%   2 - Outputs are logged in both the MATLAB command line and the log file.

% Input parser construction
ip = buildIP();
% Input parsing
ip.parse(varargin{:});
% Retrieving inputs
str = ip.Results.str;
strstat = int8(ip.Results.strstat);
fids = ip.Results.fids;
runopts = ip.Results.runopts;
assert(~(any(runopts.verbose == [1,2]) && any(runopts.log == [1,2]) && ...
      ~isfield(fids,'logfid')),'No log file ID specified!!!');
% Clearing memory
clear('ip');

if strstat <= runopts.verbose
   fmt = '%s';
   if strstat == 2
      fmt = ['\t' fmt];
   end;
   fmt = [fmt '\n'];
   
   if any(runopts.log == [0,2])
      fprintf(fmt,str);
   end;
   if any(runopts.log == [1,2])
      fprintf(fids.logfid,fmt,str);
   end;
end;
end

%% Input parser construction function

function ip = buildIP()
% The purpose of this function is to construct an input parser object for the
% log2spec function.
%
% INPUTS/OUTPUTS
% ip - Input parser object for the log2spec function.

fstrstat = @(strstat) any(strstat == [1,2]);
frunopts = @(ds) isstruct(ds) && (isfield(ds,'verbose') && ...
   ((ds.verbose == 0 && ~isfield(ds,'log')) || (any(ds.verbose == [1,2]) && ...
   (isfield(ds,'log') && any(ds.log == [0,1,2])))));

ip = inputParser;
ip.addRequired('str',@ischar);
ip.addRequired('strstat',fstrstat);
ip.addRequired('fids',@isstruct);
ip.addRequired('runopts',frunopts);
end