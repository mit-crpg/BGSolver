%% savefiledlg
% File saving dialog function.
%
% fname = savefiledlg([ext,extname,dlgtitle])
%  or
% fname = savefiledlg([{ext,extname},dlgtitle])
%
% This function acts as the MATLAB uiputfile function, but throws an error if
% the user clicks 'Cancel'. Same as uiputfile, it optionally accepts a file
% extension, a name for the file extension, and a dialog title.
%
% Package:    BGSolver v1.03
% Subpackage: Utilities
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Master file saving dialog function

function fname = savefiledlg(varargin)
% The purpose of this function is to optionally accept a file extension, a name
% for the file extension and a dialog title, and request the file name to save
% from the user. An error is thrown if the user clicks 'Cancel'.
%
% INPUTS/OUTPUTS
% ext - File extension. Default: All MATLAB files.
% extname - Name for the file extension. Default: File extension itself.
% dlgtitle - Dialog title. Default: 'Select File to Write'.
% fname - Full file name to save.

% Constructing input parser
ip = buildIP();
% Parsing input
ip.parse(varargin{:});
% Retrieving inputs
if iscellstr(ip.Results.arg1)
   ext = ip.Results.arg1{1};
   extname = ip.Results.arg1{2};
   dlgtitle = ip.Results.arg2;
else
   ext = ip.Results.arg1;
   extname = ip.Results.arg2;
   dlgtitle = ip.Results.arg3;
end;
% Clearing memory
clear('ip');

% Requesting file to save
if isempty(ext) && isempty(extname) && isempty(dlgtitle)
   [ShortFName,DirPath] = uiputfile();
elseif ~isempty(ext) && isempty(extname) && isempty(dlgtitle)
   [ShortFName,DirPath] = uiputfile(ext);
elseif ~isempty(ext) && ~isempty(extname) && isempty(dlgtitle)
   [ShortFName,DirPath] = uiputfile({ext,extname});
elseif ~isempty(ext) && ~isempty(extname) && ~isempty(dlgtitle)
   [ShortFName,DirPath] = uiputfile({ext,extname},dlgtitle);
else % Should be impossible to reach
   error('Input parsing error!!!');
end;

% Checking if user clicked 'Cancel'
assert(~all(ShortFName == 0) && ~all(DirPath == 0),...
   'File name not specified!!!');
fname = fullfile(DirPath,ShortFName);
end

%% Input parser construction function

function ip = buildIP()
% The purpose of this function is to construct an input parser object for the
% savefiledlg function.
%
% INPUTS/OUTPUTS
% ip - Input parser object for the savefiledlg function.

ip = inputParser;
ip.addOptional('arg1',[],@(arg1)(iscellstr(arg1) || ischar(arg1)));
ip.addOptional('arg2',[],@ischar);
ip.addOptional('arg3',[],@ischar);
end