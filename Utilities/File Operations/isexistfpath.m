%% isexistfpath
% File path existence checking function.
%
% ifexist = isexistfpath(fname)
%  or
% ifexist = isexistfpath({fname1,fname2...})
%
% This function checks one or more input strings for whether or not they are
% paths to existing files. A match to an existing object that is not a file
% (like a folder) returns a false.
%
% Package:    BGSolver v1.03
% Subpackage: Utilities
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Main file existence checking function

function ifexist = isexistfpath(varargin)
% The purpose of this function is to check one or more input strings for whether
% or not they are paths to existing files. A match to an existing object that is
% not a file (like a folder) returns a false.
%
% INPUTS/OUTPUTS
% fname - Possibly existing file name string.
% ifexist - Boolean array of whether or not the corresponding file exists.

% Input parser construction
ip = buildIP();
% Input parsing
ip.parse(varargin{:});
% Retrieving inputs
if iscellstr(ip.Results.arg) && ~isempty(ip.Results.arg)
   ifcell = true;
   fpaths = ip.Results.arg;
elseif ischar(ip.Results.arg) && ~isempty(ip.Results.arg)
   ifcell = false;
   fpath = ip.Results.arg;
else
   ifcell = false;
   fpath = '';
end;
% Clearing memory
clear('ip');

% Checking file existence
if ifcell
   f = @(arg) (exist(arg,'file') == 2);
   ifexist = cellfun(f,fpaths);
elseif ~isempty(fpath)
   ifexist = (exist(fpath,'file') == 2);
else
   ifexist = false;
end;
end

%% Input parser construction function

function ip = buildIP()
% The purpose of this function is to construct an input parser object for the
% isexistfpath function.
%
% INPUTS/OUTPUTS
% ip - Input parser object for the isexistfpath function.

ip = inputParser;
ip.addRequired('arg',@(arg)(iscellstr(arg) || ischar(arg)));
end