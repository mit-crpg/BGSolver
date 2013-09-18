%% mkdirwchk
% Directory creation function with existence check.
%
% dirpath = mkdirwchk(dirname,[dirrootpath])
%
% This function checks if a directory with a given name, and optionally a given
% root path, (absolute or relative), exists, and if it does not, it creates it.
% It returns the full path of this directory.
%
% Package:    BGSolver v1.03
% Subpackage: Utilities
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Master Directory Creation with Existence Check Function code

function dirpath = mkdirwchk(varargin)
% The purpose of this function is to check for a directory's existence, and if
% it does not, construct it.
%
% INPUTS/OUTPUTS
% dirname - Directory name.
% dirrootpath - Directory root path. OPTIONAL. If not given, the current
%               directory will be used as root path.
% dirpath - Full directory path.

% Input parser construction
ip = buildIP();
% Input parsing
ip.parse(varargin{:});
% Retrieving inputs
dirname = ip.Results.dirname;
dirrootpath = ip.Results.dirrootpath;
% Clearing memory
clear('ip');

% Constructing full directory path
dirpath = [dirrootpath filesep dirname];
dirpath = GetFullPath(dirpath);

% Creating directory
if ~exist(dirpath,'dir')
   mkdir(dirpath);
end;
end

%% Input parser construction function

function ip = buildIP()
% The purpose of this function is to construct an input parser object for the
% mkdirwchk function.
%
% INPUTS/OUTPUTS
% ip - Input parser object for the mkdirwchk function.

ip = inputParser;
ip.addRequired('dirname',@islegalfpath);
ip.addOptional('dirrootpath',pwd,@islegalfpath);
end