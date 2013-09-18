%% rmdirwchk
% Directory removal function with existence check.
%
% [status, message, messageid] = rmdirwchk(dirname,[dirrootpath])
%
% This function checks if a directory with a given name, and optionally a given
% root path, (absolute or relative), exists, and if it does, it removes it and
% its contents, if any. It returns the status, message and error message ID (if
% any) of the operation, identical to MATLAB built-in rmdir function.
%
% Package:    BGSolver v1.03
% Subpackage: Utilities
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Master Directory Removal with Existence Check Function code

function [status, message, messageid] = rmdirwchk(varargin)
% The purpose of this function is to check for a directory's existence, and if
% it does, remove it.
%
% INPUTS/OUTPUTS
% dirname - Directory name.
% dirrootpath - Directory root path. OPTIONAL. If not given, the current
%               directory will be used as root path.
% status - Logical scalar indicating the outcome of the operation, identical to
%          MATLAB built-in rmdir function status output.
% message - String containing the warning or error message text if the operation
%           is unsuccessful, identical to MATLAB built-in rmdir function status
%           output.
% messageid - String containing the warning or error message ID, if the
%             operation is unsuccessful, identical to MATLAB built-in rmdir
%             function status output.

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

% Removing directory
[status,message,messageid] = rmdir(dirpath,'s');
end

%% Input parser construction function

function ip = buildIP()
% The purpose of this function is to construct an input parser object for the
% rmdirwchk function.
%
% INPUTS/OUTPUTS
% ip - Input parser object for the rmdirwchk function.

ip = inputParser;
ip.addRequired('dirname',@islegalfpath);
ip.addOptional('dirrootpath',pwd,@islegalfpath);
end