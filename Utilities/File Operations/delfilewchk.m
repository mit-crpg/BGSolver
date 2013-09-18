%% delfilewchk
% File deletion function with existence check.
%
% [ifexist,ifdel] = delfilewchk(fname,[frootpath])
%
% This function checks if a file with a given name, and optionally a given root
% path, (absolute or relative), exists, and if it does, it deletes it. It
% returns the status of whether or not the file existed, and whether or not the
% deletion was successful).
%
% Package:    BGSolver v1.03
% Subpackage: Utilities
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Master File Deletion with Existence Check Function code

function [ifexist,ifdel] = delfilewchk(varargin)
% The purpose of this function is to check for a file's existence, and if it
% does, delete it.
%
% INPUTS/OUTPUTS
% fname - File name, absolute or relative.
% frootpath - File root path, absolute or relative. OPTIONAL. If not given, the
%             current directory will be used as root path.
% ifexist - Whether or not the file existed.
% ifdel - Whether or not the deletion was successful. If the file did not exist,
%         deletion is assumed not successful.

% Input parser construction
ip = buildIP();
% Input parsing
ip.parse(varargin{:});
% Retrieving inputs
fname = ip.Results.fname;
frootpath = ip.Results.frootpath;
% Clearing memory
clear('ip');

% Constructing full file path
fullfname = [frootpath filesep fname];
fullfname = GetFullPath(fullfname);

% Checking file existence
ifexist = exist(fullfname,'file') == 2;

if ifexist
   try
      delete(fullfname);
      ifdel = true;
   catch err %#ok<NASGU>
      ifdel = false;
   end;
else
   ifdel = false;
end;
end

%% Input parser construction function

function ip = buildIP()
% The purpose of this function is to construct an input parser object for the
% delfilewchk function.
%
% INPUTS/OUTPUTS
% ip - Input parser object for the delfilewchk function.

ip = inputParser;
ip.addRequired('fname',@islegalfpath);
ip.addOptional('frootpath',pwd,@islegalfpath);
end