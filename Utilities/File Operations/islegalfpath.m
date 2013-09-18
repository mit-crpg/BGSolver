%% islegalfpath
% File path legality checking function.
%
% iflegal = islegalfpath(fname,['ext',ext])
%  or
% iflegal = islegalfpath({fname1,fname2...},['ext',ext])
%
% This function checks one or more input strings for whether or not they can
% serve as legal full or local file paths. Each input is checked separately, and
% a boolean array is returned. A non-string entry returns a false. Optionally,
% an additional requirement of a specific extension can be added.
%
% Package:    BGSolver v1.03
% Subpackage: Utilities
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Single file path legality checking algorithm
%
% If PC:
%  Check if ':' is present at all. If not, conclude that the folder path
%  is partial, and skip a few steps. If yes:
%   Check if the path is 3 or less characters. If yes, instant fail.
%   If no, proceed.
%   Check if ':' is anywhere in the path other than the second symbol.
%   If yes, instant fail. If no, proceed.
%   Check if the first symbol is a letter and the third symbol is either
%   a '\' or a '/'. If yes, proceed. If no, instant fail.
%   If made it this far, conclude that the folder path is full, and
%   proceed to splitting into 3 parts.
%  Split the file path into the 3 parts.
%  Make the following checks:
%   Folder string can be empty if partial folder path is given.
%   If folder string is partial, it cannot contain: ':*?"<>|'.
%   If folder string is full, its (4,end) characters cannot contain
%   ':*?"<>|'. If the folder string is 3 characters long, proceed.
%   File name string cannot be empty.
%   File name cannot contain: '/\:*?"<>|'.
%   Extension can be empty. If it's not, extension legality is checked
%   by isext function.
%    If the extension is legal, it must also match the specified
%    extension (if any).
%  If all these checks are passed, the file path is legal. Otherwise, it
%  is not.
% If Unix or Mac:
%  Check if '~' is present at all. If not, proceed to splitting the file
%  path into the 3 parts. If yes:
%   Check if the path is 2 or less characters. If yes, instant fail. If
%   no, proceed.
%   Check if '~' is anywhere in the path other than the first character.
%   If yes, instant fail. If not, proceed.
%   If made it this far, conclude that a folder string is present, and
%   proceed to splitting the file into the 3 parts.
%  Split the file path into the 3 parts.
%  Make the following checks:
%   Folder string can be empty, unless it's supposed to be present due
%   to the presence of '~'.
%   First character of folder string may be '~', others cannot be.
%   File name string cannot be empty.
%   File name cannot contain: '~/'
%   Extension can be empty. If it's not, extension legality is checked
%   by isext function.
%    If the extension is legal, it must also match the specified
%    extension (if any).
%  If all these checks are passed, the file path is legal. Otherwise, it
%  is not.

%% Main file path legality checking function

function iflegal = islegalfpath(varargin)
% The purpose of this function is to check one or more input strings for whether
% or not they can serve as legal full or local file paths. Each input is checked
% separately, and a boolean array is returned. A non-string entry returns a
% false. Optionally, an additional requirement of a specific extension can be
% added.
%
% INPUTS/OUTPUTS
% fname - Potential file name string.
% fnames - Potential file name strings cell array.
% ext - Dotted file extension required for a legal file path. OPTIONAL.
% iflegal - Boolean array of whether or not the corresponding potential
%           file names can serve as legal full or local file paths.

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
if ~isempty(ip.Results.ext)
   ifext = true;
   ext = ip.Results.ext;
else
   ifext = false;
end;
% Clearing memory
clear('ip');

% Checking path legality
if ifcell
   if ifext
      strs = repmat({'ext'},size(fpaths));
      exts = repmat({ext},size(fpaths));
      iflegal = cellfun(@islegalfpath,fpaths,strs,exts);
   else
      iflegal = cellfun(@islegalfpath,fpaths);
   end;
elseif ~isempty(fpath)
   iflegal = true; % Assuming legal path
   
   if ispc
      % Checking if the file path is full or partial
      if any(fpath == ':')
         if length(fpath) <= 3
            iflegal = false;
         else
            if any([fpath(1),fpath(3:end)] == ':')
               iflegal = false;
            else
               if isletter(fpath(1)) && any(fpath(3) == '\/')
                  pathtype = 'full';
               else
                  iflegal = false;
               end;
            end;
         end;
      else
         pathtype = 'partial';
      end;
      
      % Checking folder path legality
      if iflegal
         % Splitting the file path
         [folderpath, fname, fext] = fileparts(fpath);
         
         illegalFolderChars = ':*?"<>|';
         if strcmp(pathtype,'full')
            if any(strcharchk(folderpath(4:end),illegalFolderChars))
               iflegal = false;
            end;
         elseif strcmp(pathtype,'partial')
            if any(strcharchk(folderpath,illegalFolderChars))
               iflegal = false;
            end;
         end;
      end;
      
      % Checking file name legality
      if iflegal
         if isempty(fname)
            iflegal = false;
         else
            illegalFNameChars = '/\:*?"<>|';
            if any(strcharchk(fname,illegalFNameChars))
               iflegal = false;
            end;
         end;
      end;
      
      % Checking extension legality
      if iflegal
         if ifext
            if ~strcmpi(fext,ext)
               iflegal = false;
            end;
         else
            if ~isempty(fext) && ~isext(fext)
               iflegal = false;
            end;
         end;
      end;
   elseif isunix || ismac
      if any(fpath == '~')
         if length(fpath) <= 2 || any(strcharchk(fpath(2:end),'~'))
            iflegal = false;
         else
            iftilde = true;
         end;
      else
         iftilde = false;
      end;
      
      % Checking folder path legality
      if iflegal
         % Splitting the file path
         [folderpath, fname, fext] = fileparts(fpath);
         
         if iftilde
            if isempty(folderpath)
               iflegal = false;
            end;
         end;
      end;
      
      % Checking file name legality
      if iflegal
         if isempty(fname)
            iflegal = false;
         else
            illegalFNameChars = '~/';
            if any(strcharchk(fname,illegalFNameChars))
               iflegal = false;
            end;
         end;
      end;
      
      % Checking extension legality
      if iflegal
         if ifext
            if ~strcmp(fext,ext)
               iflegal = false;
            end;
         else
            if ~isempty(fext) && ~isext(fext)
               iflegal = false;
            end;
         end;
      end;
   else
      error('Unknown platform!!!');
   end;
else
   iflegal = false;
end;
end

%% Input parser construction function

function ip = buildIP()
% The purpose of this function is to construct an input parser object for
% the islegalfpath function.
%
% INPUTS/OUTPUTS
% ip - Input parser object for the islegalfpath function.

ip = inputParser;
ip.addRequired('arg',@(arg)(iscellstr(arg) || ischar(arg)));
ip.addParamValue('ext',[],@isext);
end

%% File extension legality checking function

function iflegalext = isext(ext)
% The purpose of this function is to check if the supplied input can serve as a
% legal dotted extension.
%
% INPUTS/OUTPUTS
% ext - Input to check.
% iflegalext - Whether or not the input can serve as a legal dotted
%              extension.

if ~ischar(ext) || length(ext) < 2
   iflegalext = false;
else
   % An extension starts with a dot and can have anything other than
   % whitespace or punctuation in it.
   if (ext(1) == '.') && ~any(isstrprop(ext(2:end),'punct')) && ...
         ~any(isstrprop(ext(2:end),'wspace'))
      iflegalext = true;
   else
      iflegalext = false;
   end;
end;
end

%% String char-by-char checking function

function charstats = strcharchk(str,chars)
% The purpose of this function is to check if a given string contains any of a
% list of characters. A logical array is returned, with true corresponding to
% characters present in the string, and false correspondign to characters absent
% in the string.
%
% INPUTS/OUTPUTS
% str - String to check character-by-character.
% chars - Character array to check the string against.
% charstats - Logical array, true corresponds to character present in the
%             string, false corresponds to character not present.

f = @(ch) any(str == ch);
charstats = arrayfun(f,chars);
end