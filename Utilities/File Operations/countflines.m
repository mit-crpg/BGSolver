%% countflines
% File line number counting function (cross-platform).
%
% nlines = countflines(fid)
%  OR
% nlines = countflines(fname)
%
% This function counts the number of lines in either an open (with MATLAB
% built-in fopen function) or a named text file.
% If the last line is empty, it is still counted by the Windows procedure. It
% has not yet been tested if the Linux/Unix/Mac procedure counts such line.
%
% Package:    BGSolver v1.03
% Subpackage: Utilities
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Master File Line Number Counting Function code

function nlines = countflines(finput)
% The purpose of this function is to count the number of lines in a text file.
%
% INPUTS/OUTPUTS
% finput - File input variable. Can be either an open fid (numeric) or a file
%          name.
% nlines - Number of lines in the file.

% Opening file, if necessary
inpclass = class(finput);
if strcmp(inpclass,'double')
   fid = finput;
   fname = fopen(fid);
   oldpos = ftell(fid);
   frewind(fid);
elseif strcmp(inpclass,'char');
   fid = fopen(finput,'r');
   cleanupObj = onCleanup(@()cleanup(fid));
   fname = finput;
   oldpos = 0;
end;

% Checking system architecture
if ispc % Windows
   % Using MATLAB built-in functions
   nlines = 0;
   while (fgets(fid) ~= -1),
      nlines = nlines + 1;
   end;
   % Correcting for last empty line, if present
   fseek(fid,-1,'eof');
   c = fscanf(fid,'%c',1);
   n = sprintf('\n');
   if c == n
      nlines = nlines + 1;
   end;
elseif ismac || isunix % Linux/Unix/Mac
   % Using system command
   [~,nlinesStr] = system(sprintf('wc -l "%s"',fname));
   nlines = str2double(nlinesStr);
end;

% Going back to original position
fseek(fid,oldpos,'bof');
end

%% Cleanup function

function cleanup(fid)
% The purpose of this function is to run the cleanup routine when countflines
% exits, either because of an error, or regularly. The routine does the
% following:
% - Closes the file, if it was opened by countflines
%
% INPUTS/OUTPUTS
% fid - File ID of the file opened by countflines.

fclose(fid);
end