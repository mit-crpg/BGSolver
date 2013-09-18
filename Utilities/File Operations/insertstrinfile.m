%% insertstrinfile
% In-file string insertion function.
%
% newpos = insertstrinfile(fid,str,[pos])
%
% This function inserts a string into a file. The file must already be open for
% reading and writing (using MATLAB built-in fopen with r+ permission). The
% string is inserted either into the current pointer position, or into a the
% specified position, if supplied. It returns the in-file position at the ened
% of the inserted string, and leaves the pointer in that position as well. The
% function does not make any corrections (like piecewise copying) for very large
% files, so if the new file does not fit into RAM, significant slowdown may
% result.
%
% Package:    BGSolver v1.03
% Subpackage: Utilities
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Master In-file String Insertion Function code

function newpos = insertstrinfile(fid,str,pos)
% The purpose of this function is to insert a string into an open file.
%
% INPUTS/OUTPUTS
% fid - File ID to insert the string into. The file must have been open with
%       both reading and writing permission.
% str - String to insert. No special characters are added to the string.
% pos - Position to insert the string into. OPTIONAL. Default value: current
%       in-file position.
% newpos - Position at the end of the inserted string.

% Shifting the pointer, if necessary
if nargin == 3
   fseek(fid,pos,'bof');
else
   pos = ftell(fid);
end;

% Reading the rest of the file
fend = fread(fid,inf,'uchar');
fseek(fid,pos,'bof');

% Inserting the string
fprintf(fid,'%s',str);
newpos = ftell(fid);

% Appending the rest of the file
fwrite(fid,fend,'uchar');

% Seeking to the end of the inserted string
fseek(fid,newpos,'bof');
end