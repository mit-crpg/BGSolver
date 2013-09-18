%% vectorizeBVfunction
% Bond variable .m function file vectorization function.
%
% vectorizeBVfunction(fpath,scalarIDs)
%
% This function takes a BV function .m file (created by matlabFunction) that has
% at least one scalar constant/only time-dependent BV element, and makes sure
% the function is properly vectorized. The function must have at least one
% vectorized input (state or numeric variable).
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% BV .m function file vectorization function

function vectorizeBVfunction(fpath,scalarIDs)
% The purpose of this function is to vectorize a BV function .m file based on
% known information about the BV function as a whole and about the individual
% scalar (either constant or only time-dependent) elements of the BV function.
%
% INPUTS/OUTPUTS
% fpath - BV function .m file path.
% scalarIDs - Scalar element IDs vector within the BV function. These may be
%             both constant and only time-dependent. This vector can be empty,
%             in which case the BV function is not modified.

% Checking if vectorization is required
if isempty(scalarIDs)
   return;
end;

% Proceeding, assuming vectorization is required
% Opening the file
fid = fopen(fpath,'r+');
cleanupObj = onCleanup(@()cleanup(fid));

% Identifying the variable name to calculate nVect from
declStr = fgetl(fid);
varNameCell = regexp(declStr,'[\(,]\w+\)','match');
nVectVarName = varNameCell{1}(2:(end-1));

% Shifting to intended modification commentary line
flineseek(fid,7);

% Inserting the modification commentary line
modCommLine = sprintf(['%%\n%%    This function was modified by %s.\n' ...
   '%%    %s\n'],mfilename,datestr(now));
insertstrinfile(fid,modCommLine);

% Inserting the nVect calculation line
nVectLine = sprintf('\nnVect = size(%s,2); %% Inserted by %s\n',nVectVarName,...
   mfilename);
insertstrinfile(fid,nVectLine);

% Reading the output line
nlines = countflines(fid);
outpos = flineseek(fid,nlines-1);
outputline = fgetl(fid);

% Commenting out the output line
insertstrinfile(fid,sprintf('\n%% Output below replaced by %s\n%% ',...
   mfilename),outpos);

% Replacing the scalars in the output line with vectorized quantities
vectoutputline = vectorizescalars(outputline,scalarIDs);
fseek(fid,-1,'eof');
insertstrinfile(fid,vectoutputline);
end

%% Cleanup function

function cleanup(fid)
% The purpose of this function is to run the cleanup routine when
% vectorizeBVfunction exits, either because of an error, or regularly. The
% routine does the following:
% - Closes the file opened by vectorizeBVfunction
%
% INPUTS/OUTPUTS
% fid - File ID of the file opened by vectorizeBVfunction.

fclose(fid);
end

%% Scalar vectorization in the output line function

function vectoutputline = vectorizescalars(outputline,scalarIDs)
% The purpose of this function is to replace all scalar output elements with
% vectorized quantities. The scalar elements' in-output IDs must be provided.
%
% INPUTS/OUTPUTS
% outputline - Output line string to replace scalar output elements in.
% scalarIDs - Scalar element IDs vector within the BV function.
% vectoutputline - Vectorized output line.

% Extracting all elements' strings
outStrs = regexp(outputline,'(?<=(\[|;))[^\[\];]+(?=(;|\]))','match');

% Extracting output variable name
outvar = strtok(outputline);

% Vectorizing scalar elements
for oid = 1:length(scalarIDs)
   scalarID = scalarIDs(oid);
   outStrs{scalarID} = ['(' outStrs{scalarID} ')*ones(1,nVect)'];
end;

% Constructing new output line
vectoutputline = [outvar ' = [' cell2str(outStrs,';') '];'];
end