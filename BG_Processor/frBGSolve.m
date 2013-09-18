%% frBGSolve
% Bond graph processing restart function.
%
% [T,X,[E,F]] = frBGSolve(mfname,['rootpath',rootpath,'rfname',rfname,
%     'rrunopts',rrunopts,'X0',X0,'NonNegative',NonNegative,'AbsTol',AbsTol])
%
% This function restarts the processing of a bond graph system from the
% beginning of Bond Graph Process Step 7. It takes a Sorted .mat Results (SMR)
% file name, and reads the partially processed system (up to and including Step
% 6) from it.
% To correctly link to relative file and directory paths, the function needs to
% know the problem root path. Optionally, the function accepts this root path as
% an argument. If it is not given as an argument, it can infer it from the
% as-saved SMR file path (if the as-saved SMR file path was relative, or if the
% SMR file path has not changed). The problem name itself can also be modified,
% which will affect the RRM, RTRR and Log file names, if they are default.
% Optionally, the function also accepts a restart run options data structure,
% which can modify some of the options and/or file paths used in Step 7. The
% full list of adjustable Step 7 options is listed in the Adjustable options
% cell below.
% Additionally, the function can accept a new initial values vector, a new
% non-negative indices (either as numeric, or as logical indices) vector and a
% new absolute tolerances vector as parameters. These vectors must be in order
% of increasing storage element ID. Generally, in BGSolver v1.03, initial values
% are considered a part of the bond graph system, and not a run option, and so
% are not intended to be adjustable. For this reason, this option is currently
% considered in "beta" mode.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Adjustable options
%
% All adjustable options are listed below. These options are named according to
% their respective fields in the Runtime and Recording Options data structure
% (runopts), which is documented in detail in fBGSolve Variable descriptions
% cell.
% NOTE:
% 1. Almost no checking is done on these options. If they are not internally
%    consistent, a difficult-to-debug error will likely result.
% 2. While some spec-type options are available (ones that specify how an option
%    is specified - through default value, through an asked-at-runtime value, or
%    through an in-file specified value), only the default and "in-file" (i.e.,
%    in-rrunopts) options are allowed. Spec-type options that only provide the
%    options of in-file and asked-at runtime are not adjustable.
% 3. Unlike the options supplied to fBGSolve, here, only a very limited
%    completion function is executed.
%
% SECTION 2 - Runtime and Recording Options.
% .verbose - How verbose the processor should be during processing.
% .log - How to log all the outputs specified by .verbose.
% .logfnspec - How the log file name is specified.
% .logfn - Log file name.
% .rtrr - Whether or not to record the intergration results as runtime, as they
%         are being computed.
% .rtrrfq - Runtime results recording frequency.
% .rtrrfnspec - How the runtime results file name is specified.
% .rtrrfn - Runtime results file name.
% .rrm - Whether or not to record the final results in a .mat file.
% .rrmfnspec - How the final results .mat file name is specified.
% .rrmfn - Final results .mat file name.
% SECTION 3 - Time Integrator Options.
% .tit - Time integrator type.
% .ii - Integration interval.
% .tst - Time step type, either adaptive or fixed.
% .ftsspect - What the fixed time spec is specified in terms of.
% .fts - Fixed time step.
% .phxsplit - Whether, and how, the physics are split during integration.
% .addtio - Number of additional time integrator options.
% .addtiolist - [addtio x 2] cell array of additional time integrator options
%               pairs (parameter name, parameter value).
% .tidisp - Time integrator display level.
% SECTION 4 - Solver and Jumpstart Options.
% .stolx - Solver tolerance on state vector.
% .stolres - Solver tolerance on residual.
% .sdisp - Solver display level.
% .srscal - Solver residual scaling.
% .sxscal - Solver state vector scaling.
% .sxcent - Solver state vector centering.
% .spc - Solver preconditioner.
% .sigpdct - Solver initial guess predictor.
% .sjacob - Solver Jacobian construction.
% .addso - Number of additional solver options.
% .addsolist - [addso x 2] cell array of additional solver options pairs
%              (parameter name, parameter value).
% .jstype - Jumpstart type.
% .jsmsnum - Jumpstart ministep number per 1 full time step.
% .jsphxsplit - Whether, and how, the physics are split during jumpstart.
% .jstolx - Solver tolerance on state vector to use during jumpstart.
% .jstolres - Solver tolerance on residual to use during jumpstart.
% .jstit - Time integrator type to use during jumpstart.
% .jsdisp - Jumpstart solver display level.
% .jsrscal - Jumpstart residual scaling.
% .jsxscal - Jumpstart state vector scaling.
% .jsxcent - Jumpstart state vector centering.
% .jspc - Jumpstart preconditioner.
% .jsjacob - Jumpstart Jacobian construction.
% .addjso - Number of additional jumpstart solver options.
% .addjsolist - [addjso x 2] cell array of additional jumpstart solver options
%               pairs (parameter name, parameter value).
% SECTION 5 - Evaluator Options.
% - No options in this section are adjustable.
% SECTION 6 - Post-Processor Options.
% .epp - Number of effort variables to evaluate during post-processing.
% .eppIDs - List of BIDs of effort varaibles to evaluate during post-processing.
% .fpp - Number of flow variables to evaluate during post-processing.
% .fppIDs - List of BIDs of flow variables to evaluate during post-processing.
% .ppdisp - Post-processing display level.
% SECTION 7 - Linked Files and Directories.
% .lflist - List of linked file names.
% .ldlist - List of linked directory names.

%% Main BGS processing restart function

function [T,X,E,F] = frBGSolve(varargin)
% The purpose of this function is to restart the processing of a bond graph
% system from the beginning of Bond Graph Process Step 7.
%
% INPUTS/OUTPUTS
% mfname - Sorted .mat Results file name. This path must be absolute, or
%          relative off of the current directory, NOT off the problem root path.
%          REQUIRED.
% rootpath - Problem Root Path. OPTIONAL.
% rfname - Restart (problem) File Name. OPTIONAL.
% rrunopts - Restart Runtime and Recording Options data structure. OPTIONAL.
% X0 - Initial value vector. OPTIONAL.
% NonNegative - Nonnegative state indices vector, either logical indexing or
%               numeric indices. OPTIONAL.
% AbsTol - Absolute tolerances vector. OPTIONAL.
% T - Time point vector.
% X - State solution array.
% E - Effort solution array.
% F - Flow solution array.
% See the adjustable options cell for more information on rrunopts.

% Input parser construction
ip = buildIP();
% Input parsing
ip.parse(varargin{:});
% Retrieving inputs
mfname = ip.Results.mfname;
rootpath = ip.Results.rootpath;
rfname = ip.Results.rfname;
rrunopts = ip.Results.rrunopts;
X0 = ip.Results.X0;
NonNegative = ip.Results.NonNegative;
AbsTol = ip.Results.AbsTol;
% Clearing memory
clear('ip');

% Loading SMR file
load(mfname);

% Extracting problem name and as-saved root path
[srootpath,sfname,~] = fileparts(fname);
% Updating as-loaded problem name, if empty
if isempty(rfname)
   rfname = sfname;
end;

% Identifying root directory
[rootpath,rpmethod] = buildrootpath(rootpath,mfname,srootpath,...
   runopts); %#ok<NODEF>
% Checking root directory inference
assert(~isempty(rootpath),'Unable to infer root path!!!');

% Switching directory to the root
oldpath = cd(rootpath);

% Constructing temporary cleanup routine 1
cleanupObj1 = onCleanup(@()cleanup(1,'oldpath',oldpath));
% Clearing memory
clear('oldpath');

% Updating file and directory paths to account for the new root path
runopts = updatefdpaths(runopts,rootpath,srootpath);

% Storing old run options
oldrunopts = runopts;

% Overwriting run options (if necessary)
if ~isempty(rrunopts)
   runopts = dsfwovrwrite(runopts,rrunopts);
end;
% Clearing memory
clear('rrunopts');

% Expanding runopts
[runopts,ifaddtio,ifaddso,ifaddjso] = expandrunopts(runopts,rfname,oldrunopts);

% Creating required subdirectories
if ~exist('Results','dir')
   mkdir('Results');
end;
if ~exist('TmpFunctions','dir')
   mkdir('TmpFunctions');
end;
% Linking required subdirectories
addpath('Results','TmpFunctions','-begin');

% Creating/opening required files
fids = createfiles(runopts);

log2spec('BGSolver v1.03 Restart by Eugeny Sosnovsky.',1,fids,runopts);
log2spec('Restarting bond graph system processing...',1,fids,runopts);
log2spec('Inputs have been accepted.',1,fids,runopts);
switch rpmethod
   case 0
      log2spec(['Sorted .mat Results file was not moved, so the as-saved ' ...
         'root path was used.'],1,fids,runopts);
   
   case 1
      log2spec('Root path was explicitly given as argument.',1,fids,runopts);
   
   case 2
      log2spec(['Root path inferred to be a same degree parent to the ' ...
         'Sorted .mat Results file as the as-saved root path was.'],...
         1,fids,runopts);
end;

log2spec('Runtime and Recording Options data structure has been completed.',...
   2,fids,runopts);
log2spec('Results and TmpFunctions subdirectories have been added.',...
   2,fids,runopts);
log2spec('Required files have been created/opened (if any).',2,fids,runopts);
log2spec(['Cleanup routine 1 (switch to original directory) ' ...
   'has been constructed.'],2,fids,runopts);

% Constructing temporary cleanup routine 2
cleanupObj2 = onCleanup(@()cleanup(2,'fids',fids));

log2spec(['Cleanup routine 2 (closing of all open files) has been ' ...
   'constructed.'],2,fids,runopts);

log2spec('---',1,fids,runopts);
log2spec('Beginning to link required files and folders...',1,fids,runopts);

% Linking required files and folders
for lfid = 1:runopts.lfnum
   [lfpath,lfname,lfext] = fileparts(runopts.lflist{lfid});
   if strcmpi(lfext,'.m')
      addpath(lfpath,'-begin');
      eval(lfname);
      rmpath(lfpath);
   elseif strcmpi(lfext,'.mat')
      load(runopts.lflist{lfid});
   else
      error('Linked file list entry %d has invalid extension!!!',lfid);
   end;
end;
for ldid = 1:runopts.ldnum
   addpath(runopts.ldlist{ldid},'-begin');
end;
load('ti_info.mat');

log2spec('Required files and folders have been successfully linked.',1,fids,...
   runopts);

% Constructing temporary cleanup routine 3
cleanupObj3 = onCleanup(@()cleanup(3,'runopts',runopts,'rootpath',rootpath));
% Clearing memory
clear('rootpath');

log2spec(['Cleanup routine 3 (removing linked directories from path) has ' ...
   'been constructed.'],2,fids,runopts);

%
% STUB
% Pre-restart error checking procedures can be inserted here.
% /STUB
%

log2spec('---',1,fids,runopts);
log2spec('BOND GRAPH PROCESS STEP 6 (Adjustment substeps only).',...
   1,fids,runopts);

% BOND GRAPH PROCESS STEP 6
% State derivative vector is already formulated, and cannot be adjusted.
% Adjusting and reordering initial values vector
log2spec('Adjusting and reordering initial values vector...',1,fids,runopts);
if ~isempty(X0)
   X0 = X0(srtsys.newxinds); %#ok<NODEF>
   log2spec('Initial values vector adjusted and reordered.',1,fids,runopts);
else
   X0 = srtsys.x0s; %#ok<NODEF>
   log2spec(['Initial values vector did not require adjustment and ' ...
      'reordering.'],1,fids,runopts);
end;
% Adjusting and reordering nonnegative indices vector
if ~isempty(NonNegative)
   if islogical(NonNegative)
      NonNegative = NonNegative(srtsys.newxinds);
   elseif isnumeric(NonNegative)
      NonNegativeTmp = false(1,srtsys.xnum);
      NonNegativeTmp(NonNegative) = true;
      NonNegative = NonNegativeTmp(srtsys.newxinds);
   end;
   NonNegative = uint16(find(NonNegative));
end;
% Adjusting and reordering absolute tolerances vector
if ~isempty(AbsTol)
   AbsTol = AbsTol(srtsys.newxinds);
end;
% Adjusting tioptions
log2spec('Adjusting time integrator options...',1,fids,runopts);
[tioptions,iftiadj] = adjusttioptions(tioptions,fids,runopts,oldrunopts,...
   ifaddtio,ifaddso,ifaddjso,NonNegative,AbsTol); %#ok<NODEF>
if iftiadj
   log2spec('Time integrator options adjusted.',1,fids,runopts);
else
   log2spec('Time integrator options did not require adjustment.',1,fids,...
      runopts);
end;
% Reformulating post-processing vector
log2spec('Reformulating post-processing vector...',1,fids,runopts);
[ppsys,ifppadj,srtsys] = reformppsys(ppsys,srtsys,vars,refs,runopts,...
   oldrunopts); %#ok<NODEF>
if ifppadj
   log2spec('Post-processing vector reformulated.',1,fids,runopts);
else
   log2spec('Post-processing did not require reformulation.',1,fids,runopts);
end;
% Extracting reverse ordering X vector
revxinds = srtsys.revxinds;
log2spec('State vector reverse reordering vector extracted.',2,fids,runopts);
% Clearing memory
clear('fname','srtsys','vars','oldrunopts','NonNegative','NonNegativeTmp',...
   'AbsTol');

log2spec('---',1,fids,runopts);
log2spec('BOND GRAPH PROCESS STEP 7.',1,fids,runopts);
logtiinfo(fids,runopts,tioptions,ti_info);
log2spec('Integrating state derivative vector...',1,fids,runopts);

% BOND GRAPH PROCESS STEP 7
% Integrating the state derivative vector
[T,X] = runti(xdot,runopts,X0,tioptions);

log2spec('State derivative vector integrated.',1,fids,runopts);
log2spec('---',1,fids,runopts);

% POST-PROCESSING STEP
if ppsys.ifpp
   log2spec('POST-PROCESSING STEP.',1,fids,runopts);
   log2spec('Beginning post-processing array evaluation...',1,fids,runopts);
   [E,F] = evalpp(T,X,ppsys);
   log2spec('Post-processing complete.',1,fids,runopts);
else
   E = [];
   F = [];
   log2spec('No post-processing specified.',1,fids,runopts);
end;

log2spec('---',1,fids,runopts);

% Reordering the X vector
X = X(:,revxinds);
log2spec('State vector reordered back into ascending EID order.',...
   2,fids,runopts);
log2spec('---',2,fids,runopts);

% Recording results
if runopts.rrm
   save(runopts.rrmfn,'T','X','E','F');
   
   log2spec('Results recorded in .mat file:',1,fids,runopts);
   log2spec(['--> ' runopts.rrmfn],1,fids,runopts);
   log2spec('---',1,fids,runopts);
end;

log2spec('Restarted bond graph system processing complete.',1,fids,runopts);
log2spec('It''s all over. Noone will ever love you. I am quitting...',...
   1,fids,runopts);

% Cleanup routines 1-3 executed here, in arbitrary order, but not
% simultaneously.
end

%% Input parser construction function

function ip = buildIP()
% The purpose of this function is to construct an input parser object for the
% frBGSolve function.
%
% INPUTS/OUTPUTS
% ip - Input parser object for the frBGSolve function.

mfchkfun = @(mfname) islegalfpath(mfname,'ext','.mat');
rpchkfun = @(rootpath) isempty(rootpath) || ~any(ismember('*?"<>|',rootpath));
rfchkfun = @(rfname) isempty(rfname) || ~any(ismember('*?"<>|/\',rfname));
rroptschkfun = @(rrunopts) isempty(rrunopts) || isstruct(rrunopts);
X0chkfun = @(X0) isempty(X0) || (isa(X0,'double') && isvector(X0));
NNchkfun = @(NonNegative) isempty(NonNegative) || (isvector(NonNegative) && ...
   (isnumeric(NonNegative) || islogical(NonNegative)));
ATchkfun = @(AbsTol) isempty(AbsTol) || (isvector(NonNegative) && ...
   isa(AbsTol,'double') && all(AbsTol > 0));

ip = inputParser;
ip.StructExpand = false;
ip.addRequired('mfname',mfchkfun);
ip.addParamValue('rootpath','',rpchkfun);
ip.addParamValue('rfname','',rfchkfun);
ip.addParamValue('rrunopts',[],rroptschkfun);
ip.addParamValue('X0',[],X0chkfun);
ip.addParamValue('NonNegative',[],NNchkfun);
ip.addParamValue('AbsTol',[],ATchkfun);
end

%% Cleanup function

function cleanup(varargin)
% The purpose of this function is to run the cleanup routines when frBGSolve
% exits, either because of an error, or regularly. Which cleanup routine is run
% depends on the input settings.
%
% NOTE
% This implementation assumes that the cleanup routines are executed in an
% arbitrary order, but sequentially, never simultaneously. I asked a qustion
% about this on the MATLAB Answer forum:
% http://www.mathworks.com/matlabcentral/answers/...
%  2245-oncleanup-order-of-execution
% According to that link, this is the case: cleanup routines are executed in
% arbitrary order, but not simultaneously.
%
% INPUTS/OUTPUTS
% opt - Input setting which determines which cleanup routine to run. Possible
%       values:
%  1 - Cleanup routine 1: switch to original directory.
%  2 - Cleanup routine 2: closing of all open files.
%  3 - Cleanup routine 3: removing linked directories from path.
% oldpath - Original directory that frBGSolve is called from. REQUIRED.
% fids - File IDs data structure. OPTIONAL.
% runopts - Runtime and Recording options data structure. OPTIONAL.
% rootpath - Absolute root folder path. OPTIONAL.

% Constructing input parser
ip = inputParser;
ip.StructExpand = false;
fopt = @(opt) any(opt == [1,2,3,4]);
foldpath = @(oldpath) exist(oldpath,'dir');
ip.addRequired('opt',fopt);
ip.addParamValue('oldpath',pwd,foldpath);
ip.addParamValue('fids',struct(),@isstruct);
ip.addParamValue('runopts',struct(),@isstruct);
ip.addParamValue('rootpath',pwd,@ischar);
% Parsing input
ip.parse(varargin{:});
% Retrieving inputs
opt = ip.Results.opt;
oldpath = ip.Results.oldpath;
fids = ip.Results.fids;
runopts = ip.Results.runopts;
rootpath = ip.Results.rootpath;
% Clearing memory
clear('ip','fopt','foldpath');

if opt == 1
   % Cleanup routine 1: switching to original directory.
   cd(oldpath);
elseif opt == 2
   % Cleanup routine 2: closing all open files.
   fidfields = fieldnames(fids);
   f = @(fldname) fclose(fids.(fldname));
   cellfun(f,fidfields);
elseif opt == 3
   % Cleanup routine 3: removing linked directories from path.
   oldpath = cd(rootpath);
   for ldid = 1:runopts.ldnum
      rmpath(runopts.ldlist{ldid});
   end;
   rmpath('TmpFunctions','Results');
   cd(oldpath);
end;
end

%% Root directory path construction function

function [rootpath,rpmethod] = buildrootpath(rootpath,mfname,srootpath,runopts)
% The purpose of this function is to construct the root directory path by either
% using the one passed in by the user, or inferring it from the loaded SMR file
% path and saved variables. If it fails to infer a root path, a blank array is
% returned. Additionally, the function returns a description of the method
% through which it inferred the root path.
%
% INPUTS/OUTPUTS
% rootpath - Problem Root Path.
% mfname - Sorted .mat Results file name.
% srootpath - As-saved root path.
% runopts - Runtime and Recording Options data structure from the SMR file.
% rpmethod - Method through which the function inferred the root path. It can be
%            one of the following values:
%  -1 - The function failed to infer root path.
%  0 - SMR file has not been moved, so the as-saved root path is used.
%  1 - A root path was explicitly given.
%  2 - As-saved root path was a parent directory to as-saved SMR file, so the
%      function infers that the root path is a same-degree parent directory to
%      the current SMR file.

if isempty(rootpath)
   % Constructing absolute SMR file paths
   smfname = buildrdfpath(runopts.smfn,srootpath); % As-saved
   lmfname = GetFullPath(mfname);
   
   % Comparing as-saved and as-loaded SMR file paths
   if strcmp(smfname,lmfname)
      rootpath = srootpath; % SMR file has not been moved
      rpmethod = 0;
   else
      % Checking if as-saved root directory was parent to as-saved SMR
      if strcmp(smfname(1:length(srootpath)),srootpath)
         % Constructing relative path from as-saved SMR path to as-saved root
         % directory
         relpath = buildrpath(smfname,srootpath);
         
         % Constructing rootpath
         rootpath = GetFullPath(fullfile(fileparts(lmfname),relpath));
         rootpath = trimfilesep(rootpath);
         rpmethod = 2;
      else
         % Unable to infer root path
         rootpath = '';
         rpmethod = -1;
      end;
   end;
else
   rootpath = GetFullPath(rootpath);
   rootpath = trimfilesep(rootpath);
   rpmethod = 1;
end;
end

%% File and directory path root update function

function runopts = updatefdpaths(runopts,rootpath,srootpath)
% The purpose of this function is to update file and directory paths in runopts
% to account for the new root directory path.
%
% INPUTS/OUTPUTS
% runopts - Runtime and Recording Options data structure.
% rootpath - As-loaded root directory path.
% srootpath - As-saved root directory path.

% File field names to update
fldnames = {'logfn','rrmfn','rtrrfn'};
% Updating file fields
for fldid = 1:length(fldnames)
   fld = fldnames{fldid};
   if isfield(runopts,fld) && isfullpath(runopts.(fld))
      runopts.(fld) = reprootpath(runopts.(fld),srootpath,rootpath);
   end;
end;

% Linked object lists to update
nfldnames = {'lfnum','ldnum'}; % Object numbers
listfldnames = {'lflist','ldlist'}; % Object list names
% Updating linked objects
for fldid = 1:length(nfldnames)
   nfld = nfldnames{fldid};
   if runopts.(nfld) > 0
      listfld = listfldnames{fldid};
      for objid = 1:runopts.(nfld)
         if isfullpath(runopts.(listfld){objid})
            runopts.(listfld){objid} = reprootpath(runopts.(listfld){objid},...
               srootpath,rootpath);
         end;
      end;
   end;
end;
end

%% Runopts expansion and completion function

function [runopts,ifaddtio,ifaddso,ifaddjso] = expandrunopts(runopts,fname,...
   oldrunopts)
% The purpose of this function is to replace the default specifications with the
% corresponding settings and file names, to delete the logging specifications if
% they are no longer relevant, and to compute the fixed time step size, if
% required. It also evaluates the additional time integrator options, if
% necessary.
%
% INPUTS/OUTPUTS
% runopts - Current Runtime and Recording Options data structure.
% oldrunopts - As-loaded Runtime and Recording Options data structure.
% fname - Problem file name.
% ifaddtio - Whether or not additional time integrator options were adjusted.
% ifaddso - Whether or not additional solver options were adjusted.
% ifaddjso - Whether or not additional jumpstart solver options were adjusted.

% Replacing default specs with corresponding settings and file names
if ~isfield(runopts,'verbose')
   runopts.verbose = 1;
end;
if runopts.verbose == 0
   runopts = rmfieldcheck(runopts,{'log','logfnspec','logfn'});
end;
if runopts.verbose ~= 0 && ~isfield(runopts,'log')
   runopts.log = 0;
end;
if runopts.verbose ~= 0 && runopts.log ~= 0 && ~isfield(runopts,'logfnspec')
   runopts.logfnspec = 0;
end;
if runopts.verbose ~= 0 && runopts.log ~= 0 && runopts.logfnspec == 0
   runopts.logfn = ['Results' filesep fname '.log'];
end;

if ~isfield(runopts,'rtrr')
   runopts.rtrr = 0;
end;
if runopts.rtrr == 1 && ~isfield(runopts,'rtrrfnspec')
   runopts.rtrrfnspec = 0;
end;
if runopts.rtrr == 1 && runopts.rtrrfnspec == 0
   runopts.rtrrfn = ['Results' filesep fname '_RTRes.txt'];
end;

if ~isfield(runopts,'rrm')
   runopts.rrm = 1;
end;
if runopts.rrm == 1 && ~isfield(runopts,'rrmfnspec')
   runopts.rrmfnspec = 0;
end;
if runopts.rrm == 1 && runopts.rrmfnspec == 0
   runopts.rrmfn = ['Results' filesep fname '_FinRes.mat'];
end;

if ~isfield(runopts,'addtio')
   if isfield(runopts,'addtiolist')
      runopts.addtio = uint8(size(runopts.addtiolist,1));
   else
      runopts.addtio = uint8(0);
   end;
end;

if ~isfield(runopts,'addso')
   if isfield(runopts,'addsolist')
      runopts.addso = uint8(size(runopts.addsolist,1));
   else
      runopts.addso = uint8(0);
   end;
end;

if ~isfield(runopts,'addjso')
   if isfield(runopts,'addjsolist')
      runopts.addjso = uint8(size(runopts.addjsolist,1));
   else
      runopts.addjso = uint8(0);
   end;
end;

% Evaluating the additional time integrator options
ifaddtio = false;
if isfield(runopts,'addtiolist')
   if isfield(oldrunopts,'addtiolist') && runopts.addtio == oldrunopts.addtio
      for addtioID = 1:runopts.addtio
         oldaddtio = oldrunopts.addtiolist{addtioID,2};
         addtio = runopts.addtiolist{addtioID,2};
         
         if ischar(addtio) && ((ischar(oldaddtio) && ...
               ~strcmp(oldaddtio,addtio)) || ~ischar(oldaddtio))
            runopts.addtiolist{addtioID,2} = eval(addtio);
            ifaddtio = true;
         elseif ~ischar(addtio) && ischar(oldaddtio)
            error('Additional time integrator option %u must be a string!!!',...
               addtioID);
         end;
      end;
   else
      for addtioID = 1:runopts.addtio
         runopts.addtiolist{addtioID,2} = eval(runopts.addtiolist{addtioID,2});
      end;
      ifaddtio = true;
   end;
end;

% Evaluating the additional solver options
ifaddso = false;
if isfield(runopts,'addsolist')
   if isfield(oldrunopts,'addsolist') && runopts.addso == oldrunopts.addso
      for addsoID = 1:runopts.addso
         oldaddso = oldrunopts.addsolist{addsoID,2};
         addso = runopts.addsolist{addsoID,2};
         
         if ischar(addso) && ((ischar(oldaddso) && ...
               ~strcmp(oldaddso,addso)) || ~ischar(oldaddso))
            runopts.addsolist{addsoID,2} = eval(addso);
            ifaddso = true;
         elseif ~ischar(addso) && ischar(oldaddso)
            error('Additional solver option %u must be a string!!!',addsoID);
         end;
      end;
   else
      for addsoID = 1:runopts.addso
         runopts.addsolist{addsoID,2} = eval(runopts.addsolist{addsoID,2});
      end;
      ifaddso = true;
   end;
end;

% Evaluating the additional jumpstart solver options
ifaddjso = false;
if isfield(runopts,'addjsolist')
   if isfield(oldrunopts,'addjsolist') && runopts.addjso == oldrunopts.addjso
      for addjsoID = 1:runopts.addjso
         oldaddjso = oldrunopts.addjsolist{addjsoID,2};
         addjso = runopts.addjsolist{addjsoID,2};
         
         if ischar(addjso) && ((ischar(oldaddjso) && ...
               ~strcmp(oldaddjso,addjso)) || ~ischar(oldaddjso))
            runopts.addjsolist{addjsoID,2} = eval(addjso);
            ifaddjso = true;
         elseif ~ischar(addjso) && ischar(oldaddjso)
            error(['Additional jumpstart solver option %u must be a ' ...
               'string!!!'],addjsoID);
         end;
      end;
   else
      for addjsoID = 1:runopts.addjso
         runopts.addjsolist{addjsoID,2} = eval(runopts.addjsolist{addjsoID,2});
      end;
      ifaddjso = true;
   end;
end;
end

%% Root directory-based full path construction function

function fullfpath = buildrdfpath(fpath,rootdir)
% The purpose of this function is to build a full path of a file, using a root
% directory for the file path if it is relative.
%
% INPUTS/OUTPUTS
% fpath - File path (may be full or relative) to make full.
% rootdir - Root directory off of which the file path is to be computed if it is
%           relative.
% fullfpath - Full file path.

% Checking if the file path is relative
if ~isfullpath(fpath)
   rootdir = trimfilesep(rootdir);
   fpath = fullfile(rootdir,fpath);
end;

fullfpath = GetFullPath(fpath);
end

%% Relative path construction function to a parent directory

function relpath = buildrpath(cpath,ppath,iffcpath)
% The purpose of this function is to build a relative path from a child
% file/directory to a parent directory. This function will not work with two
% paths that are not a (parent,child) path pair. '..' is used as the pointer to
% a parent folder, which will work on Windows, Unix and Mac machines. The parent
% folder does not have to be an immediate parent.
%
% INPUTS/OUTPUTS
% cpath - Child path. May be file or directory.
% ppath - Parent path. Must be directory and must be (not necessarily immediate)
%         parent to cpath.
% iffcpath - Whether or not cpath is a file or a directory path. OPTIONAL.
%            Default value: true.
% relpath - Relative path from child to parent.

% Organizing input
if nargin == 2
   iffcpath = true;
elseif nargin == 3
   % iffcpath given in an argument
else
   error('Unexpected number of inputs!!!');
end;

% Converting paths to full paths
cpath = GetFullPath(cpath);
ppath = GetFullPath(ppath);

% Initializing relative path
rpnum = 0;

% Trimming paths
cpath = trimfilesep(cpath);
ppath = trimfilesep(ppath);

% Single parent directory transition sequence
pdseq = [filesep '..'];

% Checking path match
if strcmp(cpath,ppath)
   % Paths are already identical, so relative path must be empty
   relpath = '';
else
   % Constructing relative path
   while ~strcmp(cpath,ppath)
      cpath = fileparts(cpath);
      cpath = trimfilesep(cpath);
      rpnum = rpnum + 1;
   end;
   % Correcting for file name
   if iffcpath
      rpnum = rpnum - 1;
   end;
   relpath = repmat(pdseq,1,rpnum);
end;
end

%% Path separator character-trimming function

function dirpath = trimfilesep(dirpath)
% The purpose of this function is to trim the filesep character from the end of
% a directory path. This is necessary to compare two full directory paths.
%
% INPUTS/OUTPUTS
% dirpath - Directory path to trim filesep from.

% Checking platform
if ispc && (length(dirpath) == 3) && (dirpath(2) == ':') && ...
      (dirpath(3) == filesep)
   % A drive letter is the directory path, so should not strip
elseif (isunix || ismac) && (length(dirpath) == 1) && (dirpath(1) == filesep)
   % Root directory is the directory path, so should not strip
elseif dirpath(end) == filesep
   dirpath = dirpath(1:(end-1)); % Trimming dirpath
   dirpath = trimfilesep(dirpath); % Recursively trimming dirpath further
end;
end

%% In-full-path root directory path replacement function

function fullpath = reprootpath(sfullpath,srootpath,rootpath)
% The purpose of this function is to replace the as-saved root path in an
% as-saved full path with an as-loaded root path. If the as-saved root path is
% not present in the as-saved full path, the as-saved full path is turned full
% and is filesep-trimmed, but is not otherwise modified.
%
% INPUTS/OUTPUTS
% sfullpath - As-saved full file/directory path.
% srootpath - As-saved root directory path.
% rootpath - As-loaded root directory path.
% fullpath - Updated full file/directory path.

fullpath = GetFullPath(sfullpath);
fullpath = trimfilesep(fullpath);
fullpath = strrep(fullpath,srootpath,rootpath);
end

%% File creation/opening function

function fids = createfiles(runopts)
% The purpose of this function is to check which files to create and/or open,
% create/open those files for writing/appending, and return the fileIDs in a
% data structure.
%
% INPUTS/OUTPUTS
% runopts - Runtime and Recording Options data structure.
% fids - File IDs data sturcture. May be empty, if no files are required.

% Instantiating fids
fids = struct();

% Log file
if any(runopts.verbose == [1,2]) && any(runopts.log == [1,2])
   fids.logfid = fopen(runopts.logfn,'a');
end;

% Runtime results file
if runopts.rtrr
   fids.rtrfid = fopen(runopts.rtrrfn,'w');
end;
end

%% Time integrator options update function

function [tioptions,ifadj] = adjusttioptions(tioptions,fids,runopts,...
   oldrunopts,ifaddtio,ifaddso,ifaddjso,NonNegative,AbsTol)
% The purpose of this function is to adjust the time integrator options data
% structure based on the adjusted run options and the old run options. The
% adjustment logic is as follows:
%  .tdisc needs to be adjusted if one of the following fields is adjusted:
%    runopts.tst
%    runopts.fts (only exists with runopts.tst == [1,2])
%  .sjso needs to be adjusted either if ifaddso or ifaddjso parameters are true,
%  and/or if one of the following fields is adjusted:
%  (Only the corresponding field(s) need to be adjusted)
%    runopts.stolx
%    runopts.stolres
%    runopts.sdisp
%    runopts.srscal
%    runopts.sxscal
%    runopts.sxcent
%    runopts.spc
%    runopts.sigpdct
%    runopts.sjacob
%    runopts.jstype
%    runopts.jsmsnum
%    runopts.jsphxsplit
%    runopts.jstolx
%    runopts.jstolres
%    runopts.jstit
%    runopts.jsdisp
%    runopts.jsrscal
%    runopts.jsxscal
%    runopts.jsxcent
%    runopts.jspc
%    runopts.jsjacob
%  .phxsys needs to be adjusted if one of the following fields is adjusted:
%  (Only phxsys.phxsplit needs to be adjusted)
%    runopts.phxsplit
%  .fids always needs to be adjusted, unless it was and stayed empty.
%  .rtrr needs to be adjusted if one of the following fields is adjusted:
%    runopts.rtrr
%    runopts.rtrrfq (only exists with runopts.rtrr = true)
%  .NonNegative needs to be adjusted if the NonNegative parameter is non-empty.
%  .AbsTol needs to be adjusted if the AbsTol parameter is non-empty.
%  .OutputFcn needs to be adjusted if the following field is adjusted:
%    runopts.tidisp
%  Additional time integrator options need to be adjusted if the ifaddtio
%  parameter is true.
%
% INPUTS/OUTPUTS
% tioptions - Time Integrator options data structure.
% fids - File IDs data structure.
% runopts - Adjusted Runtime and Recording Options data structure.
% oldrunopts - As-saved Runtime and Recording Options data structure.
% ifaddtio - Whether or not additional time integrator options were adjusted.
% ifaddso - Whether or not additional solver options were adjusted.
% ifaddjso - Whether or not additional jumpstart solver options were adjusted.
% NonNegative - Nonnegative indices vector.
% AbsTol - Absolute tolerances vector.
% ifadj - Whether or not Time Integrator options required adjustment.

oldtioptions = tioptions;

% Adjusting time discretization options
if runopts.tst ~= oldrunopts.tst
   tioptions.tdisc.tst = runopts.tst;
end;
if any(oldtioptions.tdisc.tst == [1,2]) && (runopts.fts ~= oldrunopts.fts)
   if runopts.ftsspect == 0
      runopts.fts = diff(runopts.ii) / double(runopts.fts);
   elseif runopts.ftsspect == 1
      % Fixed time step size already specified
   end;
   tioptions.tdisc.dt = runopts.fts;
elseif any(tioptions.tdisc.tst == [1,2]) && ~any(oldtioptions.tdisc.tst == ...
      [1,2])
   if runopts.ftsspect == 0
      runopts.fts = diff(runopts.ii) / double(runopts.fts);
   elseif runopts.ftsspect == 1
      % Fixed time step size already specified
   end;
   tioptions.tdisc.dt = runopts.fts;
end;

% Adjusting solver and jumpstart options
sjsofldnames = {'stolx','stolres','sdisp','srscal','sxscal','sxcent','spc',...
   'sigpdct','sjacob','jstype','jsmsnum','jsphxsplit','jstolx','jstolres',...
   'jstit','jsdisp','jsrscal','jsxscal','jsxcent','jspc','jsjacob'};
for fldid = 1:length(sjsofldnames)
   fld = sjsofldnames{fldid};
   if ~isequal(runopts.(fld),oldrunopts.(fld))
      tioptions.sjso.(fld) = runopts.(fld);
   end;
end;
if ifaddso
   tioptions.sjso.addsolist = runopts.addsolist;
end;
if ifaddjso
   tioptions.sjso.addjsolist = runopts.addjsolist;
end;

% Adjusting physics splitting options
if runopts.phxsplit ~= oldrunopts.phxsplit
   tioptions.phxsys.phxsplit = runopts.phxsplit;
end;

% Adjusting open file options
if ~isempty(fieldnames(fids))
   tioptions.fids = fids;
end;

% Adjusting runtime recording options
if runopts.rtrr ~= oldrunopts.rtrr
   tioptions.rtrr.ifrtrr = runopts.rtrr;
end;
if oldtioptions.rtrr.ifrtrr && isfield(runopts,'rtrrfq') && ...
      isfield(oldrunopts,'rtrrfq') && (runopts.rtrrfq ~= oldrunopts.rtrrfq)
   tioptions.rtrr.rtrrfq = runopts.rtrrfq;
elseif tioptions.rtrr.ifrtrr && ~oldtioptions.rtrr.ifrtrr
   tioptions.rtrr.rtrrfq = runopts.rtrrfq;
end;

% Adjusting nonnegative indices vector
if ~isempty(NonNegative)
   tioptions.NonNegative = NonNegative;
end;

% Adjusting absolute tolerances vector
if ~isempty(AbsTol)
   tioptions.AbsTol = AbsTol;
end;

% Adjusting time integrator display level
if runopts.tidisp ~= oldrunopts.tidisp
   if runopts.tidisp == 0
      tioptions.OutputFcn = [];
   elseif runopts.tidisp == 1
      tioptions.OutputFcn = @tistepdisp;
   end;
end;

% Adjusting additional time integrator options
if ifaddtio
   for addtioID = 1:runopts.addtio
      pname = runopts.addtiolist{addtioID,1};
      pval = runopts.addtiolist{addtioID,2};
      tioptions.(pname) = pval;
   end;
end;

% Checking if adjustment was required
fldnames = {'tdisc','sjso','phxsys','fids','rtrr','NonNegative','AbsTol',...
   'OutputFcn'};
ifadj = ifaddtio;
if ~ifadj
   for fldid = 1:length(fldnames)
      fldname = fldnames{fldid};
      ifadj = ifadj || ~isequal(tioptions.(fldname),oldtioptions.(fldname));
   end;
end;
end

%% Post-processing vector reformulation function

function [ppsys,ifadj,srtsys] = reformppsys(ppsys,srtsys,vars,refs,runopts,...
   oldrunopts)
% The purpose of this function is to check if the post-processing vector data
% structure requires reformulation, and if it does, to reformulate it. The
% vector requires reformulation if any of the adjusted post-processing options
% are different from the as-saved post-processing options. If it requires
% reformulation, the sorting system data structure is first modified, with
% variables to evaluate for post-processing appended to it.
%
% INPUTS/OUTPUTS
% ppsys - Post-processing vector system data structure.
% srtsys - Sorting system data structure.
% vars - Variables data structure.
% refs - Referencing arrays data structure.
% runopts - Adjusted Runtime and Recording Options data structure.
% oldrunopts - As-saved Runtime and Recording Options data structure.
% ifadj - Whether or not post-processing vector required adjustment.

ifadj = ((runopts.epp ~= oldrunopts.epp) || ...
   (runopts.epp > 1 && ~isequal(runopts.eppIDs,oldrunopts.eppIDs)) || ...
   (runopts.fpp ~= oldrunopts.fpp) || ...
   (runopts.fpp > 1 && ~isequal(runopts.fppIDs,oldrunopts.fppIDs)));

if ifadj
   srtsys = getvars2evalpp(srtsys,runopts,refs);
   ppsys = createppsys(srtsys,vars,runopts);
end;

if runopts.ppdisp ~= oldrunopts.ppdisp
   ppsys.ppdisp = runopts.ppdisp;
   ifadj = true;
end;
end