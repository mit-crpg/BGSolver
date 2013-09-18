%% writeBGSD
% BGSD file writing function.
%
% ifwritten = writeBGSD(runopts,bgs,[fname])
%
% This function accepts the run options and bond graph system data structures,
% which are expected to have been completed by the completeBGSDdata function in
% 'writeBGSD' operations mode. It also accepts the BGSD file name. The function
% then either writes the BGSD file to the disk, or first requests the BGSD file
% name from the user (if not given). If the file with the specified name already
% exists, it is overwritten.
%
% Package:    BGSolver v1.03
% Subpackage: BGSD_Generator
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Variable descriptions
%
% All of the fields are required if and only if they are to be written into the
% BGSD file. The function checks which fields to write based on the options
% specified.
%
% INPUTS
% The data structures' specifications for this function are more stringent than
% for completeBGSDdata or for fBGSolve. Differences include data types and
% whether or not runopts' fields can be optional. The data structures'
% specifications for this function are identical to those for readBGSD.
%
% runopts - Run Options data structure. Describes sections 1-7 in the BGSD file.
%           Its fields are listed below.
%  SECTION 1 - File Header.
%  runopts.svc - Software version that created the data structures. double.
%  runopts.svf - Software version that the data structures are intended for.
%                double.
%  runopts.ifnotes - Whether or not a notes section describing the problem is
%                    present. boolean.
%  runopts.notes - A short description of the file. string.
%  SECTION 2 - Runtime and Recording Options.
%  runopts.verbose - How verbose the processor should be during processing.
%                    int8.
%  runopts.rti - Whether or not to run the time integrator after the AEs are
%                sorted. boolean.
%  runopts.log - How to log all the outputs specified by runopts.verbose. int8.
%  runopts.logfnspec - How the log file name is specified. int8.
%  runopts.logfn - Log file name. string.
%  runopts.smr - Whether or not to record the sorted expressions and supporting
%                variables into a .mat file after the algebraic equations are
%                sorted. boolean.
%  runopts.smfnspec - How the sorted .mat file name is specified. int8.
%  runopts.smfn - Sorted .mat file name. string.
%  runopts.rtrr - Whether or not to record the integration results at runtime,
%                 as they are being computed. boolean.
%  runopts.rtrrfqspec - How the runtime results recording frequency is
%                       specified. int8.
%  runopts.rtrrfq - Runtime results recording frequency. uint32.
%  runopts.rtrrfnspec - How the runtime results file name is specified. int8.
%  runopts.rtrrfn - Runtime results file name. string.
%  runopts.rrm - Whether or not to record the final results in a .mat file.
%                boolean.
%  runopts.rrmfnspec - How the final results .mat file name is specified. int8.
%  runopts.rrmfn - Final results .mat file name. string.
%  SECTION 3 - Time Integrator Options.
%  runopts.titspec - How the time integrator type is specified. int8.
%  runopts.tit - Time integrator type. string.
%  runopts.iispec - How the integration interval is specified. int8.
%  runopts.ii - Integration interval. double vector of length 2.
%  runopts.tst - Time step type, either adaptive or fixed. int8.
%  runopts.ftsspect - What the fixed time step is specified in terms of. int8.
%  runopts.ftsspec - How the fixed time step is specified. int8.
%  runopts.fts - Fixed time step. Either uint32 (number of equally-size time
%                steps to take) or double (fixed time step length).
%  runopts.phxsplit - Whether, and how, the physics are split during
%                     integration. int8.
%  runopts.addtio - Number of additional time integrator options. uint8.
%  runopts.addtiolist - [addtio x 2] cell array of additional time integrator
%                       options pairs (parameter name, parameter value).
%                       Parameter name must be a string, parameter value can be
%                       a string (gets evaluated prior to processing) or a
%                       numeric scalar.
%  runopts.tidisp - Time integrator display level. int8.
%  SECTION 4 - Solver and Jumpstart Options.
%  runopts.stolx - Solver tolerance on state vector. double.
%  runopts.stolres - Solver tolerance on residual. double.
%  runopts.sdisp - Solver display level. int8.
%  runopts.srscal - Solver residual scaling. int8.
%  runopts.sxscal - Solver state vector scaling. int8.
%  runopts.sxcent - Solver state vector centering. int8.
%  runopts.spc - Solver preconditioner. int8.
%  runopts.sigpdct - Solver initial guess predictor. int8.
%  runopts.sjacob - Solver Jacobian construction. int8.
%  runopts.sjsp - Solver Jacobian sparsity treatment. int8.
%  runopts.addso - Number of additional solver options. uint8.
%  runopts.addsolist - [addso x 2] cell array of additional solver options pairs
%                      (parameter name, parameter value). Parameter name must be
%                      a string, parameter value can be a string (gets evaluated
%                      prior to processing) or a numeric scalar.
%  runopts.jstype - Jumpstart type. int8.
%  runopts.jsmsnum - Jumpstart ministep number per 1 full time step. uint32.
%  runopts.jsphxsplit - Whether, and how, the physics are split during
%                       jumpstart. int8.
%  runopts.jstolx - Solver tolerance on state vector to use during jumpstart.
%                   double.
%  runopts.jstolres - Solver tolerance on residual to use during jumpstart.
%                     double.
%  runopts.jstit - Time integrator type to use during jumpstart. string.
%  runopts.jsdisp - Jumpstart solver display level. int8.
%  runopts.jsrscal - Jumpstart residual scaling. int8.
%  runopts.jsxscal - Jumpstart state vector scaling. int8.
%  runopts.jsxcent - Jumpstart state vector centering. int8.
%  runopts.jspc - Jumpstart preconditioner. int8.
%  runopts.jsjacob - Jumpstart Jacobian construction. int8.
%  runopts.jsjsp - Jumpstart Jacobian sparsity treatment. int8.
%  runopts.addjso - Number of additional jumpstart solver options. uint8.
%  runopts.addjsolist - [addjso x 2] cell array of additional jumpstart solver
%                       options pairs (parameter name, parameter value).
%                       Parameter name must be a string, parameter value can be
%                       a string (gets evaluated prior to processing) or a
%                       numeric scalar.
%  SECTION 5 - Evaluator Options.
%  runopts.cs2n - Whether or not to convert the symbolic expressions to numeric
%                 expressions after sorting. boolean.
%  runopts.nft - Numeric function type to convert the symbolic expressions to.
%                int8.
%  SECTION 6 - Post-Processor Options.
%  runopts.epp - Number of effort variables to evaluate during post-processing.
%                int32.
%  runopts.eppIDs - List of BIDs of effort variables to evaluate during
%                   post-processing. uint16 vector (BIDs of effort variables).
%  runopts.fpp - Number of flow variables to evaluate during post-processing.
%                int32.
%  runopts.fppIDs - List of BIDs of flow variables to evaluate during
%                   post-processing. uint16 vector (BIDs of effort variables).
%  runopts.ppdisp - Post-processing display level. int8.
%  SECTION 7 - Linked Files and Directories.
%  runopts.lfnum - Number of linked files. uint8.
%  runopts.lflist - List of linked file names. Cell vector of strings.
%  runopts.ldnum - Number of linked directories. uint8.
%  runopts.ldlist - List of linked directory names. Cell vector of strings.
%
% bgs - Bond Graph System data structure. Describes sections 8-10 in the
%       BGSD file. Its fields are listed below.
%  SECTION 8 - Additional System Information.
%  bgs.asi.sci - Specified causality information. string.
%  bgs.asi.smi - Specified modulation information. string.
%  bgs.asi.sli - Specified linearity information. string.
%  bgs.asi.sri - Specified reduction information. string.
%  bgs.asi.sei - Specified error information. string.
%  SECTION 9 - Bond Connectivity Map.
%  bgs.bnum - Number of bonds in the system. uint16.
%  bgs.bonds - Bond directionality table. uint16 [bnum x 2] array, with entries
%              (to,from) in each row.
%  SECTION 10 - Element List.
%  Storage type b is assumed.
%   bgs.enum - Number of elements in the system. uint16.
%   bgs.exprnum - Number of elements with expressions. uint16.
%   bgs.emodnum - Number of modulated elements. uint16.
%   bgs.epambnum - Number of elements with ambiguous numbers of ports. uint16.
%   bgs.ebambnum - Number of elements with ambiguous bond-to-port
%                  connectivities. uint16.
%   bgs.ecambnum - Number of elements with ambiguous causalities on their ports.
%                  uint16.
%   bgs.ennum - Number of elements with numeric expressions. uint16.
%   bgs.nvnum - Number of numeric variables. uint16.
%   bgs.xnum - Number of storage elements with initial values. uint16.
%   bgs.capnum - Number of capacitive elements. uint16.
%   bgs.inertnum - Number of inertial elements. uint16.
%   bgs.phxspecnum - Number of physics-specific elements. uint16.
%   bgs.nnspecnum - Number of non-negativity-specific elements. uint16.
%   bgs.atspecnum - Number of elements with specified absolute tolerances.
%   bgs.elements - List of element objects (data structures), in order of
%                  increasing EIDs. struct vector, one vector element per EID.
%                  An element data structure has the following fields (some of
%                  them may be empty):
%    .etype - Element type. string.
%    .exprtype - Element's expression type. string.
%    .x0 - Element's initial value. double.
%    .epnum - Element's port number. uint16.
%    .ebcon - Element's bond-to-port connectivity. uint16 vector, one BID per
%             vector element.
%    .caus - Element's causality. int8 vector, one input type per vector
%            element.
%    .phx - Element's physics. uint8.
%    .nn - Element's non-negativity status. boolean.
%    .at - Element's absolute tolerance. double.
%    .mvars - Element's modulating variables. Cell vector of strings.
%    .expr - Element's expression(s). Individual expressions may be
%            strings (for all expression types but CCs) or doubles (for CCs).
%            For multiport elements with CC expressions the coefficients are
%            stored as double matrices. For multiport elements with non-CC
%            expressions the expressions are stored as cell arrays of strings.
%            CMC and NMC expression arrays are shaped as matrices. All
%            individual expressions of an element have to be of the same class.
%
% fname - File name of the BGSD file to write the data to. This may be either
%         the full file name (with path), or only the file name itself. string.
%         OPTIONAL.
%
% OUTPUTS
% ifwritten - Whether or not the BGSD file was successfully written. boolean.

%% Master BGSD writing function

function varargout = writeBGSD(varargin)
% The purpose of this function is to take the run options and bond graph system
% data structures, and either to write them to a BGSD file with a specified
% name, or request the file name from the user and then write the data
% structures.
%
% INPUTS/OUTPUTS
% runopts - Runtime and Recording Options data structure.
% bgs - Bond Graph System data structure.
% fname - BGSD file name. OPTIONAL.
% ifwritten - Whether or not the BGSD file was successfully written.
% See the variable description cell for more information.

% Constructing input parser
ip = buildIP();
% Parsing input
ip.parse(varargin{:});
% Retrieving inputs
runopts = ip.Results.runopts;
bgs = ip.Results.bgs;
if isempty(ip.Results.fname)
   fname = savefiledlg({'*.bgsd','Bond Graph System Descriptor (*.bgsd)'},...
      'Save BGSD file as...');
else
   fname = ip.Results.fname;
end;
% Clearing memory
clear('ip');

% Opening the file for writing
fid = fopen(fname,'w');

% Constructing cleanup routine
cleanupObj = onCleanup(@()cleanup(fid));

% Writing the file sections
% SECTION 1 - File Header.
write1partline(fid,'BGSD');
write2partline(fid,'SVF',runopts.svf);
write2partline(fid,'SVC',runopts.svc);
write2partline(fid,'NOTES',runopts.ifnotes);
writeopt1partline(fid,runopts,'notes',true);
write1partline(fid,'NOTESEND');

% SECTION 2 - Runtime and Recording Options.
write1partline(fid,'RRO');
write2partline(fid,'VERBOSE',runopts.verbose);
write2partline(fid,'RTI',runopts.rti);
writeopt2partline(fid,'LOG',runopts,'log');
writeopt2partline(fid,'LOGFNSPEC',runopts,'logfnspec');
writeopt2partline(fid,'LOGFN',runopts,'logfn',true);
write2partline(fid,'SMR',runopts.smr);
writeopt2partline(fid,'SMFNSPEC',runopts,'smfnspec');
writeopt2partline(fid,'SMFN',runopts,'smfn',true);
write2partline(fid,'RTRR',runopts.rtrr);
writeopt2partline(fid,'RTRRFQSPEC',runopts,'rtrrfqspec');
writeopt2partline(fid,'RTRRFQ',runopts,'rtrrfq');
writeopt2partline(fid,'RTRRFNSPEC',runopts,'rtrrfnspec');
writeopt2partline(fid,'RTRRFN',runopts,'rtrrfn',true);
write2partline(fid,'RRM',runopts.rrm);
writeopt2partline(fid,'RRMFNSPEC',runopts,'rrmfnspec');
writeopt2partline(fid,'RRMFN',runopts,'rrmfn',true);
write1partline(fid,'RROEND');

% SECTION 3 - Time Integrator Options.
write1partline(fid,'TIO');
write2partline(fid,'TITSPEC',runopts.titspec);
writeopt2partline(fid,'TIT',runopts,'tit');
write2partline(fid,'IISPEC',runopts.iispec);
writeopt2partline(fid,'II',runopts,'ii');
write2partline(fid,'TST',runopts.tst);
writeopt2partline(fid,'FTSSPECT',runopts,'ftsspect');
writeopt2partline(fid,'FTSSPEC',runopts,'ftsspec');
writeopt2partline(fid,'FTS',runopts,'fts');
write2partline(fid,'PHXSPLIT',runopts.phxsplit);
write2partline(fid,'ADDTIO',runopts.addtio);
for addtioID = 1:runopts.addtio
   write2partline(fid,runopts.addtiolist{addtioID,1},...
      runopts.addtiolist{addtioID,2},true);
end;
write2partline(fid,'TIDISP',runopts.tidisp);
write1partline(fid,'TIOEND');

% SECTION 4 - Solver and Jumpstart Options.
write1partline(fid,'SJSO');
write2partline(fid,'STOLX',runopts.stolx);
write2partline(fid,'STOLRES',runopts.stolres);
write2partline(fid,'SDISP',runopts.sdisp);
write2partline(fid,'SRSCAL',runopts.srscal);
write2partline(fid,'SXSCAL',runopts.sxscal);
write2partline(fid,'SXCENT',runopts.sxcent);
write2partline(fid,'SPC',runopts.spc);
write2partline(fid,'SIGPDCT',runopts.sigpdct);
write2partline(fid,'SJACOB',runopts.sjacob);
write2partline(fid,'SJSP',runopts.sjsp);
write2partline(fid,'ADDSO',runopts.addso);
for addsoID = 1:runopts.addso
   write2partline(fid,runopts.addsolist{addsoID,1},...
      runopts.addsolist{addsoID,2},true);
end;
write2partline(fid,'JSTYPE',runopts.jstype);
write2partline(fid,'JSMSNUM',runopts.jsmsnum);
write2partline(fid,'JSPHXSPLIT',runopts.jsphxsplit);
write2partline(fid,'JSTOLX',runopts.jstolx);
write2partline(fid,'JSTOLRES',runopts.jstolres);
write2partline(fid,'JSTIT',runopts.jstit);
write2partline(fid,'JSDISP',runopts.jsdisp);
write2partline(fid,'JSRSCAL',runopts.jsrscal);
write2partline(fid,'JSXSCAL',runopts.jsxscal);
write2partline(fid,'JSXCENT',runopts.jsxcent);
write2partline(fid,'JSPC',runopts.jspc);
write2partline(fid,'JSJACOB',runopts.jsjacob);
write2partline(fid,'JSJSP',runopts.jsjsp);
write2partline(fid,'ADDJSO',runopts.addjso);
for addjsoID = 1:runopts.addjso
   write2partline(fid,runopts.addjsolist{addjsoID,1},...
      runopts.addjsolist{addjsoID,2},true);
end;
write1partline(fid,'SJSOEND');

% SECTION 5 - Evaluator Options.
write1partline(fid,'EO');
write2partline(fid,'CS2N',runopts.cs2n);
writeopt2partline(fid,'NFT',runopts,'nft');
write1partline(fid,'EOEND');

% SECTION 6 - Post-Processor Options.
write1partline(fid,'PPO');
write2partline(fid,'EPP',runopts.epp);
writeopt1partlist(fid,runopts,'eppIDs');
write2partline(fid,'FPP',runopts.fpp);
writeopt1partlist(fid,runopts,'fppIDs');
write2partline(fid,'PPDISP',runopts.ppdisp);
write1partline(fid,'PPOEND');

% SECTION 7 - Linked Files and Directories.
write1partline(fid,'LFD');
write2partline(fid,'LF',runopts.lfnum);
writeopt1partlist(fid,runopts,'lflist',true);
write2partline(fid,'LD',runopts.ldnum);
writeopt1partlist(fid,runopts,'ldlist',true);
write1partline(fid,'LFDEND');

% SECTION 8 - Additional System Information.
write1partline(fid,'ASI');
write2partline(fid,'SCI',bgs.asi.sci);
write2partline(fid,'SMI',bgs.asi.smi);
write2partline(fid,'SLI',bgs.asi.sli);
write2partline(fid,'SRI',bgs.asi.sri);
write2partline(fid,'SEI',bgs.asi.sei);
write1partline(fid,'ASIEND');

% SECTION 9 - Bond Connectivity Map.
write2partline(fid,'BCM',bgs.bnum);
writemat(fid,bgs.bonds);
write1partline(fid,'BCMEND');

% SECTION 10 - Element List.
write2partline(fid,'EL',bgs.enum);
write2partline(fid,'EEXPRS',bgs.exprnum);
write2partline(fid,'EMODS',bgs.emodnum);
write2partline(fid,'EPAMBS',bgs.epambnum);
write2partline(fid,'EBAMBS',bgs.ebambnum);
write2partline(fid,'ECAMBS',bgs.ecambnum);
write2partline(fid,'ENUMS',bgs.ennum);
write2partline(fid,'NVARS',bgs.nvnum);
write2partline(fid,'EV0S',bgs.xnum);
write2partline(fid,'ECAPS',bgs.capnum);
write2partline(fid,'EINERTS',bgs.inertnum);
write2partline(fid,'EPHXSPEC',bgs.phxspecnum);
write2partline(fid,'ENNSPEC',bgs.nnspecnum);
write2partline(fid,'EATSPEC',bgs.atspecnum);
for eid = 1:bgs.enum
   write2partline(fid,'EID',eid);
   elem = bgs.elements(eid);
   write2partline(fid,'ET',elem.etype);
   ifexpr = writeopt2partline(fid,'EXPRT',elem,'exprtype');
   if ifexpr
      writeopt2partline(fid,'V0',elem,'x0');
      writeopt2partline(fid,'EP',elem,'epnum');
      writeopt2partline(fid,'EB',elem,'ebcon');
      writeopt2partline(fid,'EC',elem,'caus');
      writeopt2partline(fid,'EPHX',elem,'phx');
      writeopt2partline(fid,'ENN',elem,'nn');
      writeopt2partline(fid,'EAT',elem,'at');
      writeopt2partline(fid,'EMVARS',elem,'mvars',true,', ');
      writeexpr(fid,elem.etype,elem.exprtype,elem.expr);
   end;
   write1partline(fid,'EEND');
end;
write1partline(fid,'ELEND');

% SECTION 11 - File End.
write1partline(fid,'BGSDEND');

% Assuming writing is successful
if nargout == 0
   % No return value necessary
elseif nargout == 1
   varargout = true;
else
   error('Invalid number of outputs!!!');
end;

% cleanup function executes here
end

%% Input parser construction function

function ip = buildIP()
% The purpose of this function is to construct an input parser object for the
% writeBGSD function.
%
% INPUTS/OUTPUTS
% ip - Input parser object for the writeBGSD function.

fchkfun = @(fname) islegalfpath(fname,'ext','.bgsd');

ip = inputParser;
ip.StructExpand = false;
ip.addRequired('runopts',@isstruct);
ip.addRequired('bgs',@isstruct);
ip.addOptional('fname',[],fchkfun);
end

%% Cleanup function

function cleanup(fid)
% The purpose of this function is to run the cleanup routines when writeBGSD
% exits, either because of an error, or regularly.
%
% INPUTS/OUTPUTS
% fid - BGSD file ID.

fclose(fid);
end

%% One-part line writing function

function write1partline(fid,val,ifquote,sep)
% The purpose of this function is to write a single 1-string line at the
% location of the cursor, and transition to the next line. The value to write
% can be multiple data types. An optional argument tells the function to
% surround the string in quotes. Another optional argument specifies the
% separator sequence of characters if the value is an array.
%
% INPUTS
% fid - BGSD file ID.
% val - Entry to append at the end of the BGSD file.
% ifquote - Whether or not to surround the string in quotes. OPTIONAL. Default:
%           false.
% sep - Sequence of characters to use as separator of individual elements of the
%       array. OPTIONAL. Default: ' ' (whitespace).

if nargin == 2
   ifquote = false;
   sep = ' ';
elseif nargin == 3
   % ifquote specified in the input
   sep = ' ';
elseif nargin == 4
   % ifquote specified in the input
   % sep specified in the input
else
   error('Invalid number of inputs!!!');
end;

str = field2str(val,sep);

if ifquote
   str = quote(str);
end;

fullstr = [str,'\n'];
fprintf(fid,fullstr);
end

%% Two-part line writing function

function write2partline(fid,kstr,val,ifquote,sep)
% The purpose of this function is to write a single 2-part line at the location
% of the cursor, and transition to the next line. The first entry is the
% keyword, and has to be a string. The second entry can be multiple data types.
% An optional argument tells the function to surround the second entry in
% quotes. Another optional argument specifies the separator sequence of
% characters if the second entry is an array.
%
% INPUTS
% fid - BGSD file ID.
% kstr - Keyword string.
% val - Value after the keyword (may be string, number, etc.).
% ifquote - Whether or not to surround the string in quotes. OPTIONAL. Default:
%           false.
% sep - Sequence of characters to use as separator of individual elements of the
%       array. OPTIONAL. Default: ' ' (whitespace).

if nargin == 3
   ifquote = false;
   sep = ' ';
elseif nargin == 4
   % ifquote specified in the input
   sep = ' ';
elseif nargin == 5
   % ifquote specified in the input
   % sep specified in the input
else
   error('Invalid number of inputs!!!');
end;

% Converting the value to a string
valstr = field2str(val,sep);

% Quoting the string, if necessary
if ifquote
   valstr = quote(valstr);
end;

% Constructing the full line
S = [kstr ' ' valstr];

% Writing the line
write1partline(fid,S);
end

%% Optional one-part line writing function

function ifwritten = writeopt1partline(fid,ds,fieldname,ifquote,sep)
% The purpose of this function is to check if a data structure has a given
% field, and if it does, to write the field in a single 1-string line at the
% location of the cursor, and transition to the next line. The field can be
% multiple data types. An optional argument tells the function to surround the
% field in quotes. Another optional argument specifies the separator sequence of
% characters if the field is an array.
%
% INPUTS/OUTPUTS
% fid - BGSD file ID.
% ds - Data structure which may have the specified field.
% fieldname - Name of the field to check and write, if present.
% ifquote - Whether or not to surround the string in quotes. OPTIONAL. Default:
%           false.
% sep - Sequence of characters to use as separator of individual elements of the
%       array. OPTIONAL. Default: ' ' (whitespace).
% ifwritten - Whether or not the field was present in the data structure.

if isfield(ds,fieldname) && ~isempty(ds.(fieldname))
   if nargin == 3
      ifquote = false;
      sep = ' ';
   elseif nargin == 4
      % ifquote specified in the input
      sep = ' ';
   elseif nargin == 5
      % ifquote specified in the input
      % sep specified in the input
   else
      error('Invalid number of inputs!!!');
   end;
   write1partline(fid,ds.(fieldname),ifquote,sep);
   ifwritten = true;
else
   ifwritten = false;
end;
end

%% Optional two-part line writing function

function ifwritten = writeopt2partline(fid,kstr,ds,fieldname,ifquote,sep)
% The purpose of this function is to check if a data structure has a given
% field, and if it does, to write the field, following a keyword, as a 2-part
% line at the location of the cursor, and transition to the next line. The
% keyword has to be a string. The field can be multiple data types. An optional
% argument tells the function to surround the field in quotes. Another optional
% argument specifies the separator sequence of characters if the second entry is
% an array.
%
% INPUTS/OUTPUTS
% fid - BGSD file ID.
% kstr - Keyword string.
% ds - Data structure which may have the specified field.
% fieldname - Name of the field to check and write, if present.
% ifquote - Whether or not to surround the string in quotes. OPTIONAL. Default:
%           false.
% sep - Sequence of characters to use as separator of individual elements of the
%       array. OPTIONAL. Default: ' ' (whitespace).
% ifwritten - Whether or not the field was present in the data structure.

if isfield(ds,fieldname) && ~isempty(ds.(fieldname))
   if nargin == 4
      ifquote = false;
      sep = ' ';
   elseif nargin == 5
      % ifquote specified in the input
      sep = ' ';
   elseif nargin == 6
      % ifquote specified in the input
      % sep specified in the input
   else
      error('Invalid number of inputs!!!');
   end;
   write2partline(fid,kstr,ds.(fieldname),ifquote,sep);
   ifwritten = true;
else
   ifwritten = false;
end;
end

%% One-part multiline list writing function

function write1partlist(fid,list,ifquote)
% The purpose of this function is to write a vector as a list of 1-string lines
% at the location of the cursor, and transition to the next line. The vector can
% be either a non-cell array, or a cell array of strings. An optional argument
% tells the function to surround the string in quotes.
%
% INPUTS
% fid - BGSD file ID.
% list - Data vector to append at the end of the BGSD file. It can be either a
%        non-cell array, or a cell array of strings.
% ifquote - Whether or not to surround the string in quotes. OPTIONAL. Default:
%           false.

if nargin == 2
   ifquote = false;
elseif nargin == 3
   % ifquote specified in the input
else
   error('Invalid number of inputs!!!');
end;

assert(isvector(list),'List has to be a vector!!!');

if iscellstr(list)
   for i = 1:length(list)
      write1partline(fid,list{i},ifquote);
   end;
else
   fmt = getlistfmt(list,ifquote);
   fullfmt = [fmt '\n'];
   
   fprintf(fid,fullfmt,list);
end;
end

%% Optional one-part multiline list writing function

function ifwritten = writeopt1partlist(fid,ds,fieldname,ifquote)
% The purpose of this function is to check if a data structure has a given
% field, and if it does, to write the field as a multiline list at the location
% of the cursor, and then transition to the next line. The field can be either a
% non-cell array, or a cell array of strings. An optional argument tells the
% function to surround the line entries in quotes.
%
% INPUTS/OUTPUTS
% fid - BGSD file ID.
% ds - Data structure which may have the specified field.
% fieldname - Name of the field to check and write, if present.
% ifquote - Whether or not to surround the string in quotes. OPTIONAL. Default:
%           false.
% ifwritten - Whether or not the field was present in the data structure.

if isfield(ds,fieldname) && ~isempty(ds.(fieldname))
   assert(isvector(ds.(fieldname)),'Only vector lists are supported!!!');
   if nargin == 3
      ifquote = false;
   elseif nargin == 4
      % ifquote specified in the input
   else
      error('Invalid number of inputs!!!');
   end;
   write1partlist(fid,ds.(fieldname),ifquote);
   ifwritten = true;
else
   ifwritten = false;
end;
end

%% Expression writing function

function writeexpr(fid,etype,exprtype,expr)
% The purpose of this function is to write the expression(s) of given expression
% and element types at the location of the cursor, and transition to the next
% line.
%
% INPUTS
% fid - BGSD file ID.
% etype - Element type.
% exprtype - Expression type.
% expr - Expression, can be a number, an array of numbers, a string or a cell
%        vector/array of strings.

if any(strcmp(etype,{'SE','SF','I','C','R','TF','GY','1','0','MSE','MSF',...
      'MI','MC','MR','MTF','MGY'}))
   switch exprtype
      case {'CC','NMC','NE','NME'}
         write1partline(fid,expr,false);
         
      case {'CMC','CE','CME'}
         write1partline(fid,expr,true);
         
      otherwise
         error('Unknown expression type!!!');
   end;
elseif any(strcmp(etype,{'R2','MR2','RN','MRN'}))
   switch exprtype
      case 'CC'
         writemat(fid,expr);
         
      case 'CMC'
         writemataslist(fid,expr,true);
         
      case 'NMC'
         writemataslist(fid,expr,false);
         
      case 'CE'
         write1partlist(fid,expr,true);
         
      case 'NE'
         write1partlist(fid,expr,false);
         
      case 'CME'
         write1partlist(fid,expr,true);
         
      case 'NME'
         write1partlist(fid,expr,false);
         
      otherwise
         error('Unknown expression type!!!');
   end;
else
   error('Unknown element type!!!');
end;
end

%% Matrix writing function

function writemat(fid,mat)
% The purpose of this function is to write a 2-dimensional array at the location
% of the cursor, and transition to the next line. The 2-d array has to be
% non-cell.
%
% INPUTS
% fid - BGSD file ID.
% mat - 2-dimensional array to write into the file.

fmt = getlistfmt(mat,false,true,'\t');
fullfmt = [fmt '\n'];
fprintf(fid,fullfmt,mat');
end

%% Matrix-as-list writing function

function writemataslist(fid,mat,ifquote)
% The purpose of this function is to write a 2-dimensional array as a multiline
% list at the location of the cursor, and transition to the next line. The 2-d
% array can be a cell array of strings, or a non-cell array. It is written in
% order of left to right, then top to bottom. An optional argument tells the
% function to surround each line in quotes.
%
% INPUTS
% fid - BGSD file ID.
% mat - 2-dimensional array to write into the file.
% ifquote - Whether or not to surround the string in quotes. OPTIONAL. Default:
%           false.

if nargin == 2
   ifquote = false;
elseif nargin == 3
   % ifquote specified in the input
else
   error('Invalid number of inputs!!!');
end;

mat = mat';
ifcell = false;

switch class(mat)
   case 'cell'
      assert(iscellstr(mat),'Non-string cell array cannot be written!!!');
      fmt1 = '%s';
      ifcell = true;
      
   case {'int8','int16','int32','int64'}
      fmt1 = '%ld';
      
   case {'uint8','uint16','uint32','uint64'}
      fmt1 = '%ld';
      
   case {'single','double'}
      fmt1 = '%-19.12g';
      
   case 'logical'
      fmt1 = '%d';
      
   otherwise
      error('Wrong data type!!!');
end;

if ifquote
   fmt1 = quote(fmt1);
end;
fmt = [fmt1 '\n'];

if ifcell
   fprintf(fid,fmt,mat{:});
else
   fprintf(fid,fmt,mat);
end;
end

%% Value to string conversion function

function str = field2str(val,sep)
% The purpose of this function is to identify an appropriate string format for
% the value given, based on its data type, and convert that value to a string.
% If the value is a vector, the resulting string is a separated line; the
% separator type can be specified.
%
% INPUTS/OUTPUTS
% val - Value (may be string, double, integer, etc.); can be a non-cell vector
%       array.
% sep - Separator character sequence. OPTIONAL. Default: ' ' (space).
% str - String that the value was converted to.

if nargin == 1
   sep = ' ';
elseif nargin == 2
   % sep specified in the input
else
   error('Invalid number of inputs!!!');
end;

% Constructing single-element format sequence
switch class(val)
   case {'int8','int16','int32','int64'}
      fmt = '%ld';
      
   case {'uint8','uint16','uint32','uint64'}
      fmt = '%lu';
      
   case {'single','double'}
      fmt = '%-19.12g';
      
   case 'char'
      fmt = '%s';
      
   case 'logical'
      fmt = '%d';
      
   case 'cell'
      assert(iscellstr(val),'Unknown cell array contents!!!');
      fmt = '%s';
      
   otherwise
      error('Wrong data type!!!');
end;

% Constructing full string
sepfmt = [fmt sep];
if ~ischar(val) && isvector(val) && ~iscellstr(val)
   if length(val) > 1
      str1 = sprintf(sepfmt,val(1:(end-1)));
      str2 = sprintf(fmt,val(end));
      str = [str1 str2];
   else
      str = sprintf(fmt,val);
   end;
elseif ischar(val)
   str = sprintf(fmt,val);
elseif iscellstr(val)
   if length(val) > 1
      str1 = sprintf(sepfmt,val{1:(end-1)});
      str2 = sprintf(fmt,val{end});
      str = [str1 str2];
   else
      str = sprintf(fmt,val{1});
   end;
end;
end

%% List format string construction function

function fmt = getlistfmt(val,ifquote,frcmat,sep)
% The purpose of this function is to construct the format string for 1 line of
% the list that results from the array submitted to it. The array may be a
% vector (output as a list, 1 number per line), or a 2-D array (output as a
% matrix, one row per line). An optional argument can tell the function to
% surround the string in quotes. Another optional argument can force the
% function to treat 1-dimensional arrays as matrices. A third optional argument
% specifies the separator sequence of characters, which is used if the value is
% a 2-dimensional array.
%
% INPUTS/OUTPUTS
% val - Array to construct the list format for (may be at most 2-dimensional).
% ifquote - Whether or not to surround the string in quotes. OPTIONAL. Default:
%           false.
% frcmat - Whether or not to force the function to treat 1-dimensional arrays as
%          matrices. OPTIONAL. Default: false.
% sep - Sequence of characters of characters to use as separator of individual
%       elements of the array. OPTIONAL. Default: '\t' (tab).
% fmt - Format string to be used to output the array as a list. This string does
%       not include the line termination character.

% Setting defaults
if nargin == 1
   ifquote = false;
   frcmat = false;
   sep = '\t';
elseif nargin == 2
   % ifquote specified in the input
   frcmat = false;
   sep = '\t';
elseif nargin == 3
   % ifquote specified in the input
   % frcmat specified in the input
   sep = '\t';
elseif nargin == 4
   % ifquote specified in the input
   % frcmat specified in the input
   % sep specified in the input
else
   error('Invalid number of inputs!!!');
end;

% Constructing single element's format
switch class(val)
   case {'int8','int16','int32','int64'}
      fmt1 = '%ld';
      
   case {'uint8','uint16','uint32','uint64'}
      fmt1 = '%lu';
      
   case {'single','double'}
      fmt1 = '%.15g';
      
   case 'logical'
      fmt1 = '%d';
      
   otherwise
      error('Wrong data type!!!');
end;

% Constructing full format string
if frcmat || ~isvector(val)
   nparts = size(val,2);
   fmt1sep = [fmt1,sep];
   fullfmt1 = repmat(fmt1sep,1,nparts-1);
   fullfmt2 = fmt1;
   fmt = [fullfmt1 fullfmt2];
else
   fmt = fmt1;
end;

% Surround full format string in quotes, if necessary
if ifquote
   fmt = quote(fmt);
end;
end

%% String quotation function

function qstr = quote(str)
% The purpose of this function is to surround a string with double quotes.
%
% INPUTS/OUTPUTS
% str - String.
% qstr - String surrounded with double quotation marks.

qstr = ['"' str '"'];
end