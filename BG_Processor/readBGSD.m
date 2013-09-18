%% readBGSD
% BGSD file reading function.
%
% [runopts,bgs] = readBGSD(fname)
%
% This function accepts the file name of the BGSD file to read, and returns the
% two data structures the file contains, runopts and bgs. The element list in
% bgs is returned in type b storage, which corresponds to it being completed by
% the completeBGSDdata function in 'writeBGSD' operations mode.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Variable descriptions
%
% INPUTS
% fname - File name of the BGSD file to read. This may be either the full file
%         name (with path), or only the file name itself. string. REQUIRED.
%
% OUTPUTS
% The data structures' specifications for this function are more stringent than
% for completeBGSDdata or for fBGSolve. Differences include data types and
% whether or not runopts' fields can be optional. The data structures'
% specifications for this function are identical to those for writeBGSD.
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

%% Master BGSD file reading function

function [runopts,bgs] = readBGSD(varargin)
% The purpose of this function is to accept the file name of the BGSD file to
% read, and return the two data structures the file contains, runopts and bgs.
% The element list in bgs is returned in type b storage, which corresponds to it
% being completed by the completeBGSDdata function in 'writeBGSD' operations
% mode.
%
% INPUTS/OUTPUTS
% fname - BGSD file name.
% runopts - Runtime and Recording Options data structure.
% bgs - Bond Graph System data structure.
% See the variable description cell for more information.

% Constructing input parser
ip = buildIP();
% Parsing input
ip.parse(varargin{:});
% Retrieving inputs
fname = ip.Results.fname;
% Clearing memory
clear('ip');

% Opening the file for reading
fid = fopen(fname,'r');

% Constructing cleanup routine
cleanupObj = onCleanup(@()cleanup(fid));

% Reading the file sections
% SECTION 1 - File Header.
read1partline(fid,'%s','BGSD');
runopts.svf = read2partline(fid,'%n','SVF');
runopts.svc = read2partline(fid,'%n','SVC');
runopts.ifnotes = read2partline(fid,'%d','NOTES',true);
if runopts.ifnotes
   runopts.notes = read1partline(fid,'%q');
end;
read1partline(fid,'%s','NOTESEND');

% SECTION 2 - Runtime and Recording Options.
read1partline(fid,'%s','RRO');
runopts.verbose = read2partline(fid,'%d8','VERBOSE');
runopts.rti = read2partline(fid,'%d','RTI',true);
if runopts.verbose ~= 0
   runopts.log = read2partline(fid,'%d8','LOG');
   if runopts.log ~= 0
      runopts.logfnspec = read2partline(fid,'%d8','LOGFNSPEC');
      if runopts.logfnspec == 1
         runopts.logfn = read2partline(fid,'%q','LOGFN');
      end;
   end;
end;
runopts.smr = read2partline(fid,'%d','SMR',true);
if runopts.smr
   runopts.smfnspec = read2partline(fid,'%d8','SMFNSPEC');
   if runopts.smfnspec == 1
      runopts.smfn = read2partline(fid,'%q','SMFN');
   end;
end;
runopts.rtrr = read2partline(fid,'%d','RTRR',true);
if runopts.rtrr
   runopts.rtrrfqspec = read2partline(fid,'%d8','RTRRFQSPEC');
   if runopts.rtrrfqspec == 1
      runopts.rtrrfq = read2partline(fid,'%u32','RTRRFQ');
   end;
   runopts.rtrrfnspec = read2partline(fid,'%d8','RTRRFNSPEC');
   if runopts.rtrrfnspec == 1
      runopts.rtrrfn = read2partline(fid,'%q','RTRRFN');
   end;
end;
runopts.rrm = read2partline(fid,'%d','RRM',true);
if runopts.rrm
   runopts.rrmfnspec = read2partline(fid,'%d8','RRMFNSPEC');
   if runopts.rrmfnspec == 1
      runopts.rrmfn = read2partline(fid,'%q','RRMFN');
   end;
end;
read1partline(fid,'%s','RROEND');

% SECTION 3 - Time Integrator Options.
read1partline(fid,'%s','TIO');
runopts.titspec = read2partline(fid,'%d8','TITSPEC');
if runopts.titspec == 1
   runopts.tit = read2partline(fid,'%s','TIT');
end;
runopts.iispec = read2partline(fid,'%d8','IISPEC');
if runopts.iispec == 0
   runopts.ii = read2partline(fid,'%n %n','II');
end;
runopts.tst = read2partline(fid,'%d8','TST');
if any(runopts.tst == [1,2])
   runopts.ftsspect = read2partline(fid,'%d8','FTSSPECT');
   runopts.ftsspec = read2partline(fid,'%d8','FTSSPEC');
   if runopts.ftsspec == 0
      if runopts.ftsspect == 0
         runopts.fts = read2partline(fid,'%u32','FTS');
      elseif runopts.ftsspect == 1
         runopts.fts = read2partline(fid,'%n','FTS');
      else
         error('Invalid FTSSPECT!!!');
      end;
   end;
end;
runopts.phxsplit = read2partline(fid,'%d8','PHXSPLIT');
runopts.addtio = read2partline(fid,'%u8','ADDTIO');
runopts.addtiolist = cell(runopts.addtio,2);
for addtioID = 1:runopts.addtio
   [runopts.addtiolist{addtioID,:}] = readpvalpair(fid,'%q');
end;
runopts.tidisp = read2partline(fid,'%d8','TIDISP');
read1partline(fid,'%s','TIOEND');

% SECTION 4 - Solver and Jumpstart Options.
read1partline(fid,'%s','SJSO');
runopts.stolx = read2partline(fid,'%n','STOLX');
runopts.stolres = read2partline(fid,'%n','STOLRES');
runopts.sdisp = read2partline(fid,'%d8','SDISP');
runopts.srscal = read2partline(fid,'%d8','SRSCAL');
runopts.sxscal = read2partline(fid,'%d8','SXSCAL');
runopts.sxcent = read2partline(fid,'%d8','SXCENT');
runopts.spc = read2partline(fid,'%d8','SPC');
runopts.sigpdct = read2partline(fid,'%d8','SIGPDCT');
runopts.sjacob = read2partline(fid,'%d8','SJACOB');
runopts.sjsp = read2partline(fid,'%d8','SJSP');
runopts.addso = read2partline(fid,'%u8','ADDSO');
for addsoID = 1:runopts.addso
   [runopts.addsolist{addsoID,:}] = readpvalpair(fid,'%q');
end;
runopts.jstype = read2partline(fid,'%d8','JSTYPE');
runopts.jsmsnum = read2partline(fid,'%u32','JSMSNUM');
runopts.jsphxsplit = read2partline(fid,'%d8','JSPHXSPLIT');
runopts.jstolx = read2partline(fid,'%n','JSTOLX');
runopts.jstolres = read2partline(fid,'%n','JSTOLRES');
runopts.jstit = read2partline(fid,'%s','JSTIT');
runopts.jsdisp = read2partline(fid,'%d8','JSDISP');
runopts.jsrscal = read2partline(fid,'%d8','JSRSCAL');
runopts.jsxscal = read2partline(fid,'%d8','JSXSCAL');
runopts.jsxcent = read2partline(fid,'%d8','JSXCENT');
runopts.jspc = read2partline(fid,'%d8','JSPC');
runopts.jsjacob = read2partline(fid,'%d8','JSJACOB');
runopts.jsjsp = read2partline(fid,'%d8','JSJSP');
runopts.addjso = read2partline(fid,'%u8','ADDJSO');
for addjsoID = 1:runopts.addjso
   [runopts.addsolist{addjsoID,:}] = readpvalpair(fid,'%q');
end;
read1partline(fid,'%s','SJSOEND');

% SECTION 5 - Evaluator Options.
read1partline(fid,'%s','EO');
runopts.cs2n = read2partline(fid,'%d','CS2N',true);
if runopts.cs2n
   runopts.nft = read2partline(fid,'%d8','NFT');
end;
read1partline(fid,'%s','EOEND');

% SECTION 6 - Post-Processor Options.
read1partline(fid,'%s','PPO');
runopts.epp = read2partline(fid,'%d32','EPP');
if runopts.epp > 0
   runopts.eppIDs = read1partlist(fid,'%u16',runopts.epp);
end;
runopts.fpp = read2partline(fid,'%d32','FPP');
if runopts.fpp > 0
   runopts.fppIDs = read1partlist(fid,'%u16',runopts.fpp);
end;
runopts.ppdisp = read2partline(fid,'%d8','PPDISP');
read1partline(fid,'%s','PPOEND');

% SECTION 7 - Linked Files and Directories.
read1partline(fid,'%s','LFD');
runopts.lfnum = read2partline(fid,'%u8','LF');
if runopts.lfnum > 0
   runopts.lflist = read1partlist(fid,'%q',runopts.lfnum);
end;
runopts.ldnum = read2partline(fid,'%u8','LD');
if runopts.ldnum > 0
   runopts.ldlist = read1partlist(fid,'%q',runopts.ldnum);
end;
read1partline(fid,'%s','LFDEND');

% SECTION 8 - Additional System Information.
read1partline(fid,'%s','ASI');
bgs.asi.sci = read2partline(fid,'%s','SCI');
bgs.asi.smi = read2partline(fid,'%s','SMI');
bgs.asi.sli = read2partline(fid,'%s','SLI');
bgs.asi.sri = read2partline(fid,'%s','SRI');
bgs.asi.sei = read2partline(fid,'%s','SEI');
read1partline(fid,'%s','ASIEND');

% SECTION 9 - Bond Connectivity Map.
bgs.bnum = read2partline(fid,'%u16','BCM');
bgs.bonds = readmat(fid,'%u16',2,bgs.bnum);
read1partline(fid,'%s','BCMEND');

% SECTION 10 - Element List.
bgs.enum = read2partline(fid,'%u16','EL');
bgs.exprnum = read2partline(fid,'%u16','EEXPRS');
bgs.emodnum = read2partline(fid,'%u16','EMODS');
bgs.epambnum = read2partline(fid,'%u16','EPAMBS');
bgs.ebambnum = read2partline(fid,'%u16','EBAMBS');
bgs.ecambnum = read2partline(fid,'%u16','ECAMBS');
bgs.ennum = read2partline(fid,'%u16','ENUMS');
bgs.nvnum = read2partline(fid,'%u16','NVARS');
bgs.xnum = read2partline(fid,'%u16','EV0S');
bgs.capnum = read2partline(fid,'%u16','ECAPS');
bgs.inertnum = read2partline(fid,'%u16','EINERTS');
bgs.phxspecnum = read2partline(fid,'%u16','EPHXSPEC');
bgs.nnspecnum = read2partline(fid,'%u16','ENNSPEC');
bgs.atspecnum = read2partline(fid,'%u16','EATSPEC');
bgs.elements = preallocEL(bgs.enum);
for eid = 1:bgs.enum
   eidread = read2partline(fid,'%u16','EID');
   assert(eid == eidread,'Expected EID %lu, got EID %lu instead!!!',eid,...
      eidread);
   
   bgs.elements(eid).etype = read2partline(fid,'%s','ET');
   efields = getefields(bgs.elements(eid).etype);
   
   if efields.ifexpr
      bgs.elements(eid).exprtype = read2partline(fid,'%s','EXPRT');
      efields = getefields(bgs.elements(eid).etype,bgs.elements(eid).exprtype);
      
      if efields.ifx0
         bgs.elements(eid).x0 = read2partline(fid,'%n','V0');
      end;
      
      if efields.ifpnum
         bgs.elements(eid).epnum = read2partline(fid,'%u16','EP');
      end;
      
      if efields.ifbcon
         fmt = buildfmt('%u16',bgs.elements(eid).epnum);
         bgs.elements(eid).ebcon = read2partline(fid,fmt,'EB');
      end;
      
      if efields.ifcaus
         if ~isfield(efields,'epnum')
            efields.epnum = bgs.elements(eid).epnum;
         end;
         fmt = buildfmt('%d8',efields.epnum);
         bgs.elements(eid).caus = read2partline(fid,fmt,'EC');
      end;
      
      if efields.ifphx
         bgs.elements(eid).phx = read2partline(fid,'%u8','EPHX');
      end;
      
      if efields.ifnn
         bgs.elements(eid).nn = read2partline(fid,'%d8','ENN',true);
      end;
      
      if efields.ifat
         bgs.elements(eid).at = read2partline(fid,'%n','EAT');
      end;
      
      if efields.ifmvars
         emvarsstr = read2partline(fid,'%q','EMVARS');
         bgs.elements(eid).mvars = readvec(emvarsstr,'%s',', ');
      end;
      
      if efields.ifexpr
         bgs.elements(eid).expr = readexpr(fid,bgs.elements(eid));
      end;
   end;
   
   read1partline(fid,'%s','EEND');
end;
read1partline(fid,'%s','ELEND');

% SECTION 11 - File End.
read1partline(fid,'%s','BGSDEND');

% cleanup function executes here
end

%% Input parser construction function

function ip = buildIP()
% The purpose of this function is to construct an input parser object for the
% readBGSD function.
%
% INPUTS/OUTPUTS
% ip - Input parser object for the readBGSD function.

ip = inputParser;
ip.StructExpand = false;
ip.addRequired('fname',@isexistfpath);
end

%% Cleanup function

function cleanup(fid)
% The purpose of this function is to run the cleanup routines when readBGSD
% exits, either because of an error, or regularly.
%
% INPUTS/OUTPUTS
% fid - BGSD file ID.

fclose(fid);
end

%% One-part line reading function

function val = read1partline(fid,fmt,def)
% The purpose of this function is to read a 1-part line of a specified format at
% the location of the cursor, proceed to the next line, and return the value
% that was read. Optionally, an expected value for the line can be specified,
% and an error is thrown if a different value is read.
%
% INPUTS/OUTPUTS
% fid - BGSD file ID.
% fmt - Format character sequence with which to read the line.
% def - Expected (default) value. If different from the read value, an error is
%       thrown. OPTIONAL.
% val - Read value.

C = textscan(fid,fmt,1,'CollectOutput',true);
if any(strcmp(fmt,{'%s','%q'}))
   val = C{1}{1};
else
   val = C{1};
end;

if nargin == 2
   % No expected value, do nothing
elseif nargin == 3
   if ischar(def)
      assert(strcmp(val,def),'Expected string %s, got %s instead!!!',def,val);
   elseif isnumeric(def)
      assert(val == def,'Expected value %g, got %g instead!!!',def,val);
   else
      error('Unknown expected data type. Cannot compare!!!');
   end;
else
   error('Invalid number of inputs!!!');
end;
end

%% Two-part line reading function

function val = read2partline(fid,fmt,def,ifbool)
% The purpose of this function is to read a 2-part line at the location of the
% cursor, and transition to the next line. The first entry is expected to be a
% keyword, and is read as a string. The second entry is of a specified data
% type. The second entry can be several values (like numbers), which the format
% string is supposed to reflect. Optionally, an expected value for the keyword
% can be specified, and an error is thrown if a different value is read. Another
% optional argument converts a read signed integer into a boolean.
%
% INPUTS/OUTPUTS
% fid - BGSD file ID.
% fmt - Format character sequence with which to read the line.
% def - Expected (default) keyword. If different from the read keyword, an error
%       is thrown. OPTIONAL.
% ifbool - Whether or not to convert the read integer value into a logical.
%          OPTIONAL. Default: false.
% val - Read value.

fullfmt = ['%s ' fmt];
C = textscan(fid,fullfmt,1,'CollectOutput',true);
kstr = C{1}{1};
if any(strcmp(fmt,{'%s','%q'}))
   val = C{1}{2};
else
   val = C{2};
end;

if nargin == 2
   % No expected keyword
   ifbool = false;
elseif nargin == 3
   assert(strcmp(kstr,def),'Expected keyword %s, got %s instead!!!',def,kstr);
   ifbool = false;
elseif nargin == 4
   assert(strcmp(kstr,def),'Expected keyword %s, got %s instead!!!',def,kstr);
   % ifbool specified in the input
else
   error('Invalid number of inputs!!!');
end;

if ifbool
   val = logical(val);
end;
end

%% One-part multiline list writing function

function list = read1partlist(fid,fmt,len)
% The purpose of this function is to read a 1-part multiline list of a specified
% format at the location of the cursor, proceed to the next line, and return the
% read list as an array. The list is returned either as a cell array of strings,
% or as an array of specified data type, depending on the format. Optionally,
% the list length can be specified.
%
% INPUTS/OUTPUTS
% fid - BGSD file ID.
% fmt - Format character sequence of a single line of the list.
% len - List length. OPTIONAL.
% list - Read list array.

if nargin == 2
   C = textscan(fid,fmt,'CollectOutput',true);
elseif nargin == 3
   C = textscan(fid,fmt,len,'CollectOutput',true);
else
   error('Invalid number of inputs!!!');
end;
list = C{1};
end

%% Matrix reading function

function mat = readmat(fid,fmt,ncols,nrows)
% The purpose of this function is to read a 2-dimensional array of a specified
% format at the location of the cursor, proceed to the next line, and return the
% array as either a cell array of strings, or as an array of specified data
% type, depending on the format. The number of columns has to be specified;
% optionally, the number of rows can be specified as well.
%
% INPUTS/OUTPUTS
% fid - BGSD file ID.
% fmt - Format character sequence of a single element of the array.
% ncols - Number of columns in the array. REQUIRED.
% nrows - Number of rows in the array. OPTIONAL.
% mat - Read 2-dimensional array.

fullfmt = repmat([fmt ' '],1,ncols-1);
fullfmt = [fullfmt fmt];

if nargin == 3
   C = textscan(fid,fullfmt,'CollectOutput',true);
elseif nargin == 4
   C = textscan(fid,fullfmt,nrows,'CollectOutput',true);
else
   error('Invalid number of inputs!!!');
end;
mat = C{1};
end

%% List-as-matrix reading function

function mat = readlistasmat(fid,fmt,ncols,nrows)
% The purpose of this function is to read a multiline list as a 2-dimensional
% array at the location of the cursor, proceed to the next line, and return the
% array as either a cell array of strings, or as an array of specified data
% type, depending on the format. The number of columns and rows in the resulting
% array has to be specified.
%
% INPUTS/OUTPUTS
% fid - BGSD file ID.
% fmt - Format character sequence of a single element of the array.
% ncols - Number of columns in the array. REQUIRED.
% nrows - Number of rows in the array. REQUIRED.
% mat - Read 2-dimensional array.

siz = ncols*nrows;
C = textscan(fid,fmt,siz,'CollectOutput',true);
mat = reshape(C{1},ncols,nrows).';
end

%% Element fields retrieval function

function efields = getefields(etype,exprtype)
% The purpose of this function is to take a given element type, and either state
% whether or not the element has an expression, or, if the expression type is
% given, to state what other fields the element has.
%
% INPUTS/OUTPUTS
% etype - Element type.
% exprtype - Expression type. OPTIONAL.
% efields - Element field status data structure.

if nargin == 1
   switch etype
      case {'SE','SF','I','C','R','TF','GY','MSE','MSF','MI','MC','MR','MTF',...
            'MGY','R2','MR2','RN','MRN'}
         efields.ifexpr = true;
         
      case {'1','0'}
         efields.ifexpr = false;
         
      otherwise
         error('Unknown element type!!!');
   end;
elseif nargin == 2
   switch etype
      case {'SE','SF','TF','GY'}
         efields.ifexpr = true;
         efields.ifx0 = false;
         efields.ifpnum = false;
         efields.ifbcon = false;
         efields.ifcaus = false;
         efields.ifphx = false;
         efields.ifnn = false;
         efields.ifat = false;
         efields.ifmvars = false;
         
      case {'I','C'}
         efields.ifexpr = true;
         efields.ifx0 = true;
         efields.ifpnum = false;
         efields.ifbcon = false;
         efields.ifcaus = false;
         efields.ifphx = true;
         efields.ifnn = true;
         efields.ifat = true;
         efields.ifmvars = false;
         
      case 'R'
         efields.ifexpr = true;
         efields.ifx0 = false;
         efields.ifpnum = false;
         efields.ifbcon = false;
         if strcmp(exprtype,'CC')
            efields.ifcaus = false;
         elseif any(strcmp(exprtype,{'CE','NE'}))
            efields.ifcaus = true;
            efields.epnum = uint16(1);
         elseif any(strcmp(exprtype,{'CMC','NMC','CME','NME'}))
            error('Modulated expression on an unmodulated element!!!');
         else
            error('Unknown expression type!!!');
         end;
         efields.ifphx = false;
         efields.ifnn = false;
         efields.ifat = false;
         efields.ifmvars = false;
         
      case {'1','0'}
         error('Expression type specified for a junction!!!');
         
      case {'MSE','MSF','MTF','MGY'}
         efields.ifexpr = true;
         efields.ifx0 = false;
         efields.ifpnum = false;
         efields.ifbcon = false;
         efields.ifcaus = false;
         efields.ifphx = false;
         efields.ifnn = false;
         efields.ifat = false;
         efields.ifmvars = true;
         
      case {'MI','MC'}
         efields.ifexpr = true;
         efields.ifx0 = true;
         efields.ifpnum = false;
         efields.ifbcon = false;
         efields.ifcaus = false;
         efields.ifphx = false;
         efields.ifnn = true;
         efields.ifat = true;
         efields.ifmvars = true;
         
      case 'MR'
         efields.ifexpr = true;
         efields.ifx0 = false;
         efields.ifpnum = false;
         efields.ifbcon = false;
         if any(strcmp(exprtype,{'CMC','NMC'}))
            efields.ifcaus = false;
         elseif any(strcmp(exprtype,{'CME','NME'}))
            efields.ifcaus = true;
            efields.epnum = uint16(1);
         elseif any(strcmp(exprtype,{'CC','CE','NE'}))
            error('Unmodulated expression on a modulated element!!!');
         else
            error('Uknown expression type!!!');
         end;
         efields.ifphx = false;
         efields.ifnn = false;
         efields.ifat = false;
         efields.ifmvars = true;
         
      case 'R2'
         efields.ifexpr = true;
         efields.ifx0 = false;
         efields.ifpnum = false;
         efields.ifbcon = false;
         efields.ifcaus = true;
         efields.epnum = uint16(2);
         efields.ifphx = false;
         efields.ifnn = false;
         efields.ifat = false;
         efields.ifmvars = false;
         
      case 'MR2'
         efields.ifexpr = true;
         efields.ifx0 = false;
         efields.ifpnum = false;
         efields.ifbcon = false;
         efields.ifcaus = true;
         efields.epnum = uint16(2);
         efields.ifphx = false;
         efields.ifnn = false;
         efields.ifat = false;
         efields.ifmvars = true;
         
      case 'RN'
         efields.ifexpr = true;
         efields.ifx0 = false;
         efields.ifpnum = true;
         efields.ifbcon = true;
         efields.ifcaus = true;
         efields.ifphx = false;
         efields.ifnn = false;
         efields.ifat = false;
         efields.ifmvars = false;
         
      case 'MRN'
         efields.ifexpr = true;
         efields.ifx0 = false;
         efields.ifpnum = true;
         efields.ifbcon = true;
         efields.ifcaus = true;
         efields.ifphx = false;
         efields.ifnn = false;
         efields.ifat = false;
         efields.ifmvars = true;
         
      otherwise
         error('Unknown element type!!!');
   end;
else
   error('Invalid number of inputs!!!');
end;
end

%% Vector reading function

function vec = readvec(str,fmt1,sep)
% The purpose of this function is to read a vector from a string, based on the
% format of a single element and an optional separator string.
%
% INPUTS/OUTPUTS
% str - String containing the vector.
% fmt1 - Format of a single element.
% sep - Separator character sequence. OPTIONAL. Default: ' ' (whitespace).
% vec - Vector read from the string. May be a cell vector of strings.

if nargin == 2
   sep = ' ';
elseif nargin == 3
   % sep specified in the input
else
   error('Invalid number of inputs!!!');
end;

C = textscan(str,fmt1,'WhiteSpace',sep);
vec = C{1};
end

%% Expression reading function

function expr = readexpr(fid,elem)
% The purpose of this function is to read the element expression based on all
% other known properties of the element.
%
% INPUTS/OUTPUTS
% fid - BGSD file ID.
% elem - Element data structure, excluding the element expression.
% expr - Element expression.

switch elem.etype
   case {'SE','SF','I','C','R','TF','GY','MSE','MSF','MI','MC','MR',...
         'MTF','MGY'}
      exprnum = uint16(1);
      
   case {'1','0'}
      error('Attempting to read an expression on a junction element!!!');
      
   case {'R2','MR2'}
      exprnum = uint16(2);
      
   case {'RN','MRN'}
      exprnum = elem.epnum;
      
   otherwise
      error('Unknown element type!!!');
end;

if exprnum == 1
   switch elem.exprtype
      case 'CC'
         expr = read1partline(fid,'%n');
         
      case {'CMC','CE','CME'}
         expr = read1partline(fid,'%q');
         
      case {'NMC','NE','NME'}
         expr = read1partline(fid,'%s');
         
      otherwise
         error('Unknown expression type!!!');
   end;
else
   switch elem.exprtype
      case 'CC'
         expr = readmat(fid,'%n',exprnum,exprnum);
         
      case 'CMC'
         expr = readlistasmat(fid,'%q',exprnum,exprnum);
         
      case 'NMC'
         expr = readlistasmat(fid,'%s',exprnum,exprnum);
         
      case 'CE'
         expr = readlistasmat(fid,'%q',1,exprnum);
         
      case 'NE'
         expr = readlistasmat(fid,'%s',1,exprnum);
         
      case 'CME'
         expr = readlistasmat(fid,'%q',1,exprnum);
         
      case 'NME'
         expr = readlistasmat(fid,'%s',1,exprnum);
         
      otherwise
         error('Unknown expression type!!!');
   end;
end;
end

%% Vector format string construction function

function fmt = buildfmt(fmt1,len,sep)
% The purpose of this function is to construct a format character sequence for
% reading a vector of values, based on a specified single-element format, the
% length of the vector and, optionally, the separator character(s).
%
% INPUTS/OUTPUTS
% fmt1 - Format string of a single element of the vector.
% len - Number of elements in the vector.
% sep - Separator string. OPTIONAL.
% fmt - Full format string of the vector.

if nargin == 2
   sep = ' ';
elseif nargin == 3
   % sep specified in the input
else
   error('Invalid number of inputs!!!');
end;

fmt = repmat([fmt1 sep],1,len-1);
fmt = [fmt fmt1]; 
end

%% (Parameter name, parameter value) pair reading function

function [pname,pval] = readpvalpair(fid,fmt)
% The purpose of this function is to read a (parameter name, parameter value)
% from the location of the cursor, according to the optionally specified
% value format. It then proceeds to the next line.
%
% INPUTS/OUTPUTS
% fid - BGSD file ID.
% fmtval - Format string of the parameter value. OPTIONAL. Default value: '%q'.
% pname - Parameter name read.
% pval - Parameter value read.

if nargin == 1
   fmt = '%q';
elseif nargin == 2
   % fmt specified in input
else
   error('Invalid number of inputs!!!');
end;

fullfmt = ['%s ' fmt];
C = textscan(fid,fullfmt,1,'CollectOutput',true);

pname = C{1}{1};
if any(strcmp(fmt,{'%s','%q'}))
   pval = C{1}{2};
else
   pval = C{2};
end;
end