%% fBGSolve
% Bond graph processing function.
%
% [T,X,[E,F]] = fBGSolve(fname,[runopts,bgs],['ovrrunopts',ovrrunopts])
%
% This function processes a bond graph system according to a specified set of
% run options. It takes a BGSD file name, and either reads it, or uses it as the
% root file name if the BGSD data structures are already provided as inputs. If
% the BGSD file exists, AND the BGSD data structures are specified, an error is
% thrown. Optionally, the function also accepts an overwriting runtime and
% recording options data structure, to adjust the run with.
% 
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       April 18, 2013
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Variable descriptions
%
% INPUTS
% BGSD file name.
% The BGSD file name specifies the BGSD file to process, or, if BGSD data
% structures are also input, to use as root file during the run. If the BGSD
% file name specified already exists, and the BGSD data structures are input, an
% error is thrown.
%
% fname - BGSD file name. This may be either the full file name (with path), or
%         only the file name itself. string. REQUIRED.
%
% BGSD data structures.
% The BGSD data structures are optional, but if they are given, both of them
% must be given. In that case, the specified file name will be used as the root
% file to process the system with. If the BGSD file name specified already
% exists, and the BGSD data structures are input, an error is thrown. The BGSD
% data structures get completed by completeBGSDdata prior to processing, so
% their specifications are the same as for completeBGSDdata. These
% specifications less stringent than those for writeBGSD (or that are returned
% by readBGSD).
%
% runopts - Run Options data structure. Describes sections 1-7 in the BGSD file.
%           All of the fields are optional (although with certain input values
%           in some fields, other fields can be made required). Its fields are
%           listed below.
%  SECTION 1 - File Header.
%  runopts.svc - Software version that created the data structures. double.
%  runopts.svf - Software version that the data structures are intended for.
%                double.
%  runopts.ifnotes - Whether or not a notes section describing the problem is
%                    present. Either boolean or int8.
%  runopts.notes - A short description of the file. string.
%  SECTION 2 - Runtime and Recording Options.
%  runopts.verbose - How verbose the processor should be during processing.
%                    int8.
%  runopts.rti - Whether or not to run the time integrator after the AEs are
%                sorted. Either boolean or int8.
%  runopts.log - How to log all the outputs specified by runopts.verbose. int8.
%  runopts.logfnspec - How the log file name is specified. int8.
%  runopts.logfn - Log file name. string.
%  runopts.smr - Whether or not to record the sorted expressions and supporting
%                variables into a .mat file after the algebraic equations are
%                sorted. Either boolean or int8.
%  runopts.smfnspec - How the sorted .mat file name is specified. int8.
%  runopts.smfn - Sorted .mat file name. string.
%  runopts.rtrr - Whether or not to record the integration results at runtime,
%                 as they are being computed. Either boolean or int8.
%  runopts.rtrrfqspec - How the runtime results recording frequency is
%                       specified. int8.
%  runopts.rtrrfq - Runtime results recording frequency. uint32.
%  runopts.rtrrfnspec - How the runtime results file name is specified. int8.
%  runopts.rtrrfn - Runtime results file name. string.
%  runopts.rrm - Whether or not to record the final results in a .mat file.
%                Either boolean or int8.
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
%                       a string (gets evaluated prior to processing), numeric
%                       scalar, or a function handle.
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
%                      prior to processing), numeric scalar or a function
%                      handle.
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
%                       a string (gets evaluated prior to processing), numeric
%                       scalar or a function handle.
%  SECTION 5 - Evaluator Options.
%  runopts.cs2n - Whether or not to convert the symbolic expressions to numeric
%                 expressions after sorting. Either boolean or int8.
%  runopts.nft - Numeric function type to convert the symbolic expressions to.
%                int8.
%  SECTION 6 - Post-Processor Options.
%  runopts.epp - Number of effort variables to evaluate during post-processing.
%                int32.
%  runopts.eppIDs - List of BIDs of effort variables to evaluate during
%                   post-processing. Either uint16 vector (BIDs of effort
%                   variables) or boolean vector (each element corresponds to a
%                   BID, true to evaluate the effort).
%  runopts.fpp - Number of flow variables to evaluate during post-processing.
%                int32.
%  runopts.fppIDs - List of BIDs of flow variables to evaluate during
%                   post-processing. Either uint16 vector (BIDs of flow
%                   variables) or boolean vector (each element corresponds to a
%                   BID, true to evaluate the flow).
%  runopts.ppdisp - Post-processing display level. int8.
%  SECTION 7 - Linked Files and Directories.
%  runopts.lfnum - Number of linked files. uint8.
%  runopts.lflist - List of linked file names. Cell vector of strings.
%  runopts.ldnum - Number of linked directories. uint8.
%  runopts.ldlist - List of linked directory names. Cell vector of strings.
%
% bgs - Bond Graph System data structure. Describes sections 8-10 in the BGSD
%       file. Its fields are listed below.
%  SECTION 8 - Additional System Information.
%  bgs.asi.sci - Specified causality information. string. OPTIONAL.
%  bgs.asi.smi - Specified modulation information. string. OPTIONAL.
%  bgs.asi.sli - Specified linearity information. string. OPTIONAL.
%  bgs.asi.sri - Specified reduction information. string. OPTIONAL.
%  bgs.asi.sei - Specified error information. string. OPTIONAL.
%  SECTION 9 - Bond Connectivity Map.
%  bgs.bnum - Number of bonds in the system. uint16. OPTIONAL.
%  bgs.bonds - Bond directionality table. Either data structure with fields
%              .from and .to (uint16 vectors, one element per bond), or uint16
%              [bnum x 2] array, with entries (to,from) in each row. REQUIRED.
%  SECTION 10 - Element List.
%  The fields in this section can be input in one of two ways:
%   a). As vectors of element types, expressions, etc.
%   b). As a vector of element data structures, with each element of the vector
%       containing all of the corresponding element's information.
%  Both ways are described below.
%  Element list storage type (a).
%   bgs.enum - Number of elements in the system. uint16. OPTIONAL.
%   bgs.exprnum - Number of elements with expressions. uint16. OPTIONAL.
%   bgs.emodnum - Number of modulated elements. uint16. OPTIONAL.
%   bgs.epambnum - Number of elements with ambiguous numbers of ports. uint16.
%                  OPTIONAL.
%   bgs.ebambnum - Number of elements with ambiguous bond-to-port
%                  connectivities. uint16. OPTIONAL.
%   bgs.ecambnum - Number of elements with ambiguous causalities on their ports.
%                  uint16. OPTIONAL.
%   bgs.ennum - Number of elements with numeric expressions. uint16. OPTIONAL.
%   bgs.nvnum - Number of numeric variables. uint16. OPTIONAL.
%   bgs.xnum - Number of storage elements with initial values. uint16. OPTIONAL.
%   bgs.capnum - Number of capacitive elements. uint16. OPTIONAL.
%   bgs.inertnum - Number of inertial elements. uint16. OPTIONAL.
%   bgs.phxspecnum - Number of physics-specific elements. uint16. OPTIONAL.
%   bgs.nnspecnum - Number of non-negativity-specific elements. uint16.
%                   OPTIONAL.
%   bgs.atspecnum - Number of elements with specified absolute tolerances.
%                   uint16. OPTIONAL.
%   bgs.etypes - List of element types, in order of increasing EIDs. Cell vector
%                of strings. REQUIRED.
%   bgs.exprIDs - List of EIDs with expressions, in order of increasing EIDs.
%                 Either uint16 vector of EIDs, one EID with expression per
%                 element, or boolean vector with entries corresponding to all
%                 EIDs (true if EID has an expression). OPTIONAL.
%   bgs.exprtypes - List of expresion types. Cell vector of strings, one vector
%                   element per expression. REQUIRED.
%   bgs.exprs - List of elements' expressions. Cell vector of various objects,
%               one element's expression per cell vector element. The individual
%               expressions may be strings (for all expression types but CCs),
%               doubles (for CCs) or symbolic expressions (for CMCs, CEs, CMEs).
%               For multiport elements with CC expressions the coefficients are
%               stored as double matrices. For multiport elements with non-CC
%               expressions the expressions are stored as either cell arrays
%               (for strings), or symbolic arrays. CMC and NMC expression arrays
%               can be entered either as matrices, or as vectors (ordering
%               coefficients left to right, then top to bottom). They are
%               returned as matrices. All individual expressions of an element
%               have to be of the same class. REQUIRED.
%   bgs.modIDs - List of EIDs with modulated expressions, in order of increasing
%                EIDs. Either uint16 vector of EIDs, one EID with modulated
%                expression per vector element, or boolean vector with entries
%                corresponding to all EIDs (true if EID has a modulated
%                expression). OPTIONAL.
%   bgs.mvars - List of elements' modulating variable lists. Cell vector of
%               either cell vectors of strings, or symbolic vectors, one
%               modulating variable list per (symbolic or string cell) vector
%               element. REQUIRED.
%   bgs.pambIDs - List of EIDs with ambiguous port numbers, in order of
%                 increasing EIDs. Either uint16 vector of EIDs, one EID with
%                 ambiguous port number per vector element, or boolean vector
%                 with entries corresponding to all EIDs (true if EID has
%                 ambiguous port number). OPTIONAL.
%   bgs.epnums - List of port numbers of the elements with ambiguous numbers of
%                ports. uint16 vector, one number of ports per vector element.
%                OPTIONAL.
%   bgs.bambIDs - List of EIDs with ambiguous bond-to-port connectivities, in
%                 order of increasing EIDs. Either uint16 vector of EIDs, one
%                 EID with ambiguous bond-to-port connectivity per vector
%                 element, or boolean vector with entries corresponding to all
%                 EIDs (true if EID has ambiguous bond-to-port connectivity).
%                 OPTIONAL.
%   bgs.ebcons - List of bond-to-port connectivities of the elements with
%                ambiguous bond-to-port connectivities. Cell vector of uint16
%                vectors, one bond-to-port connectivity (uint16 vector of BIDs)
%                per cell vector element. REQUIRED.
%   bgs.cambIDs - List of EIDs with ambiguous causalities, in order of
%                 increasing EIDs. Either uint16 vector of EIDs, one EID with
%                 ambiguous causality per vector element, or boolean vector with
%                 entries corresponding to all EIDs (true if EID has ambiguous
%                 port number). OPTIONAL.
%   bgs.causes - List of causalities of elements with ambiguous causalities.
%                Cell vector of int8 vectors, one element's causality (int8
%                vector of element ports' input types) per cell vector element.
%                REQUIRED.
%   bgs.numIDs - List of EIDs with numeric expressions, in order of increasing
%                EIDs. Either uint16 vector of EIDs, one EID with numeric
%                expression per vector element, or boolean vector with entries
%                corresponding to all EIDs (true if EID has numeric expression).
%                OPTIONAL.
%   bgs.x0IDs - List of EIDs of storage elements with initial values. Either
%               uint16 vector, one storage EID per vector element, or boolean
%               vector with entries corresponding to all EIDs (true if EID is a
%               storage element). OPTIONAL.
%   bgs.capIDs - List of EIDs of capacitive elements with initial values. Either
%                uint16 vector, one capacitive EID per vector element, or
%                boolean vector with entries corresponding to all EIDs (true if
%                EID is a capacitive element). OPTIONAL.
%   bgs.inertIDs - List of EIDs of inertial elements with initial values. Either
%                  uint16 vector, one inertial EID per vector element, or
%                  boolean vector with entries corresponding to all EIDs (true
%                  if EID is an inertial element). OPTIONAL.
%   bgs.x0s - List of initial values of storage elements. double vector, one
%             entry per storage element. REQUIRED.
%   bgs.phxIDs - List of EIDs of physics-specific elements. Either uint16
%                vector, one physics-specific EID per vector element, or boolean
%                vector with entries corresponding to all EIDs (true if EID is
%                a physics-specific element). OPTIONAL.
%   bgs.phxs - List of physics-specific elements' physics. uint8 vector, one
%              entry per physics-specific element.
%   bgs.nnIDs - List of EIDs of non-negativity-specific elements'. Either uint16
%               vector, one non-negativity-specific EID per vector element, or
%               boolean vector with entries corresponding to all EIDs (true if
%               EID is a non-negativity-specific element). OPTIONAL.
%   bgs.nns - List of non-negativity specifications of non-negativity-specific
%             elements. Boolean vector, one entry per non-negativity-specific
%             element.
%   bgs.atIDs - List of EIDs of elements with specified absolute tolerances.
%               Either uint16 vector, one absolute tolerance-specific EID per
%               vector element, or boolean vector with entries corresponding to
%               all EIDs (true if EID is an absolute tolerance-specific
%               element). OPTIONAL.
%   bgs.ats - List of absolute tolerances of absolute tolerance-specific
%             elements. double vector, one entry per absolute tolerance-specific
%             element.
%  Element list storage type (b).
%   bgs.enum - Number of elements in the system. uint16. OPTIONAL.
%   bgs.exprnum - Number of elements with expressions. uint16. OPTIONAL.
%   bgs.emodnum - Number of modulated elements. uint16. OPTIONAL.
%   bgs.epambnum - Number of elements with ambiguous numbers of ports. uint16.
%                  OPTIONAL.
%   bgs.ebambnum - Number of elements with ambiguous bond-to-port
%                  connectivities. uint16. OPTIONAL.
%   bgs.ecambnum - Number of elements with ambiguous causalities on their ports.
%                  uint16. OPTIONAL.
%   bgs.ennum - Number of elements with numeric expressions. uint16. OPTIONAL.
%   bgs.nvnum - Number of numeric variables. uint16. OPTIONAL.
%   bgs.xnum - Number of storage elements with initial values. uint16. OPTIONAL.
%   bgs.capnum - Number of capacitive elements. uint16. OPTIONAL.
%   bgs.inertnum - Number of inertial elements. uint16. OPTIONAL.
%   bgs.phxspecnum - Number of physics-specific elements. uint16. OPTIONAL.
%   bgs.nnspecnum - Number of non-negativity-specific elements. uint16.
%                   OPTIONAL.
%   bgs.atspecnum - Number of elements with specified absolute tolerances.
%                   uint16. OPTIONAL.
%   bgs.elements - List of element objects (data structures), in order of
%                  increasing EIDs. struct vector, one vector element per EID.
%                  REQUIRED. An element data structure has the following fields
%                  (some of them may be empty):
%    .etype - Element type. string. REQUIRED.
%    .exprtype - Element's expression type. string. REQUIRED if the element has
%                an expression.
%    .x0 - Element's initial value. double. REQUIRED if the element is a storage
%          element.
%    .epnum - Element's port number. uint16. OPTIONAL.
%    .ebcon - Element's bond-to-port connectivity. uint16 vector, one BID per
%             vector element. REQUIRED if the element has ambiguous bond-to-port
%             connectivity.
%    .caus - Element's causality. int8 vector, one input type per vector
%            element. REQUIRED if the element's causality is ambiguous.
%    .phx - Element's physics. uint8. OPTIONAL.
%    .nn - Element's non-negativity status. boolean. OPTIONAL.
%    .at - Element's absolute tolerance. double. OPTIONAL.
%    .mvars - Element's modulating variables. Either a cell vector of strings or
%             a symbolic vector, one modulating variable per (symbolic or string
%             cell) vector element. REQUIRED for modulated elements.
%    .expr - Element's expression(s). Individual expressions may be strings (for
%            all expression types but CCs), doubles (for CCs) or symbolic
%            expressions (for CMCs, CEs, CMEs). For multiport elements with CC
%            expressions the coefficients are stored as double matrices. For
%            multiport elements with non-CC expressions the expressions are
%            stored as either cell arrays (for strings), or symbolic arrays. CMC
%            and NMC expression arrays can be entered either as matrices, or as
%            vectors (ordering coefficients left to right, then top to bottom).
%            All individual expressions of an element have to be of the same
%            class. REQUIRED if the element has an expression.
%
% OPTIONAL INPUTS
% ('ovrrunopts',ovrrunopts)
%  Runtime and Recording Options data structure whose fields (all of them
%  optional) overwrite the runopts in the BGSD file or supplied to the function.
%  All fields in ovrrunopts are optional; after its fields overwrite the
%  corresponding fields in runopts, the data structures are completed and then
%  processed.
%  It is recommended to use ovrrunopts, primarily, to turn on silent processing,
%  logging, and to modify integrator specs for comparison studies.
%
% OUTPUTS
% T - Vertical vector of time points. double.
% X - 2-dimensional array of horizontal state vectors, one horizontal state
%     vector per time point. double.
% E - 2-dimensional array of horizontal effort vectors, one horizontal vector
%     per time point. double.
% F - 2-dimensional array of horizontal flow vectors, one horizontal vector per
%     time point. double.
%
% STANDARD VARIABLES
% eid - Element ID; used to refer to elementwise vectors.
% bid - Bond ID; used to refer to bondwise vectors.
% pid - Port ID; used to refer to elements' port vectors.
%
% PROCESSING VARIABLES
% refs - Referencing arrays data structure. Contains index arrays that are used
%        to reference quantities by EID (or BID), and not by quantity ID.
%        refs.Quantity(EID) is the quantity ID for quantity Quantity of element
%        EID. Referenced quantities are:
%  - Expressions (for elements with expressions).
%  - Expression types (for elements with expressions).
%  - Modulating variable lists (for modulated elements).
%  - Port numbers (for elements with ambiguous numbers of ports).
%  - Bond-to-port connectivities (for elements with ambiguous bond-to-port
%    connectivities).
%  - Causalities (for elements with ambiguous causalities).
%  - State variables' XIDs (state vector indices).
%  - Effort variables' BVIDs.
%  - Flow varaibles' BVIDs.
%        Its fields are:
%  .exprs - Expression referencing array.
%  .mvars - Modulating variable list referencing array.
%  .pnums - Port number referencing array.
%  .bcons - Bond-to-port connectivity referencing array.
%  .causes - Causality referencing array.
%  .nums - Numeric element referencing array.
%  .xids - State variables' state vector index array.
%  .ebvid - Effort BVID calculation function (of bid).
%  .fbvid - Flow BVID calculation function (of bid).
%  .bidebv - BID based on Effort BVID calculation function.
%  .bidfbv - BID based on Flow BVID calculation function.
% vars - Symbolic variable data structure. Contains symbolic variables that
%        are used for sorting and formulation of algebraic equations. Its
%        fields are:
%  .ui.t - Time (unindexed).
%  .ui.e - Effort (unindexed).
%  .ui.f - Flow (unindexed).
%  .ui.q - Displacement (unindexed).
%  .ui.p - Momentum (unindexed).
%  .ui.ei - Input effort (unindexed).
%  .ui.eo - Output effort (unindexed).
%  .ui.fi - Input flow (unindexed).
%  .ui.fo - Output flow (unindexed).
%  .e - Efforts vector.
%  .f - Flows vector.
%  .b - Bond variables vector (all efforts, all flows).
%  .q - Displacements vector.
%  .p - Momenta vector.
%  .ep - Port efforts vector.
%  .fp - Port flows vector.
%  .x - State vector (indexed in order of increasing EIDs).
%  .n - Numeric variables vector.
%
% ebcm - Element-indexed Bond Connectivity Map data structure. Its fields
%        are:
%  .to{EID} - Horizontal vector of BIDs of bonds pointing to element EID.
%             uint16 vector.
%  .from{EID} - Horizontal vector of BIDs of bonds pointing to element EID.
%               uint16 vector.
%  .conn{EID} - Horizontal vector of BIDs of bonds connected to element
%               EID. uint16 vector.

%% Master BGS processing function

function [T,X,E,F] = fBGSolve(varargin)
% The purpose of this function is to process a bond graph system according to a
% specified set of run options. It takes a BGSD file name, and either reads it,
% or uses it as the root file name if the BGSD data structures are already
% provided as inputs. If the BGSD file exists, AND the BGSD data structures are
% specified, an error is thrown. Optionally, the function also accepts an
% overwriting runtime and recording options data structure, to adjust the run
% with.
%
% INPUTS/OUTPUTS
% fname - BGSD file name.
% runopts - Runtime and Recording Options data structure.
% bgs - Bond Graph System data structure. Element list storage type b is
%       expected.
% ovrrunopts - Runtime and Recording Options data structure to overwrite runopts
%              with.
% T - Time point vector.
% X - State solution array.
% E - Effort solution array.
% F - Flow solution array.
% See the variable description cell for more information.

% Input parser construction
ip = buildIP();
% Input parsing
ip.parse(varargin{:});
% Retrieving inputs
fname = ip.Results.fname;
if isexistfpath(fname)
   assert(isempty(ip.Results.runopts) && isempty(ip.Results.bgs),...
      'Cannot use an existing BGSD file name as root file!!!');
   [runopts,bgs] = readBGSD(fname);
else
   assert(~(isempty(ip.Results.runopts) || isempty(ip.Results.bgs)),...
      'No existing file or data structures input!!!');
   runopts = ip.Results.runopts;
   bgs = ip.Results.bgs;
end;
ovrrunopts = ip.Results.ovrrunopts;
% Clearing memory
clear('ip');

% Switching directory to the BGSD file root
[fpath,sfname,~] = fileparts(fname);
if isempty(fpath)
   oldpath = pwd;
   rootpath = pwd;
else
   oldpath = cd(fpath);
   rootpath = pwd;
end;

% Constructing temporary cleanup routine 1
cleanupObj1 = onCleanup(@()cleanup(1,'oldpath',oldpath));
% Clearing memory
clear('oldpath');

% Overwriting run options (if necessary)
if ~isempty(ovrrunopts)
   runopts = dsfwovrwrite(runopts,ovrrunopts);
end;
% Clearing memory
clear('ovrrunopts');

% Completing and converting data structures
[runopts,bgs] = completeBGSDdata(runopts,bgs,'opmode','BGSolve');

% Checking software versions
assert(runopts.svc == 1.03 && runopts.svf == 1.03,'Wrong software version!!!');

% Expanding runopts
runopts = expandrunopts(runopts,sfname);

% Creating required subdirectories
if ~exist('Results','dir')
   mkdir('Results');
end;
if ~exist('TmpFunctions','dir')
   mkdir('TmpFunctions');
end;
% Linking required subdirectories
addpath('Results','TmpFunctions','-begin');

% Creating required files
fids = createfiles(runopts);

log2spec('BGSolver v1.03 by Eugeny Sosnovsky.',1,fids,runopts);
log2spec('Beginning bond graph system processing...',1,fids,runopts);
log2spec('Inputs have been accepted.',1,fids,runopts);

log2spec('Data structures have been completed.',2,fids,runopts);
log2spec('Results and TmpFunctions subdirectories have been added.',...
   2,fids,runopts);
log2spec('Required files have been created (if any).',2,fids,runopts);
log2spec(['Cleanup routine 1 (switch to original directory) ' ...
   'has been constructed.'],2,fids,runopts);

% Constructing temporary cleanup routine 2
cleanupObj2 = onCleanup(@()cleanup(2,'fids',fids));

log2spec(['Cleanup routine 2 (closing of all open files) ' ...
   'has been constructed.'],2,fids,runopts);

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

% Checking specified system information
assert(strcmp(bgs.asi.sci,'FC'),...
   'Can only process fully causal systems => conflict with ASI.SCI!!!');
assert(strcmp(bgs.asi.sei,'DNCE'),...
   'Can only process error-free systems => conflict with ASI.SEI!!!');

%
% STUB
% Other error checking procedures can be inserted here.
% /STUB
%

log2spec('---',1,fids,runopts);
log2spec('Error checking procedures complete.',1,fids,runopts);

% Constructing referencing arrays
refs = createrefs(bgs);

log2spec('Referencing arrays created.',1,fids,runopts);

% Instantiating symbolic variables
vars = createvars(bgs,refs);

log2spec('Symbolic variables instantiated.',1,fids,runopts);

% Constructing element-indexed bond connectivity map
ebcm = createebcm(bgs);

log2spec('Element-indexed bond connectivity map created.',1,fids,runopts);
log2spec('---',1,fids,runopts);
log2spec('Converting numeric functions to function handles...',1,fids,runopts);

% Converting numeric functions to function handles
for nEID = bgs.numIDs
   if any(strcmp(bgs.etypes{nEID},{'R2','MR2','RN','MRN'}))
      if iscellstr(bgs.exprs{refs.exprs(nEID)})
         bgs.exprs{refs.exprs(nEID)} = ...
            cellfun(@eval,bgs.exprs{refs.exprs(nEID)},'UniformOutput',false);
      end;
   else
      if ischar(bgs.exprs{refs.exprs(nEID)})
         bgs.exprs{refs.exprs(nEID)} = eval(bgs.exprs{refs.exprs(nEID)});
      end;
   end;
end;

log2spec('Numeric function conversion complete.',1,fids,runopts);
log2spec('Converting additional time integrator options...',1,fids,runopts);

% Converting additional time integrator, solver and jumpstart options
for addtioID = 1:runopts.addtio
   runopts.addtiolist{addtioID,2} = eval(runopts.addtiolist{addtioID,2});
end;
for addsoID = 1:runopts.addso
   runopts.addsolist{addsoID,2} = eval(runopts.addsolist{addsoID,2});
end;
for addjsoID = 1:runopts.addjso
   runopts.addjsolist{addjsoID,2} = eval(runopts.addjsolist{addjsoID,2});
end;

log2spec('Additional TI, solver and jumpstart options conversion complete.',...
   1,fids,runopts);
log2spec('---',1,fids,runopts);
log2spec('BOND GRAPH PROCESS STEP 4.',1,fids,runopts);
log2spec('Formulating algebraic equations...',1,fids,runopts);

% BOND GRAPH PROCESS STEP 4
% Formulating algebraic equations
srtsys = formulateAEs(bgs,refs,vars,ebcm);
% Clearing memory
clear('bgs','ebcm');

log2spec('Algebraic equations formulation complete.',1,fids,runopts);

% Preparing warning controls for sorting
[msgstr,msgid] = lastwarn('','');
warnstate = warning('off','symbolic:solve:warnmsg3');

log2spec('Warning controls prepared for sorting.',2,fids,runopts);

% Constructing temporary cleanup routine 4
cleanupObj4 = onCleanup(@()cleanup(4,'msgstr',msgstr,'msgid',msgid,...
   'warnstate',warnstate));
% Clearing memory
clear('msgstr','msgid','warnstate');

log2spec(['Cleanup routine 4 (warning control state initialization) has ' ...
   'been constructed.'],2,fids,runopts);
log2spec('---',1,fids,runopts);
log2spec('BOND GRAPH PROCESS STEP 5.',1,fids,runopts);
log2spec('Sorting algebraic equations...',1,fids,runopts);

% BOND GRAPH PROCESS STEP 5
% Reordering state derivative vector
[srtsys,vars,refs] = reorderXasBVIDs(srtsys,vars,refs);
srtsys.xnns = uint16(find(srtsys.xnns));
log2spec('State derivative vector reordered in order of increasing BVID.',2,...
   fids,runopts);
% Sorting algebraic equations
[rASbvids,rNSbvids] = calcIDsForASBVs(srtsys,refs);
log2spec('BV reordering vectors constructed.',2,fids,runopts);
bc = num2cell(vars.b(rNSbvids));
srtsys.b = cell(srtsys.bvnum,1);
log2spec('Sorting arrays preallocated.',2,fids,runopts);
[srtsys.b{:}] = solve(srtsys.aes,bc{:});
assert(strcmp(lastwarn,''),'BG Process Step 5 - Sorting - failed!!!');
log2spec('Algebraic equations successfully sorted.',1,fids,runopts);
srtsys.b = sym(srtsys.b(rASbvids));
log2spec('Bond variable vector formulated.',1,fids,runopts);
% Identifying variable dependencies in the sorted bond variable vector
[srtsys.bvitids,srtsys.bvinpxids,srtsys.bvinpnids] = getevids(srtsys.b,vars,...
   refs,'vars2check','txn');
log2spec('Bond variable vector variable dependencies identified.',2,fids,...
   runopts);
% Assigning numeric layers
srtsys = assignlayers(srtsys);
log2spec('Numeric layers assigned.',2,fids,runopts);
% Identifying BVs and NVs required for xdot evaluation and post-processing
srtsys = getvars2evalxdot(srtsys);
log2spec('BVs and NVs required for xdot identified.',2,fids,runopts);
srtsys = getvars2evalpp(srtsys,runopts,refs);
log2spec('BVs and NVs required for post-processing identified.',2,fids,runopts);
% Clearing memory
clear('rASbvids','rNSbvids');

log2spec('---',1,fids,runopts);
log2spec('BOND GRAPH PROCESS STEP 6.',1,fids,runopts);
log2spec('Formulating state derivative vector...',1,fids,runopts);

% BOND GRAPH PROCESS STEP 6
% Formulating state derivative vector for time integrator
xdotsys = createxdotsys(srtsys,vars,runopts);
xdot = @(t,x) evalxdot(t,x,xdotsys);
% Formulating Jacobian sparsity pattern
srtsys = getspatternxdot(srtsys);
log2spec('State derivative Jacobian sparsity pattern formulated.',...
   1,fids,runopts);
% Creating time integrator options
tioptions = createtioptions(srtsys,runopts,fids,'Vectorized','on');
log2spec('State derivative vector formulated.',1,fids,runopts);
log2spec('Time integrator options formulated.',1,fids,runopts);
log2spec('Formulating post-processing vector...',1,fids,runopts);
% Formulating post-processing vector
ppsys = createppsys(srtsys,vars,runopts);
log2spec('Post-processing vector formulated.',1,fids,runopts);
X0 = srtsys.x0s;
% Extracting reverse ordering X vector
revxinds = srtsys.revxinds;
log2spec('State vector reverse reordering vector extracted.',2,fids,runopts);

% Recording sorted expressions and variables
if runopts.smr
   fname = GetFullPath(fname); %#ok<NASGU>
   save(runopts.smfn,'fname','runopts','srtsys','vars','refs','xdotsys',...
      'xdot','tioptions','ppsys');
   
   log2spec('Sorted expressions and variables recorded in .mat file:',...
      1,fids,runopts);
   log2spec(['--> ' runopts.smfn],1,fids,runopts);
end;

% Checking Run Time Integrator flag
if ~runopts.rti
   T = [];
   X = [];
   E = [];
   F = [];
   
   log2spec('---',1,fids,runopts);
   log2spec('Processing stopped before Bond Graph Process Step 7.',1,fids,...
      runopts);
   log2spec('Returned [T,X,E,F] arrays are empty dummy variables.',2,fids,...
      runopts);
   log2spec('It''s all over. Noone will ever love you. I am quitting...',...
      1,fids,runopts);
   
   return;
end;

% Run Time Integrator flag is on
% Clearing memory
clear('srtsys','vars','refs','fname');

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

log2spec('Bond graph system processing complete.',1,fids,runopts);
log2spec('It''s all over. Noone will ever love you. I am quitting...',...
   1,fids,runopts);

% Cleanup routines 1-4 are executed here, in arbitrary order, but not
% simultaneously.
end

%% Input parser construction function

function ip = buildIP()
% The purpose of this function is to construct an input parser object for the
% fBGSolve function.
%
% INPUTS/OUTPUTS
% ip - Input parser object for the fBGSolve function.

fchkfun = @(fname) islegalfpath(fname,'ext','.bgsd');

ip = inputParser;
ip.StructExpand = false;
ip.addRequired('fname',fchkfun);
ip.addOptional('runopts',[],@isstruct);
ip.addOptional('bgs',[],@isstruct);
ip.addParamValue('ovrrunopts',[],@isstruct);
end

%% Cleanup function

function cleanup(varargin)
% The purpose of this function is to run the cleanup routines when fBGSolve
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
%  4 - Cleanup routine 4: restores warning control state.
% oldpath - Original directory that fBGSolve is called from. REQUIRED.
% fids - File IDs data structure. OPTIONAL.
% runopts - Runtime and Recording options data structure. OPTIONAL.
% rootpath - Absolute root folder path. OPTIONAL.
% msgstr - Last warning message string. OPTIONAL.
% msgid - Last warning message identifier. OPTIONAL.
% warnstate - Explicit symbolic solution warning state. OPTIONAL.

% Constructing input parser
ip = inputParser;
ip.StructExpand = false;
fopt = @(opt) any(opt == [1,2,3,4]);
foldpath = @(oldpath) exist(oldpath,'dir');
fwarnstate = @(warnstate) isfield(warnstate,'state');
ip.addRequired('opt',fopt);
ip.addParamValue('oldpath',pwd,foldpath);
ip.addParamValue('fids',struct(),@isstruct);
ip.addParamValue('runopts',struct(),@isstruct);
ip.addParamValue('rootpath',pwd,@ischar);
ip.addParamValue('msgstr','',@ischar);
ip.addParamValue('msgid','',@ischar);
ip.addParamValue('warnstate',struct('state','on'),fwarnstate);
% Parsing input
ip.parse(varargin{:});
% Retrieving inputs
opt = ip.Results.opt;
oldpath = ip.Results.oldpath;
fids = ip.Results.fids;
runopts = ip.Results.runopts;
rootpath = ip.Results.rootpath;
msgstr = ip.Results.msgstr;
msgid = ip.Results.msgid;
warnstate = ip.Results.warnstate;
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
elseif opt == 4
   % Cleanup routine 4: restoring warning control state.
   lastwarn(msgstr,msgid);
   warning(warnstate);
end;
end

%% Runopts expansion function

function runopts = expandrunopts(runopts,sfname)
% The purpose of this function is to replace the default specifications with the
% corresponding settings and file names, to obtain all user input required, and
% to compute the fixed time step size, if required.
%
% INPUTS/OUTPUTS
% runopts - Runtime and Recording Options data structure.
% sfname - BGSD file name.

% Instantiating user input settings data structure
inpset.logfn = false;
inpset.smfn = false;
inpset.rtrrfq = false;
inpset.rtrrfn = false;
inpset.rrmfn = false;
inpset.tit = false;
inpset.ii = false;
inpset.fts = false;

% Replacing default specs with corresponding settings and file names
if runopts.verbose ~= 0 && runopts.log ~= 0
   if runopts.logfnspec == 0
      runopts.logfn = ['Results' filesep sfname '.log'];
   elseif runopts.logfnspec == 2
      inpset.logfn = true;
   end;
end;

if runopts.smr
   if runopts.smfnspec == 0
      runopts.smfn = ['Results' filesep sfname '_SrtRes.mat'];
   elseif runopts.smfnspec == 2
      inpset.smfn = true;
   end;
end;

if runopts.rtrr
   if runopts.rtrrfqspec == 2
      inpset.rtrrfq = true;
   end;
   if runopts.rtrrfnspec == 0
      runopts.rtrrfn = ['Results' filesep sfname '_RTRes.txt'];
   elseif runopts.rtrrfnspec == 2
      inpset.rtrrfn = true;
   end;
end;

if runopts.rrm
   if runopts.rrmfnspec == 0
      runopts.rrmfn = ['Results' filesep sfname '_FinRes.mat'];
   elseif runopts.rrmfnspec == 2
      inpset.rrmfn = true;
   end;
end;

if runopts.titspec == 2
   inpset.tit = true;
end;

if runopts.iispec == 1
   inpset.ii = true;
end;

if any(runopts.tst == [1,2])
   if runopts.ftsspec == 1
      inpset.fts = true;
   end;
end;

for lfid = 1:runopts.lfnum
   if exist(runopts.lflist{lfid},'file') == 2
      % The file is in root or in full path, do nothing
   elseif exist(['LinkedResources' filesep runopts.lflist{lfid}],'file') == 2
      runopts.lflist{lfid} = ['LinkedResources' filesep runopts.lflist{lfid}];
   else
      error('Unable to find file %d in linked files list!!!',lfid);
   end;
end;

for ldid = 1:runopts.ldnum
   if exist(runopts.ldlist{ldid},'dir')
      % The folder is in root or in full path, do nothing
   elseif exist(['LinkedResources' filesep runopts.ldlist{ldid}],'dir')
      runopts.ldlist{ldid} = ['LinkedResources' filesep runopts.ldlist{ldid}];
   else
      error('Unable to find directory %d in linked directories list!!!',ldid);
   end;
end;

% Checking if input requirements violate the verbose requirements
assert(runopts.verbose ~= 0 || ~any(structfun(@(f) f,inpset)),...
   'User input required during silent processing!!!');

% Getting required user input
runopts = getusrinput(runopts,inpset);

% Computing fixed time step size, if necessary
if any(runopts.tst == [1,2])
   if runopts.ftsspect == 0
      runopts.fts = diff(runopts.ii) / double(runopts.fts);
   elseif runopts.ftsspect == 1
      % Fixed time step size already specified
   end;
end;
end

%% User input collection function

function runopts = getusrinput(runopts,inpset)
% The purpose of this function is to get all necessary inputs from the user that
% the inpset data structure specifies. If the user clicks cancel on any of the
% dialog windows, an error will be thrown.
%
% INPUTS/OUTPUTS
% runopts - Runtime and Recording Options data structure.
% inpset - Input Settings data structure.

% File name input
% runopts.logfn
if inpset.logfn
   runopts.logfn = savefiledlg({'*.log','Log file (*.log)'},...
      'Choose the log file to log the runtime in');
end;

% runopts.smfn
if inpset.smfn
   runopts.smfn = savefiledlg({'*.mat','Sorted expressions file (*.mat)'},...
      'Choose the .mat file to save sorted expressions in');
end;

% runopts.rtrrfn
if inpset.rtrrfn
   runopts.rtrrfn = savefiledlg({'*.txt','Runtime results file (*.txt)'},...
      'Choose the .txt file to save runtime results in');
end;

% runopts.rrmfn
if inpset.rrmfn
   runopts.rrmfn = savefiledlg({'*.mat','Results file (*.mat)'},...
      'Choose the .mat file to save results in');
end;

% Parameter input
dlginps = {'rtrrfq','tit','ii','fts'};
inpstats = cellfun(@(str)inpset.(str),dlginps);
inpnums = [1 1 2 1];
inpnum = sum(inpnums(inpstats));

prompts = cell(inpnum,1);
dlgtitle = 'User Input';
nlines = 1;
definps = repmat({''},[inpnum,1]);
options.Resize = 'on';
options.WindowStyle = 'normal';
options.Interpreter = 'tex';
pmtid = 1;

% runopts.rtrrfq
if inpset.rtrrfq
   if any(runopts.tst == [1,2])
      prompts{pmtid} = ...
         'Number of fixed time steps to take between each runtime recording:';
   else
      error('Fixed steps required for runtime recording!!!');
   end;
   pmtid = pmtid + 1;
end;

% runopts.tit
if inpset.tit
   prompts{pmtid} = 'Time integrator type:';
   pmtid = pmtid + 1;
end;

% runopts.ii
if inpset.ii
   prompts{pmtid} = 'Beginning of integration interval:';
   pmtid = pmtid + 1;
   prompts{pmtid} = 'End of integration interval:';
   pmtid = pmtid + 1;
end;

% runopts.fts
if inpset.fts
   if runopts.ftsspect == 0
      prompts{pmtid} = 'Number of equal fixed time steps to take:';
   elseif runopts.ftsspect == 1
      prompts{pmtid} = 'Size of one fixed time step:';
   end;
end;

% Collecting and interpreting user input
if inpnum > 0
   inps = inputdlg(prompts,dlgtitle,nlines,definps,options);
   
   pmtid = 1;
   if inpset.rtrrfq
      runopts.rtrrfq = uint32(str2double(inps{pmtid}));
      pmtid = pmtid + 1;
   end;
   if inpset.tit
      runopts.tit = inps{pmtid};
      pmtid = pmtid + 1;
   end;
   if inpset.ii
      runopts.ii = [str2double(inps{pmtid}) ...
         str2double(inps{pmtid+1})];
      pmtid = pmtid + 2;
   end;
   if inpset.fts
      if runopts.ftsspect == 0
         runopts.fts = uint32(str2double(inps{pmtid}));
      elseif runopts.ftsspect == 1
         runopts.fts = str2double(inps{pmtid});
      end;
   end;
end;
end

%% File creation function

function fids = createfiles(runopts)
% The purpose of this function is to check which files to create, create those
% files, and return the file IDs in a data structure.
%
% INPUTS/OUTPUTS
% runopts - Runtime and Recording Options data structure.
% fids - File IDs data structure. May be empty, if no files are required.

% Instantiating fids
fids = struct();

% Log file
if any(runopts.verbose == [1,2]) && any(runopts.log == [1,2])
   fids.logfid = fopen(runopts.logfn,'w');
end;

% Runtime results file
if runopts.rti && runopts.rtrr
   fids.rtrfid = fopen(runopts.rtrrfn,'w');
end;
end