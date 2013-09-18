%% completeBGSDdata
% BGSD completion function.
%
% [runopts,bgs] = completeBGSDdata(runopts,bgs,['opmode',opmode,
%     'valcheck',valcheck])
%
% This function accepts the run options and bond graph system data structures,
% and either converts them into a form fit for writing the BGSD file, or into a
% form fit for processing.
%
% See individual function descriptions or the documentation for usage syntax.
%
% Package:    BGSolver v1.03
% Subpackage: BGSD_Generator
% Date:       January 14, 2013
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Variable descriptions
%
% - On input, multiple fields can be of several different data types (i.e.,
%   boolean or int8). On output, their data types will be converted by the
%   convertDT function, depending on the operation mode.
% - Certain required variables may be empty for certain BGSs. In these cases,
%   the empty variables must still be present in the bgs input.
%
% INPUTS/OUTPUTS
% The data structures' specifications for this function are more flexible than
% for readBGSD or writeBGSD. The input specifications for this function are
% identical to those for fBGSolve. More data types are allowed for this function
% than for readBGSD/writeBGSD, and more fields can be optional. The returned
% data structures fit the more stringent specifications of readBGSD and
% writeBGSD.
% The possible field values are listed in the BGSD file format documentation.
% The default values are listed there as well. Refer to the documentation for
% more information. Note, that all multi-option fields, like runopts.verbose,
% are specified using int8. This is because integers are not subject to the
% roundoff error, and so integers of different types can be compared.
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
%  runopts.fts - Fixed time step. Either uint32 (number of equally-sized time
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
% ('opmode',opmode)
% Operation mode for this function. string. Possible values are:
% 	'writeBGSD' - Mode 1 (Default mode).
%  'BGSolve' - Mode 2.
%
% ('valcheck',valcheck)
% Whether or not to check input for validity. boolean. Possible values are:
%  false - Do not check the data structures for validity (Default).
%  true - Check the data structures for validity.
%
% STANDARD VARIABLES
% eid - Element ID; used to refer to elementwise vectors.
% bid - Bond ID; used to refer to bondwise vectors.
% pid - Port ID; used to refer to elements' port vectors.

%% Master BGSD completion function

function [runopts,bgs] = completeBGSDdata(varargin)
% The purpose of this function is to complete the data structures supplied to
% it, in one of two modes. In both modes, missing optional information is first
% filled into the data structures, using either the default settings, or by
% deducing the missing information from the information available in the data
% structures. If any required information is missing, an error is thrown.
% Optionally, the function can also check for validity of information input.
%
% MODES
% Mode 1: The data structures are formatted to be passed to the writeBGSD
%         function which creates a BGSD file based on the two data structures.
%         Element list storage type b is used.
% Mode 2: The data structures are formatted to be passed to the fBGSolve
%         function which processes the problem described by the two data
%         structures. Element list storage type a is used.
%
% INPUTS/OUTPUTS
% runopts - Runtime and Recording Options data structure (may be empty).
% bgs - Bond Graph System data structure.
% opmode - Mode of operations.
% valcheck - Whether or not to check data structures' validity.
% See the variable description cell for more information.

% Constructing input parser
ip = buildIP();
% Parsing input
ip.parse(varargin{:});
% Retrieving inputs
runopts = ip.Results.runopts;
bgs = ip.Results.bgs;
opmode = ip.Results.opmode;
valcheck = ip.Results.valcheck;
% Clearing memory
clear('ip');

% Loading BGSD defaults
load('BGSD_defaults.mat');

% Making sure validity check is not required
assert(~valcheck,'Current revision does not support validity checking!!!');

% Checking element list storage type
if isfield(bgs,'etypes')
   eltype = 'a';
elseif isfield(bgs,'elements')
   eltype = 'b';
else
   error('Unknown element list storage type!!!');
end;

% Making sure appropriate element list storage type is used
assert(strcmp(eltype,'b'),...
   'Current revision does not support element list storage type a!!!');

% Making bgs.elements vertical
bgs.elements = reshape(bgs.elements,[],1);

% Setting element defaults
bgs = setedefaults(bgs,defaults);

% Converting data types
[runopts,bgs] = convertDTs(runopts,bgs,opmode);

% Converting element list storage type
% If processing for writing, EL stays in type b.
% If processing for processing, EL converts to type a.
if strcmp(opmode,'writeBGSD')
   bgs = convertELtype(bgs,eltype,'b');
elseif strcmp(opmode,'BGSolve')
   bgs = convertELtype(bgs,eltype,'a');
end;

% Setting defaults
[runopts,bgs] = setdefaults(runopts,bgs,defaults);

% Trimming unneeded fields in runopts
runopts = trimrunopts(runopts);
end

%% Input parser construction function

function ip = buildIP()
% The purpose of this function is to construct an input parser object for the
% completeBGSDdata function.
%
% INPUTS/OUTPUTS
% ip - Input parser object for the completeBGSDdata function.

validOpModes = {'writeBGSD','BGSolve'};

ip = inputParser;
ip.StructExpand = false;
ip.addRequired('runopts',@isstruct);
ip.addRequired('bgs',@isstruct);
ip.addParamValue('opmode','writeBGSD',@(str)any(strcmp(str,validOpModes)));
ip.addParamValue('valcheck',false,@islogical);
end

%% Data type conversion function

function [runopts,bgs] = convertDTs(runopts,bgs,opmode)
% The purpose of this function is to convert the data types of the fields in the
% data structures supplied to it, depending on the mode of operation. Function
% handles are all converted to strings, which may produce errors if the handles
% are handled incorrectly.
%
% INPUTS/OUTPUTS
% runopts - Runtime and Recording Options data structure.
% bgs - Bond Graph System data structure.
% opmode - Mode of operation.
% Type b element list storage type and valid data structures are assumed.
% See the variable description cell for more information.

% Mode-independent data type conversion
% SECTION 1 - File Header.
fieldnames1 = {'svc','svf','ifnotes'};
dtypes1 = {'double','double','logical'};
runopts = convertFieldDTs(runopts,fieldnames1,dtypes1);

% SECTION 2 - Runtime and Recording Options.
fieldnames2 = {'verbose','rti','log','logfnspec','smr','smfnspec','rtrr',...
   'rtrrfqspec','rtrrfq','rtrrfnspec','rrm','rrmfnspec'};
dtypes2 = {'int8','logical','int8','int8','logical','int8','logical','int8',...
   'uint32','int8','logical','int8'};
runopts = convertFieldDTs(runopts,fieldnames2,dtypes2);

% SECTION 3 - Time Integration Options.
fieldnames3 = {'titspec','iispec','ii','tst','ftsspect','ftsspec','phxsplit',...
   'addtio','tidisp'};
dtypes3 = {'int8','int8','double','int8','int8','int8','int8','uint8','int8'};
runopts = convertFieldDTs(runopts,fieldnames3,dtypes3);

% SECTION 4 - Solver and Jumpstart Options.
fieldnames4 = {'stolx','stolres','sdisp','srscal','sxscal','sxcent','spc',...
   'sigpdct','sjacob','sjsp','addso','jstype','jsmsnum','jsphxsplit',...
   'jstolx','jstolres','jsdisp','jsrscal','jsxscal','jsxcent','jspc',...
   'jsjacob','jsjsp','addjso'};
dtypes4 = {'double','double','int8','int8','int8','int8','int8','int8',...
   'int8','int8','uint8','int8','uint32','int8','double','double','int8',...
   'int8','int8','int8','int8','int8','int8','uint8'};
runopts = convertFieldDTs(runopts,fieldnames4,dtypes4);

% SECTION 5 - Evaluator Options.
fieldnames5 = {'cs2n','nft8'};
dtypes5 = {'logical','int8'};
runopts = convertFieldDTs(runopts,fieldnames5,dtypes5);

% SECTION 6 - Post-Processor Options.
fieldnames6 = {'epp','fpp','ppdisp'};
dtypes6 = {'int32','int32','int8'};
runopts = convertFieldDTs(runopts,fieldnames6,dtypes6);
if isfield(runopts,'eppIDs') && islogical(runopts.eppIDs)
   runopts.eppIDs = uint16(find(runopts.eppIDs));
elseif isfield(runopts,'eppIDs') && ~islogical(runopts.eppIDs)
   runopts.eppIDs = uint16(runopts.eppIDs);
end;
if isfield(runopts,'fppIDs') && islogical(runopts.fppIDs)
   runopts.fppIDs = uint16(find(runopts.fppIDs));
elseif isfield(runopts,'fppIDs') && ~islogical(runopts.fppIDs)
   runopts.fppIDs = uint16(runopts.fppIDs);
end;

% SECTION 7 - Linked Files and Directories.
fieldnames7 = {'lfnum','ldnum'};
dtypes7 = {'uint8','uint8'};
runopts = convertFieldDTs(runopts,fieldnames7,dtypes7);

% SECTION 8 - Additional System Information.

% SECTION 9 - Bond Connectivity Map.
fieldnames9 = 'bnum';
dtypes9 = 'uint16';
bgs = convertFieldDTs(bgs,fieldnames9,dtypes9);
if isstruct(bgs.bonds)
   bgs.bonds = uint16([reshape(bgs.bonds.to,[],1),...
      reshape(bgs.bonds.from,[],1)]);
else
   bgs.bonds = uint16(bgs.bonds);
end;

% SECTION 10 - Element List.
% Assuming element list storage type b.
bgs.enum = uint16(bgs.enum);
fieldnames10 = {'x0','epnum','ebcon','caus','phx','nn','at'};
dtypes10 = {'double','uint16','uint16','int8','uint8','logical','double'};
for eid = 1:bgs.enum
   bgs.elements(eid) = convertFieldDTs(bgs.elements(eid),fieldnames10,dtypes10);
   % bgs.elements(eid).expr reshape
   if any(strcmp(bgs.elements(eid).etype,{'R2','MR2','RN','MRN'})) && ...
         any(strcmp(bgs.elements(eid).exprtype,{'CC','CMC','NMC'}))
      bgs.elements(eid).expr = reshape(bgs.elements(eid).expr,...
         sqrt(round(numel(bgs.elements(eid).expr))),[]).';
   end;
end;

% Mode-dependent conversions.
% Assuming element list storage type b.
if strcmp(opmode,'writeBGSD')
   % runopts.tit
   if isfield(runopts,'tit') && isa(runopts.tit,'function_handle')
      runopts.tit = func2str(runopts.tit);
   end;
   % runopts.addtiolist
   runopts = convertAddOptDTs(runopts,'addtiolist');
   % runopts.addsolist
   runopts = convertAddOptDTs(runopts,'addsolist');
   % runopts.addjsolist
   runopts = convertAddOptDTs(runopts,'addjsolist');
   
   % bgs.elements(eid).mvars
   for eid = 1:bgs.enum
      if ~isempty(bgs.elements(eid).mvars) && isa(bgs.elements(eid).mvars,'sym')
         bgs.elements(eid).mvars = ...
            reshape(convertSym2Str(bgs.elements(eid).mvars),[],1);
      end;
   end;
   % bgs.elements(eid).expr
   for eid = 1:bgs.enum
      if ~isempty(bgs.elements(eid).expr)
         if any(strcmp(bgs.elements(eid).exprtype,{'NMC','NE','NME'}))
            % NMC, NE, NME
            if iscell(bgs.elements(eid).expr) && ...
                  ~iscellstr(bgs.elements(eid).expr)
               bgs.elements(eid).expr = cellfun(@func2str,...
                  bgs.elements(eid).expr,'UniformOutput',false);
            elseif ~iscell(bgs.elements(eid).expr) && ...
                  isa(bgs.elements(eid).expr,'function_handle')
               bgs.elements(eid).expr = func2str(bgs.elements(eid).expr);
            end;
         elseif any(strcmp(bgs.elements(eid).exprtype,{'CMC','CE','CME'}))
            % CMC, CE, CME
            if isa(bgs.elements(eid).expr,'sym')
               bgs.elements(eid).expr = ...
                  convertSym2Str(bgs.elements(eid).expr,true);
            end;
         else
            % CC
            bgs.elements(eid).expr = double(bgs.elements(eid).expr);
         end;
      end;
   end;
elseif strcmp(opmode,'BGSolve')
   % bgs.elements(eid).mvars
   for eid = 1:bgs.enum
      if ~isempty(bgs.elements(eid).mvars) && ...
            iscellstr(bgs.elements(eid).mvars)
         bgs.elements(eid).mvars = ...
            reshape(sym(bgs.elements(eid).mvars),[],1);
      end;
   end;
   % bgs.elements(eid).expr
   for eid = 1:bgs.enum
      if ~isempty(bgs.elements(eid).expr)
         if any(strcmp(bgs.elements(eid).exprtype,{'CMC','CE','CME'}))
            % CMC, CE, CME
            if iscellstr(bgs.elements(eid).expr) || ...
                  ischar(bgs.elements(eid).expr)
               bgs.elements(eid).expr = sym(bgs.elements(eid).expr);
            end;
         else
            % CC - do nothing
         end;
      end;
   end;
else
   error('Invalid operation mode!!!');
end;
end

%% Field data type conversion function

function ds = convertFieldDTs(ds,fieldnames,dtypes)
% The purpose of this function is to convert the data types of an array of
% fields of a data structure. The function is vectorized.
%
% INPUTS/OUTPUTS
% ds - Data structure which contains the fields to be converted.
% fieldnames - Cell vector of field names whose data types to convert, OR an
%              individual field;
%              Some of these fields may not exist.
% dtypes - Cell vector of the new data types OR an individual data type.

if iscell(fieldnames)
   F = length(fieldnames); % Number of field names
   for f = 1:F
      fieldname = fieldnames{f};
      if isfield(ds,fieldname)
         ds.(fieldname) = feval(dtypes{f},ds.(fieldname));
      end;
   end;
else
   if isfield(ds,fieldnames)
      ds.(fieldnames) = feval(dtypes,ds.(fieldnames));
   end;
end;
end

%% Additional options list conversion function

function runopts = convertAddOptDTs(runopts,fieldname)
% The purpose of this function is to convert the parameter values in an
% additional options list to data types appropriate for writeBGSD operation
% mode.
%
% INPUTS/OUTPUTS
% runopts - Runtime and Recording Options data structure.
% fieldname - Additional options list field name.

if isfield(runopts,fieldname) && size(runopts.(fieldname),1) > 0
   for addoID = 1:size(runopts.(fieldname),1)
      addoval = runopts.(fieldname){addoID,2};
      
      if ischar(addoval)
         % Do nothing
      elseif isnumeric(addoval) && isscalar(addoval)
         % Do nothing
      elseif isa(addoval,'function_handle')
         runopts.(fieldname){addoID,2} = func2str(addoval);
      else
         if strcmp(fieldname,'addtiolist')
            optname = 'TIO';
         elseif strcmp(fieldname,'addsolist')
            optname = 'SO';
         elseif strcmp(fieldname,'addjsolist')
            optname = 'JSO';
         end;
         
         error('Additional %s %u is of unsupported type!!!',optname,addoID);
      end;
   end;
end;
end

%% Symbolic array to string cell array conversion function

function C = convertSym2Str(S,opt)
% The purpose of this function is to convert a symbolic array to an equally
% sized string cell array. Optionally, if the symbolic array is 1x1, the outcome
% can be converted to a single string.
%
% INPUTS/OUTPUTS
% S - Symbolic array.
% opt - Whether or not to convert a 1x1 cell array into a string. Default:
%       false.
% C - String cell array (or single string).

if nargin == 1
   opt = false;
end;

C = arrayfun(@char,S,'UniformOutput',false);

if numel(C) == 1 && opt
   C = C{1};
end;
end

%% Element list storage type conversion function

function bgs = convertELtype(bgs,oldtype,newtype)
% The purpose of this function is to convert the element list storage type in a
% given BGS. If the old and new types are the same, nothing happens. Data types
% are assumed to have been converted. See the variable description cell for more
% information about the EL storage types.
%
% INPUTS/OUTPUTS
% bgs - Bond Graph System data structure.
% oldtype - Old EL storage type.
% newtype - New EL storage type.
% See the variable description cell for more information.

if strcmp(oldtype,'a') && strcmp(newtype,'b')
   % Type a EL storage input is not supported in this revision
elseif strcmp(oldtype,'b') && strcmp(newtype,'a')
   % Numbers and IDs
   if isfield(bgs,'enum')
      bgs.enum = uint16(length(bgs.elements));
   end;
   [bgs.exprnum,bgs.exprIDs] = getSIDs(bgs.elements,'expr');
   [bgs.emodnum,bgs.modIDs] = getSIDs(bgs.elements,'mvars');
   [bgs.epambnum,bgs.pambIDs] = getSIDs(bgs.elements,'epnum');
   [bgs.ebambnum,bgs.bambIDs] = getSIDs(bgs.elements,'ebcon');
   [bgs.ecambnum,bgs.cambIDs] = getSIDs(bgs.elements,'caus');
   [bgs.ennum,bgs.numIDs,bgs.nvnum] = getNEIDs(bgs.elements);
   [bgs.xnum,bgs.x0IDs] = getSIDs(bgs.elements,'x0');
   [bgs.capnum,bgs.capIDs] = getValSIDs(bgs.elements,'etype',{'C','MC'});
   [bgs.inertnum,bgs.inertIDs] = getValSIDs(bgs.elements,'etype',{'I','MI'});
   [bgs.phxspecnum,bgs.phxIDs] = getSIDs(bgs.elements,'phx');
   [bgs.nnspecnum,bgs.nnIDs] = getSIDs(bgs.elements,'nn');
   [bgs.atspecnum,bgs.atIDs] = getSIDs(bgs.elements,'at');
   
   % Element data
   bgs.etypes = {bgs.elements.etype}.';
   bgs.exprtypes = {bgs.elements(bgs.exprIDs).exprtype};
   bgs.exprs = {bgs.elements(bgs.exprIDs).expr};
   bgs.mvars = {bgs.elements(bgs.modIDs).mvars};
   bgs.epnums = [bgs.elements(bgs.pambIDs).epnum];
   bgs.ebcons = {bgs.elements(bgs.bambIDs).ebcon};
   bgs.causes = {bgs.elements(bgs.cambIDs).caus};
   bgs.x0s = [bgs.elements(bgs.x0IDs).x0];
   bgs.phxs = [bgs.elements(bgs.phxIDs).phx];
   bgs.nns = [bgs.elements(bgs.nnIDs).nn];
   bgs.ats = [bgs.elements(bgs.atIDs).at];
   
   % Making element data vectors vertical
   bgs.etypes = reshape(bgs.etypes,[],1);
   bgs.exprtypes = reshape(bgs.exprtypes,[],1);
   bgs.exprs = reshape(bgs.exprs,[],1);
   bgs.mvars = reshape(bgs.mvars,[],1);
   bgs.epnums = reshape(bgs.epnums,[],1);
   bgs.ebcons = reshape(bgs.ebcons,[],1);
   bgs.causes = reshape(bgs.causes,[],1);
   bgs.x0s = reshape(bgs.x0s,[],1);
   bgs.phxs = reshape(bgs.phxs,[],1);
   bgs.nns = reshape(bgs.nns,[],1);
   bgs.ats = reshape(bgs.ats,[],1);
   
   % Clearing type b storage
   bgs = rmfield(bgs,'elements');
else
   % No conversion, do nothing
end;
end

%% Field-based struct ID retrieval function

function [SIDnum,SIDs] = getSIDs(S,fieldname)
% The purpose of this function is to retrieve the number of nonempty elements of
% a struct array whose fields of a given name are nonempty, and, optionally,
% their indices. The size of the SID vector is returned as uint16. The struct
% indices (denoted SIDs, Struct IDs) are optionally returned, also as uint16
% vector.
%
% INPUTS/OUTPUTS
% S - Struct array.
% fieldname - Name of the field to check for emptiness.
% SIDnum - Number of SIDs returned, uint16.
% SIDs - Struct ID vector, uint16.

assert(isfield(S,fieldname),'Non-existant field being retrieved!!!');
SIDsBoolean = ~cellfun('isempty',{S.(fieldname)});
SIDs = uint16(find(SIDsBoolean));
SIDnum = uint16(length(SIDs));
end

%% Numeric element and variable ID retrieval function

function [ennum,numIDs,nvnum] = getNEIDs(elements)
% The purpose of this function is to retrieve the number of numeric elements,
% their EIDs, and the number of numeric variables in a BGS. The numbers are
% returned as uint16s, and the indices are returned as uint16 vector. Element
% list storage type b is assumed.
%
% INPUTS/OUTPUTS
% elements - Element List data structure in a BGS, storage tybe b.
% ennum - Numeric element number, uint16.
% numIDs - Numeric elements' EID vector, uint16.
% nvnum - Numeric variable number, uint16.

% Extracting fields
etypes = {elements.etype};
exprtypes = {elements.exprtype};
pnums = {elements.epnum};

% Constructing evaluation functions
isNMC = @(expr) strcmp(expr,'NMC');
isNE = @(expr) any(strcmp(expr,{'NE','NME'}));
isR2 = @(etype) any(strcmp(etype,{'R2','MR2'}));
isRN = @(etype) any(strcmp(etype,{'RN','MRN'}));

% Constructing Boolean index vectors
NMCIDsBoolean = cellfun(isNMC,exprtypes);
NEIDsBoolean = cellfun(isNE,exprtypes);
numIDsBoolean = NMCIDsBoolean | NEIDsBoolean;
R2IDsBoolean = cellfun(isR2,etypes);
RNIDsBoolean = cellfun(isRN,etypes);
R2NIDsBoolean = R2IDsBoolean | RNIDsBoolean;

% Correcting port number vector
pnums(R2IDsBoolean) = {uint16(2)};
pnums(~R2NIDsBoolean) = {uint16(0)};
pnums = uint16(cell2mat(pnums));

% Evaluating numbers and ID vector
numIDs = uint16(find(numIDsBoolean));
ennum = uint16(length(numIDs));
nvnum = uint16(sum(pnums(R2NIDsBoolean & NEIDsBoolean)) + ...
   sum(pnums(R2NIDsBoolean & NMCIDsBoolean).^2) + ...
   nnz(~R2NIDsBoolean & numIDsBoolean));
end

%% Field value-based struct ID retrieval function

function [SIDnum,SIDs] = getValSIDs(S,fieldname,vals,fcomp)
% The purpose of this function is to retrieve the number of elements of a struct
% array whose fields of a given name have one of the given values, and,
% optionally, their indices. Optionally, the equality checking function for the
% values is input as well. By default, case-sensitive string equality checking
% is used. The size of the SID vector is returned as uint16. The struct indices
% (denoted SIDs, Struct IDs) are optionally returned, also as uint16 vector.
%
% INPUTS/OUTPUTS
% S - Struct array.
% fieldname - Name of the field to check the value of.
% vals - Possible field values to include in the SID list, vector of
%        corresponding type. Default: string cell array.
% fcomp - Function to check the field value equality. OPTIONAL.
% SIDnum - Number of SIDs returned, uint16.
% SIDs - Struct ID vector, uint16.

if nargin == 3
   fcomp = @(str) any(strcmp(str,vals));
elseif nargin == 4
   % fcomp specified in the input
else
   error('Invalid number of inputs!!!');
end;

assert(isfield(S,fieldname),'Non-existant field being retrieved!!!');
SIDsBoolean = cellfun(fcomp,{S.(fieldname)});
SIDs = uint16(find(SIDsBoolean));
SIDnum = uint16(length(SIDs));
end

%% Elementwise default value setting function

function bgs = setedefaults(bgs,defaults)
% The purpose of this function is to set several empty element fields in the
% bgs.elements data structure using default and deduced values. Element list
% storage type b is assumed.
%
% INPUTS/OUTPUTS
% bgs - Bond Graph System data structure.
% defaults - Default values data structure.
% See the variable description cell for more information.

% Setting element number
if ~isfield(bgs,'enum')
   bgs.enum = uint16(length(bgs.elements));
end;

% Setting element fields
% .epnum, .phx, .nn, .at
for eid = 1:bgs.enum
   etype = bgs.elements(eid).etype;
   
   % .epnum
   if isempty(bgs.elements(eid).epnum) && any(strcmp(etype,{'RN','MRN'}))
      exprtype = bgs.elements(eid).exprtype;
      expr = bgs.elements(eid).expr;
      if (strcmp(etype,'RN') && strcmp(exprtype,'CC')) || ...
            (strcmp(etype,'MRN') && any(strcmp(exprtype,{'CMC','NMC'})))
         bgs.elements(eid).epnum = uint16(round(sqrt(numel(expr))));
      elseif (strcmp(etype,'RN')) && any(strcmp(exprtype,{'CE','NE'})) || ...
            (strcmp(etype,'MRN') && any(strcmp(exptype,{'CME','NME'})))
         bgs.elements(eid).epnum = uint16(numel(expr));
      end;
   end;
   
   % .phx
   if isempty(bgs.elements(eid).phx) && any(strcmp(etype,{'C','I'}))
      bgs.elements(eid).phx = defaults.bgs.element.phx;
   end;
   
   % .nn
   if isempty(bgs.elements(eid).nn) && any(strcmp(etype,{'C','I','MC','MI'}))
      bgs.elements(eid).nn = defaults.bgs.element.nn;
   end;
   
   % .at
   if isempty(bgs.elements(eid).at) && any(strcmp(etype,{'C','I','MC','MI'}))
      bgs.elements(eid).at = defaults.bgs.element.at;
   end;
end;
end

%% Default value setting function

function [runopts,bgs] = setdefaults(runopts,bgs,defaults)
% The purpose of this function is to set several empty fields in the data
% structures using default and deduced values. Data types and the element list
% storage type are assumed to have been converted according to the mode of
% operations, and the remaining data types are converted after all leftover
% default values are filled in.
%
% INPUTS/OUTPUTS
% runopts - Runtime and Recording Options data structure.
% bgs - Bond Graph System data structure.
% defaults - Default values data structure.
% Valid data structures and EL storage type appropriate for the mode are
% assumed. See the variable description cell for more information.

% Mode-independent settings.
%
% SECTION 1 - File Header.
%
% runopts.svc and runopts svf
if ~isfield(runopts,'svc') && ~isfield(runopts,'svf')
   runopts.svc = defaults.runopts.svc;
   runopts.svf = defaults.runopts.svf;
elseif ~isfield(runopts,'svc') && isfield(runopts,'svf')
   runopts.svc = runopts.svf;
elseif isfield(runopts,'svc') && ~isfield(runopts,'svf')
   runopts.svf = runopts.svc;
end;
%
% runopts.ifnotes, runopts.notes
if ~isfield(runopts,'ifnotes') && isfield(runopts,'notes')
   runopts.ifnotes = true;
elseif isfield(runopts,'ifnotes') && ~isfield(runopts,'notes')
   runopts.notes = defaults.runopts.notes;
elseif ~isfield(runopts,'ifnotes') && ~isfield(runopts,'notes')
   runopts.ifnotes = defaults.runopts.ifnotes;
   runopts.notes = defaults.runopts.notes;
end;

% SECTION 2 - Runtime and Recording Options.
%
% runopts.verbose
if ~isfield(runopts,'verbose')
   runopts.verbose = defaults.runopts.verbose;
end;
%
% runopts.rti
if ~isfield(runopts,'rti')
   runopts.rti = defaults.runopts.rti;
end;
%
% runopts.log
if ~isfield(runopts,'log')
   runopts.log = defaults.runopts.log;
end;
%
% runopts.logfnspec
if ~isfield(runopts,'logfnspec')
   if isfield(runopts,'logfn')
      runopts.logfnspec = int8(1);
   else
      runopts.logfnspec = defaults.runopts.logfnspec;
   end;
end;
%
% runopts.smr
if ~isfield(runopts,'smr')
   runopts.smr = defaults.runopts.smr;
end;
%
% runopts.smfnspec
if ~isfield(runopts,'smfnspec')
   if isfield(runopts,'smfn')
      runopts.smfnspec = int8(1);
   else
      runopts.smfnspec = defaults.runopts.smfnspec;
   end;
end;
%
% runopts.rtrr
if ~isfield(runopts,'rtrr')
   runopts.rtrr = defaults.runopts.rtrr;
end;
%
% runopts.rtrrfqspec
if ~isfield(runopts,'rtrrfqspec')
   if isfield(runopts,'rtrrfq')
      runopts.rtrrfqspec = int8(1);
   else
      runopts.rtrrfqspec = defaults.runopts.rtrrfqspec;
   end;
end;
%
% runopts.rtrrfnspec
if ~isfield(runopts,'rtrrfnspec')
   if isfield(runopts,'rtrrfn')
      runopts.rtrrfnspec = int8(1);
   else
      runopts.rtrrfnspec = defaults.runopts.rtrrfnspec;
   end;
end;
%
% runopts.rrm
if ~isfield(runopts,'rrm')
   runopts.rrm = defaults.runopts.rrm;
end;
%
% runopts.rrmfnspec
if ~isfield(runopts,'rrmfnspec')
   if isfield(runopts,'rrmfn')
      runopts.rrmfnspec = int8(1);
   else
      runopts.rrmfnspec = defaults.runopts.rrmfnspec;
   end;
end;

% SECTION 3 - Time Integrator Options.
%
% runopts.titspec
if ~isfield(runopts,'titspec')
   runopts.titspec = defaults.runopts.titspec;
end;
%
% runopts.tit
if ~isfield(runopts,'tit')
   runopts.tit = defaults.runopts.tit;
end;
%
% runopts.iispec
if ~isfield(runopts,'iispec')
   if isfield(runopts,'ii')
      runopts.iispec = int8(0);
   else
      runopts.iispec = defaults.runopts.iispec;
   end;
end;
%
% runopts.tst
if ~isfield(runopts,'tst')
   runopts.tst = int8(0);
end;
%
% runopts.ftsspect
if ~isfield(runopts,'ftsspect')
   runopts.ftsspect = defaults.runopts.ftsspect;
end;
%
% runopts.ftsspec
if ~isfield(runopts,'ftsspec')
   if isfield(runopts,'fts')
      runopts.ftsspec = int8(0);
   else
      runopts.ftsspec = defaults.runopts.ftsspec;
   end;
end;
%
% runopts.phxsplit
if ~isfield(runopts,'phxsplit')
   runopts.phxsplit = defaults.runopts.phxsplit;
end;
%
% runopts.addtio
if ~isfield(runopts,'addtio')
   if isfield(runopts,'addtiolist')
      runopts.addtio = uint8(size(runopts.addtiolist,1));
   else
      runopts.addtio = uint8(0);
   end;
end;
%
% runopts.tidisp
if ~isfield(runopts,'tidisp')
   runopts.tidisp = defaults.runopts.tidisp;
end;

% SECTION 4 - Solver and Jumpstart Options.
%
% runopts.stolx
if ~isfield(runopts,'stolx')
   runopts.stolx = defaults.runopts.stolx;
end;
%
% runopts.stolres
if ~isfield(runopts,'stolres')
   runopts.stolres = defaults.runopts.stolres;
end;
%
% runopts.sdisp
if ~isfield(runopts,'sdisp')
   runopts.sdisp = defaults.runopts.sdisp;
end;
%
% runopts.srscal
if ~isfield(runopts,'srscal')
   runopts.srscal = defaults.runopts.srscal;
end;
%
% runopts.sxscal
if ~isfield(runopts,'sxscal')
   runopts.sxscal = defaults.runopts.sxscal;
end;
%
% runopts.sxcent
if ~isfield(runopts,'sxcent')
   runopts.sxcent = defaults.runopts.sxcent;
end;
%
% runopts.spc
if ~isfield(runopts,'spc')
   runopts.spc = defaults.runopts.spc;
end;
%
% runopts.sigpdct
if ~isfield(runopts,'sigpdct')
   runopts.sigpdct = defaults.runopts.sigpdct;
end;
%
% runopts.sjacob
if ~isfield(runopts,'sjacob')
   runopts.sjacob = defaults.runopts.sjacob;
end;
%
% runopts.sjsp
if ~isfield(runopts,'sjsp')
   runopts.sjsp = defaults.runopts.sjsp;
end;
%
% runopts.addso
if ~isfield(runopts,'addso')
   if isfield(runopts,'addsolist')
      runopts.addso = uint8(size(runopts.addsolist,1));
   else
      runopts.addso = uint8(0);
   end;
end;
%
% runopts.jstype
if ~isfield(runopts,'jstype')
   runopts.jstype = defaults.runopts.jstype;
end;
%
% runopts.jsmsnum
if ~isfield(runopts,'jsmsnum')
   runopts.jsmsnum = defaults.runopts.jsmsnum;
end;
%
% runopts.jsphxsplit
if ~isfield(runopts,'jsphxsplit')
   runopts.jsphxsplit = defaults.runopts.jsphxsplit;
end;
%
% runopts.jstolx
if ~isfield(runopts,'jstolx')
   runopts.jstolx = defaults.runopts.jstolx;
end;
%
% runopts.jstolres
if ~isfield(runopts,'jstolres')
   runopts.jstolres = defaults.runopts.jstolres;
end;
%
% runopts.jstit
if ~isfield(runopts,'jstit')
   runopts.jstit = defaults.runopts.jstit;
end;
%
% runopts.jsdisp
if ~isfield(runopts,'jsdisp')
   runopts.jsdisp = defaults.runopts.jsdisp;
end;
%
% runopts.jsrscal
if ~isfield(runopts,'jsrscal')
   runopts.jsrscal = defaults.runopts.jsrscal;
end;
%
% runopts.jsxscal
if ~isfield(runopts,'jsxscal')
   runopts.jsxscal = defaults.runopts.jsxscal;
end;
%
% runopts.jsxcent
if ~isfield(runopts,'jsxcent')
   runopts.jsxcent = defaults.runopts.jsxcent;
end;
%
% runopts.jspc
if ~isfield(runopts,'jspc')
   runopts.jspc = defaults.runopts.jspc;
end;
%
% runopts.jsjacob
if ~isfield(runopts,'jsjacob')
   runopts.jsjacob = defaults.runopts.jsjacob;
end;
%
% runopts.jsjsp
if ~isfield(runopts,'jsjsp')
   runopts.jsjsp = defaults.runopts.jsjsp;
end;
%
% runopts.addjso
if ~isfield(runopts,'addjso')
   if isfield(runopts,'addjsolist')
      runopts.addjso = uint8(size(runopts.addjsolist,1));
   else
      runopts.addjso = uint8(0);
   end;
end;

% SECTION 5 - Evaluator Options.
%
% runopts.cs2n
if ~isfield(runopts,'cs2n')
   runopts.cs2n = defaults.runopts.cs2n;
end;
%
% runopts.nft
if ~isfield(runopts,'nft')
   runopts.nft = defaults.runopts.nft;
end;

% SECTION 6 - Post-Processor Options.
%
% runopts.epp
if ~isfield(runopts,'epp')
   if isfield(runopts,'eppIDs')
      runopts.epp = int32(length(runopts.eppIDs));
   else
      runopts.epp = defaults.runopts.epp;
   end;
end;
% runopts.fpp
if ~isfield(runopts,'fpp')
   if isfield(runopts,'fppIDs')
      runopts.fpp = int32(length(runopts.fppIDs));
   else
      runopts.fpp = defaults.runopts.fpp;
   end;
end;
%
% runopts.ppdisp
if ~isfield(runopts,'ppdisp')
   runopts.ppdisp = defaults.runopts.ppdisp;
end;

% SECTION 7 - Linked Files and Directories.
%
% runopts.lfnum
if ~isfield(runopts,'lfnum')
   if isfield(runopts,'lflist')
      runopts.lfnum = uint8(length(runopts.lflist));
   else
      runopts.lfnum = defaults.runopts.lfnum;
   end;
end;
%
% runopts.ldnum
if ~isfield(runopts,'ldnum')
   if isfield(runopts,'ldlist')
      runopts.ldnum = uint8(length(runopts.ldlist));
   else
      runopts.ldnum = defaults.runopts.ldnum;
   end;
end;

% SECTION 8 - Additional System Information.
%
% bgs.asi
if ~isfield(bgs,'asi')
   bgs.asi = defaults.bgs.asi;
else
   if ~isfield(bgs.asi,'sci')
      bgs.asi.sci = defaults.bgs.asi.sci;
   end;
   if ~isfield(bgs.asi,'smi')
      bgs.asi.smi = defaults.bgs.asi.smi;
   end;
   if ~isfield(bgs.asi,'sli')
      bgs.asi.sli = defaults.bgs.asi.sli;
   end;
   if ~isfield(bgs.asi,'sri')
      bgs.asi.sri = defaults.bgs.asi.sri;
   end;
   if ~isfield(bgs.asi,'sei')
      bgs.asi.sei = defaults.bgs.asi.sei;
   end;
end;

% SECTION 9 - Bond Connectivity Map.
%
% bgs.bnum
if ~isfield(bgs,'bnum')
   bgs.bnum = uint16(size(bgs.bonds,1));
end;

% SECTION 10 - Element List.
%
% bgs.exprnum, bgs.emodnum, bgs.epambnum, bgs.ebambnum, bgs.ecambnum,
% bgs.ennum, bgs.nvnum, bgs.xnum, bgs.capnum, bgs.inertnum, bgs.phxspecnum,
% bgs.nnspecnum, bgs.atspecnum
if isfield(bgs,'elements')
   if ~isfield(bgs,'exprnum')
      bgs.exprnum = getSIDs(bgs.elements,'expr');
   end;
   if ~isfield(bgs,'emodnum')
      bgs.emodnum = getSIDs(bgs.elements,'mvars');
   end;
   if ~isfield(bgs,'epambnum')
      bgs.epambnum = getSIDs(bgs.elements,'epnum');
   end;
   if ~isfield(bgs,'ebambnum')
      bgs.ebambnum = getSIDs(bgs.elements,'ebcon');
   end;
   if ~isfield(bgs,'ecambnum')
      bgs.ecambnum = getSIDs(bgs.elements,'caus');
   end;
   if ~isfield(bgs,'ennum')
      [bgs.ennum,~,~] = getNEIDs(bgs.elements);
   end;
   if ~isfield(bgs,'nvnum')
      [~,~,bgs.nvnum] = getNEIDs(bgs.elements);
   end;
   if ~isfield(bgs,'xnum')
      bgs.xnum = getSIDs(bgs.elements,'x0');
   end;
   if ~isfield(bgs,'capnum')
      bgs.capnum = getValSIDs(bgs.elements,'etype',{'C','MC'});
   end;
   if ~isfield(bgs,'inertnum')
      bgs.inertnum = getValSIDs(bgs.elements,'etype',{'I','MI'});
   end;
   if ~isfield(bgs,'phxspecnum')
      bgs.phxspecnum = getSIDs(bgs.elements,'phx');
   end;
   if ~isfield(bgs,'nnspecnum')
      bgs.nnspecnum = getSIDs(bgs.elements,'nn');
   end;
   if ~isfield(bgs,'atspecnum')
      bgs.atspecnum = getSIDs(bgs.elements,'at');
   end;
end;

% Converting data types that are dependent on defaults.
%
% runopts.fts
if isfield(runopts,'fts')
   if runopts.ftsspect == 0
      runopts.fts = uint32(runopts.fts);
   elseif runopts.ftsspect == 1
      runopts.fts = double(runopts.fts);
   end;
end;
end

%% Runopts trimming function

function runopts = trimrunopts(runopts)
% The purpose of this function is to remove unneeded fields from the runopts
% data structure. Which fields are not needed is decided based on other fields'
% values; see the code documentation and the function code for more information.
%
% INPUTS/OUTPUTS
% runopts - Runtime and Recording Options data structure.

% SECTION 1 - File Header.
%
% runopts.notes
if ~runopts.ifnotes
   runopts = rmfield(runopts,'notes');
end;

% SECTION 2 - Runtime and Recording Options.
%
% runopts.log, runopts.logfnspec, runopts.logfn
if runopts.verbose == 0
   runopts = rmfieldcheck(runopts,{'log','logfnspec','logfn'});
else
   if runopts.log == 0
      runopts = rmfieldcheck(runopts,{'logfnspec','logfn'});
   elseif runopts.logfnspec ~= 1
      runopts = rmfieldcheck(runopts,'logfn');
   end;
end;
%
% runopts.smfnspec, runopts.smfn
if ~runopts.smr
   runopts = rmfieldcheck(runopts,{'smfnspec','smfn'});
elseif runopts.smfnspec ~= 1
   runopts = rmfieldcheck(runopts,'smfn');
end;
%
% runopts.rtrrfqspec, runopts.rtrrfq
if ~runopts.rtrr
   runopts = rmfieldcheck(runopts,{'rtrrfqspec','rtrrfq'});
elseif runopts.rtrrfqspec ~= 1
	runopts = rmfieldcheck(runopts,'rtrrfq');
end;
%
% runopts.rtrrfnspec, runopts.rtrrfn
if ~runopts.rtrr
   runopts = rmfieldcheck(runopts,{'rtrrfnspec','rtrrfn'});
elseif runopts.rtrrfnspec ~= 1
   runopts = rmfieldcheck(runopts,'rtrrfn');
end;
%
% runopts.rrmfnspec, runopts.rrmfn
if ~runopts.rrm
   runopts = rmfieldcheck(runopts,{'rrmfnspec','rrmfn'});
elseif runopts.rrmfnspec ~= 1
   runopts = rmfieldcheck(runopts,'rrmfn');
end;

% SECTION 3 - Time Integrator Options.
%
% runopts.tit
if runopts.titspec ~= 1
   runopts = rmfieldcheck(runopts,'tit');
end;
%
% runopts.ii
if runopts.iispec ~= 0
   runopts = rmfieldcheck(runopts,'ii');
end
%
% runopts.ftsspect, runopts.ftsspec, runopts.fts
if all(runopts.tst ~= [1,2])
   runopts = rmfieldcheck(runopts,{'ftsspect','ftsspec','fts'});
elseif runopts.ftsspec ~= 0
   runopts = rmfieldcheck(runopts,'fts');
end;

% SECTION 4 - Solver and Jumpstart Options.
%

% SECTION 5 - Evaluator Options.
%
% runopts.nft
if runopts.cs2n ~= 1
   runopts = rmfieldcheck(runopts,'nft');
end;

% SECTION 6 - Post-Processor Options.
%
% runopts.eppIDs
if runopts.epp <= 0
   runopts = rmfieldcheck(runopts,'eppIDs');
end;
%
% runopts.fppIDs
if runopts.fpp <= 0
   runopts = rmfieldcheck(runopts,'fppIDs');
end;

% SECTION 7 - Linked Files and Directories.
%
% runopts.lflist
if runopts.lfnum == 0
   runopts = rmfieldcheck(runopts,'lflist');
end;
%
% runopts.ldlist
if runopts.ldnum == 0
   runopts = rmfieldcheck(runopts,'ldlist');
end;
end