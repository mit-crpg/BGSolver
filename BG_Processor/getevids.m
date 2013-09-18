%% getevids
% Variable IDs extraction function for expression vectors.
%
% [[ift,xids,bvids,nids]] = getevids(E,vars,refs,
%              ['vars2check',vars2check,'opmode',opmode,'outform',outform])
%
% This function extracts variable IDs out of a vertical vector of symbolic
% expressions. The following variable IDs can be checked: time, state variables,
% bond variables and numeric variables. By default, time, state and numeric
% variables are checked. Variable order tracking cannot be enforced, because
% variables in symbolic expressions cannot have set orders. This function
% utilizes getvvids on each expression in the expression vector, and optional
% arguments that getvvids takes (other than vars2order) can be passed to it
% through this function.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Master variable ID extraction function for expression vectors

function varargout = getevids(varargin)
% The purpose of this function is to extract the corresponding variable IDs of
% variables present in the symbolic expression vector input to the function.
% Depending on requested outputs, not all variables are checked. The function
% works by running getvvids on every element of the expression vector.
%
% INPUTS/OUTPUTS
% E - Vertical symbolic expression vector to extract variable IDs from.
% vars - Symbolic variable data structure. Contains symbolic variables that are
%        used for comparison with variables in E and extraction.
% refs - Referencing arrays data structure.
% The optional input parameters are passed directly to getvvids, and are
% documented in getvvids documentation and comments. Order tracking is disabled,
% so vars2order cannot be used.
%
% ift - Full vertical boolean vector, with elements corresponding to whether or
%       not time is present in each expression in E.
% xids - State variable IDs present in E, formatted according to outform.
% bvids - Bond variable IDs present in E, formatted according to outform.
% nids - Numeric variable IDs present in E, formatted according to outform.
%
% Outputs are ordered the same way vars2check is ordered.

% Input parser construction
ip = buildIP();
% Input parsing
ip.parse(varargin{:});
% Retrieving inputs
E = ip.Results.E;
vars = ip.Results.vars;
refs = ip.Results.refs;
vars2check = ip.Results.vars2check;
outform = ip.Results.outform;
% Clearing memory
clear('ip');

% Constructing cell vector of variable vectors
exprnum = uint16(length(E));
V = cell(exprnum,1);
for exprid = 1:exprnum
   V{exprid} = symvar(E(exprid));
end;

% Getting variable type options
topts = getvtopts('t',vars2check);
xopts = getvtopts('x',vars2check);
bopts = getvtopts('b',vars2check);
nopts = getvtopts('n',vars2check);
vtnum = uint8(topts.ifcheck + xopts.ifcheck + bopts.ifcheck + nopts.ifcheck);

% Constructing output array
% Preallocating output array
varargout = cell(1,vtnum);

% Extracting variable IDs
[varargout{:}] = getvvids(V,vars,refs,'vars2check',vars2check,...
   'outform',outform);
end

%% Input parser construction function

function ip = buildIP()
% The purpose of this function is to construct an input parser object for
% the getevids function.
%
% INPUTS/OUTPUTS
% ip - Input parser object for the getvvids function.

fEcheck = @(E) isvector(E) && isa(E,'sym');
fVars2Check = @(vars2check) all(ismember(vars2check,'txbn'));
fOutFormCheck = @(retform) any(strcmp(retform,{'default','mat','inds'}));

ip = inputParser;
ip.addRequired('E',fEcheck);
ip.addRequired('vars',@isstruct);
ip.addRequired('refs',@isstruct);
ip.addParamValue('vars2check','txn',fVars2Check);
ip.addParamValue('outform','default',fOutFormCheck);
end