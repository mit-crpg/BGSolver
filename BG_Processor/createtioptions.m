%% createtioptions
% Time integrator options creation function.
%
% tioptions = createtioptions(srtsys,runopts,fids,['Vectorized',Vectorized])
%
% This function creates the custom time integration options data structure for
% use with both custom-written and MATLAB built-in time integrators.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Variable descriptions
%
% The time integrator options (tioptions) data structure is passed to both
% custom-written and MATLAB built-in time integrators. As such, it has all the
% fields that a MATLAB ODE options data structure has; the non-default ones are:
% .Vectorized - Set to 'on', unless overwritten by input.
% 
% The custom fields are:
% .tdisc - Time discretization options data structure.
% .sjso - Solver and jumpstart options data structure.
% .phxsys - Physics splitting options data structure.
% .fids - Open file IDs data structure.
% .rtrr - Runtime results recording options data structure.

%% Time integrator options creation function

function tioptions = createtioptions(varargin)
% The purpose of this function is to create the time integrator options
% data structure, based on several run options and the physics splitting options
% in the sorting system.
%
% INPUTS/OUTPUTS
% srtsys - Sorting system data structure.
% runopts - Runtime and Recording Options data structure.
% fids - File IDs data structure.
% Vectorized - State derivative vector function state vectorization option.
% tioptions - Time Integrator options data structure.

% Input parser construction
ip = buildIP();
% Input parsing
ip.parse(varargin{:});
% Retrieving inputs
srtsys = ip.Results.srtsys;
runopts = ip.Results.runopts;
fids = ip.Results.fids;
Vectorized = ip.Results.Vectorized;
% Clearing memory
clear('ip');

% Creating standard MATLAB ODE options data structure
% Assuming Jacobian sparsity pattern has already been constructed
if runopts.tidisp == 0
   tioptions = odeset('Vectorized',Vectorized,'NonNegative',srtsys.xnns,...
      'AbsTol',srtsys.xats,'JPattern',srtsys.xdotjsp);
elseif runopts.tidisp == 1
   tioptions = odeset('Vectorized',Vectorized,'NonNegative',srtsys.xnns,...
      'AbsTol',srtsys.xats,'JPattern',srtsys.xdotjsp,'OutputFcn',@tistepdisp);
end;
for addtioID = 1:runopts.addtio
   pname = runopts.addtiolist{addtioID,1};
   pval = runopts.addtiolist{addtioID,2};
   tioptions.(pname) = pval;
end;

% Appending time discretization options
tioptions.tdisc.tst = runopts.tst;
if any(runopts.tst == [1,2])
   tioptions.tdisc.dt = runopts.fts;
end;

% Appending solver and jumpstart options
tioptions.sjso.stolx = runopts.stolx;
tioptions.sjso.stolres = runopts.stolres;
tioptions.sjso.sdisp = runopts.sdisp;
tioptions.sjso.srscal = runopts.srscal;
tioptions.sjso.sxscal = runopts.sxscal;
tioptions.sjso.sxcent = runopts.sxcent;
tioptions.sjso.spc = runopts.spc;
tioptions.sjso.sigpdct = runopts.sigpdct;
tioptions.sjso.sjacob = runopts.sjacob;
tioptions.sjso.sjsp = runopts.sjsp;
tioptions.sjso.addso = runopts.addso;
if runopts.addso > 0
   tioptions.sjso.addsolist = runopts.addsolist;
end;
tioptions.sjso.jstype = runopts.jstype;
tioptions.sjso.jsmsnum = runopts.jsmsnum;
tioptions.sjso.jsphxsplit = runopts.jsphxsplit;
tioptions.sjso.jstolx = runopts.jstolx;
tioptions.sjso.jstolres = runopts.jstolres;
tioptions.sjso.jstit = runopts.jstit;
tioptions.sjso.jsdisp = runopts.jsdisp;
tioptions.sjso.jsrscal = runopts.jsrscal;
tioptions.sjso.jsxscal = runopts.jsxscal;
tioptions.sjso.jsxcent = runopts.jsxcent;
tioptions.sjso.jspc = runopts.jspc;
tioptions.sjso.jsjacob = runopts.jsjacob;
tioptions.sjso.jsjsp = runopts.jsjsp;
tioptions.sjso.addjso = runopts.addjso;
if runopts.addjso > 0
   tioptions.sjso.addjsolist = runopts.addjsolist;
end;

% Appending physics splitting options
tioptions.phxsys = createphxsys(srtsys,runopts);

% Appending open file options
tioptions.fids = fids;

% Appending runtime results recording options
tioptions.rtrr.ifrtrr = runopts.rtrr;
if tioptions.rtrr.ifrtrr && isfield(runopts,'rtrrfq')
   tioptions.rtrr.rtrrfq = runopts.rtrrfq;
end;
end

%% Input parser construction function

function ip = buildIP()
% The purpose of this function is to construct an input parser object for the
% createtioptions function.
%
% INPUTS/OUTPUTS
% ip - Input parser object for the createtioptions function.

vchkfun = @(Vectorized) any(strcmp(Vectorized,{'on','off'}));

ip = inputParser;
ip.StructExpand = false;
ip.addRequired('srtsys',@isstruct);
ip.addRequired('runopts',@isstruct);
ip.addRequired('fids',@isstruct);
ip.addParamValue('Vectorized','on',vchkfun);
end