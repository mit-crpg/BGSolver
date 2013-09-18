%% createvars
% Symbolic variable creation function.
%
% vars = createvars(bgs,[refs])
%
% This function instantiates all the symbolic variables which are used for
% sorting and formulation of AEs and returns them as symbolic vector fields of a
% variable data structure. If a referencing data structure is supplied, the bond
% variables vector is ordered according to the EBVIDs and FBVIDs calculated by
% refs, otherwise, the bond variables vector is ordered as [efforts,flows].
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Symbolic variable creation function

function vars = createvars(bgs,refs)
% The purpose of this function is to instantiate all the symbolic variables
% which are used for sorting and formulation of AEs and return them as symbolic
% vector fields of a variables data structure. If a referencing data structure
% is supplied, the bond variables vector is ordered according to the EBVIDs and
% FBVIDs calculated by refs, otherwise, the bond variables vector is ordered as
% [efforts,flows].
%
% INPUTS/OUTPUTS
% bgs - Bond Graph System data structure. EL storage type a.
% refs - Referencing arrays data structure. OPTIONAL. Important fields are:
%  .ebvid - Effort BVID calculation function (of bid).
%  .fbvid - Flow BVID calculation function (of bid).
% vars - Variables data structure. Its fields are:
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
% See the variable description cell in fBGSolve for more information.

% Instantiating unindexed variables
vars.ui.t = sym('t');
vars.ui.e = sym('e');
vars.ui.f = sym('f');
vars.ui.q = sym('q');
vars.ui.p = sym('p');
vars.ui.ei = sym('ei');
vars.ui.eo = sym('eo');
vars.ui.fi = sym('fi');
vars.ui.fo = sym('fo');

% Instantiating bond variables
% Effort and flow vectors
C.e = cell(bgs.bnum,1);
C.f = cell(bgs.bnum,1);
for bid = 1:bgs.bnum
   C.e{bid} = sprintf('e%u',bid);
   C.f{bid} = sprintf('f%u',bid);
end;
vars.e = sym(C.e);
vars.f = sym(C.f);
% Bond variable vector
if nargin == 1 % No BVID functions supplied in refs, so default ordering
   vars.b = [vars.e;vars.f];
elseif nargin == 2 % BVID functions supplied in refs
   ebvids = refs.ebvid(1:bgs.bnum);
   fbvids = refs.fbvid(1:bgs.bnum);
   efs = [vars.e;vars.f];
   efbvids = [ebvids,fbvids];
   [~,bvids] = sort(efbvids);
   vars.b = efs(bvids);
else
   error('Invalid number of inputs!!!');
end;

% Instantiating state variables
if bgs.capnum > 0
   C.q = cell(bgs.capnum,1);
   for capid = 1:bgs.capnum
      C.q{capid} = sprintf('q%u',bgs.capIDs(capid));
   end;
   vars.q = sym(C.q);
else
   vars.q = sym([]);
end;
if bgs.inertnum > 0
   C.p = cell(bgs.inertnum,1);
   for inertid = 1:bgs.inertnum
      C.p{inertid} = sprintf('p%u',bgs.inertIDs(inertid));
   end;
   vars.p = sym(C.p);
else
   vars.p = sym([]);
end;
xIDs = [bgs.capIDs,bgs.inertIDs];
[~,xIDsInds] = sort(xIDs);
vars.x = [vars.q;vars.p];
vars.x = vars.x(xIDsInds);

% Instantiating port variables
pnum = max(bgs.epnums);
if pnum > 0
   C.ep = cell(pnum,1);
   C.fp = cell(pnum,1);
   for pid = 1:pnum
      C.ep{pid} = sprintf('ep%u',pid);
      C.fp{pid} = sprintf('fp%u',pid);
   end;
   vars.ep = sym(C.ep);
   vars.fp = sym(C.fp);
else
   vars.ep = sym([]);
   vars.fp = sym([]);
end;

% Instantiating numeric variables
if bgs.nvnum > 0
   C.n = cell(bgs.nvnum,1);
   for nvid = 1:bgs.nvnum
      C.n{nvid} = sprintf('n%u',nvid);
   end;
   vars.n = sym(C.n);
else
   vars.n = sym([]);
end;
end