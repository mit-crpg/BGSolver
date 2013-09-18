%% createrefs
% Reference arrays creation function.
%
% refs = createrefs(bgs)
%
% This function creates a referencing array data structure, for referencing
% quantities like expressions by EID (or BID), and not by the quantity ID.
% Referenced quantities are:
% - Expressions (for elements with expressions).
% - Expression types (for elements with expressions).
% - Modulating variable lists (for modulated elements).
% - Port numbers (for elements with ambiguous numbers of ports).
% - Bond-to-port connectivities (for elements with ambiguous bond-to-port
%   connectivities).
% - Causalities (for elements with ambiguous causalities).
% - Physics (for physics-specific elements).
% - Non-negativity statuses (for elements whose state variables can be enforced
%   to be non-negative).
% - Absolute tolerances (for elements whose state variables can have specific
%   absolute tolerances assigned to them).
% - Numeric expressions (for numeric elements).
% - State variables' XIDs (state vector indices).
% - Effort variables' BVIDs.
% - Flow varaibles' BVIDs.
% - BIDs based on Effort BVIDs.
% - BIDs based on Flow BVIDs.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Master referencing arrays creation function

function refs = createrefs(bgs)
% The purpose of this function is to create a referencing array data structure,
% for referencing quantities like expressions by EID (or BID), and not by the
% quantity ID. Referenced quantities are:
% - Expressions (for elements with expressions).
% - Expression types (for elements with expressions).
% - Modulating variable lists (for modulated elements).
% - Port numbers (for elements with ambiguous numbers of ports).
% - Bond-to-port connectivities (for elements with ambiguous bond-to-port
%   connectivities).
% - Causalities (for elements with ambiguous causalities).
% - Physics (for physics-specific elements).
% - Non-negativity statuses (for elements whose state variables can be enforced
%   to be non-negative).
% - Absolute tolerances (for elements whose state variables can have specific
%   absolute tolerances assigned to them).
% - Numeric expressions (for numeric elements).
% - State variables' XIDs (state vector indices).
% - Effort variables' BVIDs.
% - Flow variables' BVIDs.
% - BIDs based on Effort BVIDs.
% - BIDs based on Flow BVIDs.
%
% INPUTS/OUTPUTS
% bgs - Bond Graph System data structure. EL storage type a.
% refs - Referencing arrays data structure. Its fields are:
%  .exprs - Expression referencing array.
%  .mvars - Modulating variable list referencing array.
%  .pnums - Port number referencing array.
%  .bcons - Bond-to-port connectivity referencing array.
%  .causes - Causality referencing array.
%  .phxs - Physics referencing array.
%  .nns - Non-negativity statuses referencing array.
%  .ats - Absolute tolerances referencing array.
%  .nums - Numeric element referencing array.
%  .xids - State variables' state vector index array.
%  .ebvid - Effort BVID calculation function (of bid).
%  .fbvid - Flow BVID calculation function (of bid).
%  .bidebv - BID based on Effort BVID calculation function.
%  .bidfbv - BID based on Flow BVID calculation function.
% See the variable description cell in fBGSolve for more information.

% Expression referencing array
refs.exprs = zeros(1,bgs.enum,'uint16');
refs.exprs(bgs.exprIDs) = 1:bgs.exprnum;

% Modulating variable list referencing array
refs.mvars = zeros(1,bgs.enum,'uint16');
refs.mvars(bgs.modIDs) = 1:bgs.emodnum;

% Port number referencing array
refs.pnums = zeros(1,bgs.enum,'uint16');
refs.pnums(bgs.pambIDs) = 1:bgs.epambnum;

% Bond-to-port connectivity referencing array
refs.bcons = zeros(1,bgs.enum,'uint16');
refs.bcons(bgs.bambIDs) = 1:bgs.ebambnum;

% Causality referencing array
refs.causes = zeros(1,bgs.enum,'uint16');
refs.causes(bgs.cambIDs) = 1:bgs.ecambnum;

% Physics referencing array
refs.phxs = zeros(1,bgs.enum,'uint16');
refs.phxs(bgs.phxIDs) = 1:bgs.phxspecnum;

% Non-negativity statuses referencing array
refs.nns = zeros(1,bgs.enum,'uint16');
refs.nns(bgs.nnIDs) = 1:bgs.nnspecnum;

% Absolute tolerances referencing array
refs.ats = zeros(1,bgs.enum,'uint16');
refs.ats(bgs.atIDs) = 1:bgs.atspecnum;

% Numeric element referencing array
refs.nums = zeros(1,bgs.enum,'uint16');
refs.nums(bgs.numIDs) = 1:bgs.ennum;

% State variables' state vector index referencing array
refs.xids = zeros(1,bgs.enum,'uint16');
refs.xids(bgs.x0IDs) = 1:bgs.xnum;

% Effort BVID calculation function.
refs.ebvid = @(bid) 2*bid - 1;

% Flow BVID calculation function.
refs.fbvid = @(bid) 2*bid;

% BID based on Effort BVID calculation function.
refs.bidebv = @(ebvid) (ebvid+1) / 2;

% BID based on Flow BVID calculation function.
refs.bidfbv = @(fbvid) fbvid / 2;
end