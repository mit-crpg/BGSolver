%% preallocEL
% Element List data structure preallocation function.
%
% elements = preallocEL(enum)
%
% This function accepts the element number, and constructs an array of element
% data structures with all of the appropriate fields. This array can be used as
% a field in the Bond Graph System (bgs) data structure. It does not instantiate
% any of the fields' values.
%
% Package:    BGSolver v1.03
% Subpackage: BGSD_Generator
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Element list data structure preallocation function.

function elements = preallocEL(enum)
% The purpose of this function is to take the element number and construct
% an empty element data structure array with all of the appropriate fields,
% without instantiating any of them.
%
% INPUTS/OUTPUTS
% enum - Element number.
% elements - Preallocated element list data structure.

elements(1:enum) = struct('etype',[],'exprtype',[],'x0',[],'epnum',[],...
                           'ebcon',[],'caus',[],'phx',[],'nn',[],'at',[],...
                           'mvars',[],'expr',[]);
elements = elements.';
end