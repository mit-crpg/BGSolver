%% getrightindvects
% Right indexing vector construction function.
%
% rightinds = getrightindvects(inpids,rowids,[outform])
%
% This function constructs either one or a cell vector of right indexing
% vector(s) based on one or more given row(s) of a sparse ordered input tracking
% matrix. By default, a single right indexing vector is still enclosed in a cell
% vector. An output format command can force the right indexing vector itself to
% be returned.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Main right indexing vectors construction function code.

function rightinds = getrightindvects(inpids,rowids,outform)
% The purpose of this function is to construct either one or a cell vector of
% right indexing vector(s) based on one or more given row(s) of a sparse ordered
% input tracking matrix. By default, a single right indexing vector is still
% enclosed in a cell vector. An output format command can force the right
% indexing vector itself to be returned.
%
% INPUTS/OUTPUTS
% inpids - Sparse ordered input tracking matrix, similar in format to the one
%          output by getevids. Sparse double matrix.
% rowids - Vector of one or more row numbers of the input tracking matrix to
%          construct the right indexing vector(s) for. Horizontal uint16
%          scalar/vector.
% outform - Output format descriptor. OPTIONAL. Possible values:
%  'cell' - If there is one right indexing vector, it is enclosed in a cell
%           vector. If there are multiple, they are all contained in the cell
%           vector. This is the default option.
%  'snglvect' - If there is one right indexing vector, it is returned as is,
%               without being enclosed in a cell vector. If there are multiple,
%               they are all contained in the cell vector.
% rightinds - Either a single right indexing vector based on a given row of
%             inpids, or a vertical cell vector of right indexing vectors. The
%             output format is set by outform. Horizontal uint16 vector OR
%             vertical cell vector of horizontal uint16 vectors.

if nargin == 2
   outform = 'cell';
end;

rownum = uint16(length(rowids));
rightinds = cell(rownum,1);

for row = 1:rownum
   rowid = rowids(row);
   
   [~,nzinpIDs,inporder] = find(inpids(rowid,:));
   nzinpIDs = uint16(nzinpIDs);
   inporder = uint8(inporder);
   [~,nzinpIDsOrder] = sort(inporder);
   nzinpIDsOrder = uint8(nzinpIDsOrder);
   
   rightinds{row} = nzinpIDs(nzinpIDsOrder);
end;

if isscalar(rowids) && strcmp(outform,'snglvect')
   rightinds = rightinds{1};
end;
end