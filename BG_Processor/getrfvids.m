%% getrfvids
% Relative (in-RFX/RFPP) VID calculation function.
%
% [vidsrf,[vidsrfBool]] = getrfvids(vids,reqvids,absvnum,[ornt])
%
% This function calculates the relative (in-RFX/RFPP) variable IDs of a subset
% of RF VIDs. Optionally, it also returns the relative logical indexing variable
% ID vector. By default, the variable IDs are returned in horizontal vectors,
% but can optionally be returned as vertical vectors.
% See createxdotsys() variable descriptions cell for the definition of RFX and
% createppsys() variable descriptions cell for the definition of RFPP.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Master relative VID calculation function code

function [vidsrf,vidsrfBool] = getrfvids(vids,reqvids,absvnum,ornt)
% The purpose of this function is to calculate the relative (in-RFX/RFPP)
% variable IDs of a subset of RF VIDs. Optionally, it also returns the relative
% logical indexing variable ID vector. By default, the variable IDs are returned
% in horizontal vectors, but can optionally be returned as vertical vectors. The
% function does using temporary logical arrays, which, as shown by tests, is
% significantly faster than intersect()-based equivalent VID calculation.
%
% INPUTS/OUTPUTS
% vids - Vector of absolute variable IDs to calculate the relative VIDs of.
%        Horizontal uint16 vector OR horizontal boolean (logical indexing)
%        vector.
% reqvids - Vector of all absolute RF variable IDs, out which vids is a subset.
%           Horizontal uint16 vector OR horizontal boolean (logical indexing)
%           vector.
% absvnum - Number of absolute variable IDs. uint16 scalar.
% ornt - Output VID vector orientation. OPTIONAL. Possible values:
%  'horiz' or 'horizontal' - Output VID vector(s) is/are horizontal. This is the
%                            default option.
%  'vert' or 'vertical' - Output VID vector(s) is/are vertical.
% vidsrf - Vector of relative VIDs which correspond to vids' positions in
%          reqvids. Horizontal uint16 vector.
% vidsrfBool - Logical relative indexing vector of vidsrf. Horizontal boolean
%          vector. OPTIONAL.

% Reordering vids to be in ascending order (if numeric)
if isnumeric(vids)
   [vids,vidsOrder] = sort(vids);
   [~,vidsRevOrder] = sort(vidsOrder);
   vidsRevOrder = uint16(vidsRevOrder);
end;

vidsBool = false(1,absvnum);
vidsBool(vids) = true;
vidsrf = uint16(find(vidsBool(reqvids)));

% Reordering vids back
if isnumeric(vids)
   vidsrf = vidsrf(vidsRevOrder);
end;

if nargout == 2
   vidsrfBool = false(1,nnz(reqvids));
   vidsrfBool(vidsrf) = true;
end;

if nargin == 3 || (nargin == 4 && any(strcmp(ornt,{'horiz','horizontal'})))
	% VIDs are already horizontal, do nothing.
elseif nargin == 4 && any(strcmp(ornt,{'vert','vertical'}))
   vidsrf = vidsrf.';
   if nargout == 2
      vidsrfBool = vidsrfBool.';
   end;
end;
end