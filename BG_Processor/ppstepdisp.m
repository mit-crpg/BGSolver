%% ppstepdisp
% Time step information display function for use with post-processor.
%
% status = ppstepdisp(tn)
%
% This function outputs the information about the time step after every time
% step taken by the post-processor.
% Currently, the function only outputs the time. It may, in the future, be
% configured to output additional information.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Post-processor time step information display function

function status = ppstepdisp(tn)
% The purpose of this function is to output the information about a time step
% after every time step the post-processor takes.
%
% INPUTS/OUTPUTS
% tn - Time of the last post-processing evaluation.
% status - Evaluation status. If the evaluation is successful, status is 0.
%          Otherwise, it is non-zero. Currently, it is always set to 0.

fprintf('Post-processing, t = %19.12g\n',tn);

status = 0;
end