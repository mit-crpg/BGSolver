%% evalxdot1phx
% Single-physics state derivative evaluation function.
%
% xdot1phx = evalxdot1phx(t,x1phx,xdot,phxid,xotherphx,phxsys)
%
% This function evaluates the state derivative of a single physics at a given
% time, with a given (possibly vectorized) x vector for that physics, based on a
% given state derivative vector, the single physics ID, an array of state
% vectors (not vectorized) for the other physics and the physics splitting
% system data structure.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Single-physics state derivative evaluation function.

function xdot1phx = evalxdot1phx(t,x1phx,xdot,phxid,xotherphx,phxsys)
% The purpose of this function is to evaluate the state derivative of a single
% physics at a given time.
%
% INPUTS/OUTPUTS
% t - Time to input to xdot. double scalar.
% x1phx - Single physics state vector to evaluate xdot at, possibly vectorized
%         such that x1phx(:,m) is the m'th vertical single physics state vector
%         to input to xdot. Vertical double vector OR a matrix of vertical
%         double vectors.
% xdot - Full state derivative vector function handle.
% phxid - Single physics ID being evaluated.
% xotherphx - Other physics' state vectors to pass to xdot. Cell vector of
%             vertical double state vectors, sorted in ascending physics order.
%             xotherphx{phxid} is not used, and should be empty (but does not
%             have to be). The individual state vectors are ordered in ascending
%             full XID order.
% phxsys - Physics splitting system data structure.
% xdot1phx - Single physics state derivative vector evaluated at time t, with
%            single physics state vector x1phx and with other physics' state
%            values xotherphx. If x1phx was vectorized, xdot1phx is a matrix,
%            such that xdot(:,m) is the m'th vertical single physics state
%            derivative vector at time t and single physics state x1phx(:,m).

% Preallocating state vector
x = zeros(phxsys.xnum,1);

% Filling the other physics' state vectors into the state vector
ophxids = 1:phxsys.phxnum;
ophxids(phxid) = [];
for ophxid = ophxids
   x(phxsys.xphxs{ophxid}) = xotherphx{ophxid};
end;

% Sizing the problem
vectnum = uint32(size(x1phx,2));
x = repmat(x,1,vectnum);

% Retrieving single physics XIDs
xids1phx = phxsys.xphxs{phxid};

% Filling the state vector with single physics
x(xids1phx,:) = x1phx;

% Evaluating the full state
xdotfull = xdot(t,x);

% Trimming the full state
xdot1phx = xdotfull(xids1phx,:);
end