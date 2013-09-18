%% runti
% Time integrator execution function.
%
% [T,X] = runti(xdot,runopts,X0,tioptions)
%
% This function executes the time integrator on the state derivative vector
% according to time step type, time integrator type, integration interval,
% initial values vector and time integrator options. Prior to doing so, it
% evaluates the initial state derivative, to ensure that the referenced
% functions are compiled by the Just-in-Time compiler.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Time integrator execution function

function [T,X] = runti(xdot,runopts,X0,tioptions,varargin)
% The purpose of this function is to execute the time integrator on the state
% derivative system.
%
% INPUTS/OUTPUTS
% xdot - State derivative vector function handle.
% runopts - Runtime and Recording Options data structure.
% X0 - Initial values vector.
% tioptions - Time integrator options data structure.
% T - Vertical vector of time values.
% X - Vertical array of horizontal state vectors corresponding to the time
%     values in T.

% Checking time integrator
assert(runopts.tst ~= 2,...
   'Fixed time-stepping with adaptive time integrator not supported!');

% Preliminary evaluation
xdot0 = xdot(runopts.ii(1),X0); %#ok<NASGU>
% Clearing memory
clear('xdot0');

if any(runopts.tst == [0,1])
   [T,X] = feval(runopts.tit,xdot,runopts.ii,X0,tioptions);
elseif runopts.tst == 2
   %
   % STUB
   % Using an adaptive time integrator to take small fixed size steps.
   % Potentially runtime-recording every few time steps.
   % Here the code should involve a loop, time stepping through the fixed size
   % steps, running the adaptive time integrator from t_n to t_n+1, then
   % advancing n.
   % Write this section if ever decide to use TST = 2. Until then, stub.
   % /STUB
   %
end;

% Clearing memory
clear('evalxdot');
end