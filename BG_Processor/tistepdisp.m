%% tistepdisp
% Time step information display function for use with MATLAB built-in time
% integrators.
%
% status = tistepdisp(t,x,[flag])
%
% This function works similarly to MATLAB built-in odeprint, without clearing
% the screen. It is called by the MATLAB built-in adaptive time integrator prior
% to starting the integration, every time it takes a time step, and after all
% time steps have been taken. The evaluation type is specified by the flag
% argument. The function outputs the information about the time step, and is
% intended to be used by enabling the 'OutputFcn' field in time integrator
% options data structure.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Time integrator time step information display function

function status = tistepdisp(t,x,flag) %#ok<INUSL>
% The purpose of this function is to output the information about a time step
% the time integrator takes.
%
% INPUTS/OUTPUTS
% t - During the initialization, t is the [t_0,t_f] vector. During subsequent
%     time steps, t is the time to which the time step was just taken. After the
%     final time step, t is empty.
% x - During the initialization, x is the initial state vector. During
%     subsequent time steps, x is the state vector to which the time step was
%     just taken. After final time step, x is empty.
% flag - Evaluation type of the time step. OPTIONAL. Can be one of the following
%        options:
%  'init' - Initialization. Called prior to taking any time steps.
%  [] (or not given) - Regular time step. Called every time a time step is
%                      taken. This is the default option.
%  'done' - Completion. Called after the final time step.
% status - Evaluation status. If the evaluation is successful, status is 0.
%          Otherwise, it is non-zero. Currently, it is always set to 0.

if nargin == 2 || isempty(flag)
   fprintf('Integrating, t = %19.12g\n',t);
else
   if strcmp(flag,'init')
      fprintf('Initializing time integrator...\n');
      fprintf('Initial time = %19.12g\n',t(1));
      fprintf('Final time   = %19.12g\n',t(2));
   elseif strcmp(flag,'done')
      fprintf('Time integration completed.\n');
   end;
end;

status = 0;
end