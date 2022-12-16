% EXPERIMENT FOR THE SDP FORMULATION: 

% In this file, we solve the SDP formulation of the OPF problem when the
% mpc file case3sc has a value of 60 MVA for the line-flow limit on the line from
% bus 3 to bus 2.
% The semidefinite relaxation of the OPF problem successfully solves this
% problem.
% However, solutions are different than those in the reference paper.

clearvars;
clc;

%% Load mpc data
mpc = loadcase('case3sc_60.m');

% MATPOWER results:
results_matpower = runopf(mpc);

% Create PSD matrices
[Y_k, Y_bar_k, Y_lm, Y_bar_lm, M_k] = sdp_matrices(mpc);

% Solve SDP using YALMIP 
[W_opt, cost_SDP] = SDP_formulation_primal(mpc, Y_k, Y_bar_k, Y_lm, Y_bar_lm, M_k);

% Decompose W_opt and check if V_opt holds the constraints
[V_opt_abs,V_opt_ang] = decompose_W(mpc, Y_k, Y_bar_k, Y_lm, Y_bar_lm, M_k, W_opt)
