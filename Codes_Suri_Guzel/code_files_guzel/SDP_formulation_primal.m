function [W_opt, cost_SDP] = SDP_formulation_primal(mpc, Y_k, Y_bar_k, Y_lm, Y_bar_lm, M_k)
%% Primal SDP Relaxation Problem (Optimization 3 in Lavaei,Low)

N_bus = size(mpc.bus,1);
N_gen = size(mpc.gen,1);
N_branch = size(mpc.branch,1);
Sbase = mpc.baseMVA;

% Clear YALMIPs internal database
yalmip('clear')

% Define variables
W = sdpvar(2*N_bus,2*N_bus, 'symmetric');    % Symmetric 2N_bus matrix
alpha = sdpvar(N_gen, 1);  % Size N_gen vector: From the schur complement of generation cost C(P_k) <= alpha_k

% Define objective
Objective = sum(alpha); % C(P_k) <= alpha_k

Constraints = [W >= 0];

% Define constraints:

for i=1:N_bus
    % BUS DATA COLUMNS: BUS_I:1, PD:3, QD:4, VMIN: 13, VMAX: 12, BUS_TYPE:2 (1 = PQ, 2 = PV, 3 = ref, 4 = isolated)
    % GEN DATA COLUMNS: GEN_BUS: 1, PMAX:9, PMIN:10, QMAX:4, QMIN:5
    % GEN COST DATA: Only MODEL = 2 works with this code. 

    % Bus Voltage Magnitudes
    Vmin = mpc.bus(i,13);
    Vmax = mpc.bus(i,12);
    Constraints = [Constraints, W >= 0,(Vmin)^2 <= trace(M_k{i}*W) ];
    Constraints = [Constraints, W >= 0, (Vmax)^2 >= trace(M_k{i}*W)];

    % P and Q injections 
    Pd = mpc.bus(i,3) / Sbase;
    Qd = mpc.bus(i,4) / Sbase;
    if mpc.bus(i,2) == 1 % Check PQ bus
        Constraints = [Constraints, W >= 0, -Pd <= trace(Y_k{i}*W)];
        Constraints = [Constraints, W >= 0, -Pd >= trace(Y_k{i}*W)];
        Constraints = [Constraints, W >= 0, -Qd <= trace(Y_bar_k{i}*W)];
        Constraints = [Constraints, -Qd >= trace(Y_bar_k{i}*W)];
    end
    if mpc.bus(i,2) == 2 || mpc.bus(i,2) == 3  % Check PV or slack 
        gen_ind = find(mpc.gen(:,1) == i);
        Pmin = mpc.gen(gen_ind,10) / Sbase;
        Pmax = mpc.gen(gen_ind,9) / Sbase;
        Qmin = mpc.gen(gen_ind,5) / Sbase;
        Qmax = mpc.gen(gen_ind,4) / Sbase;
        % P and Q injection constraints for generators
            if abs(Pmin-Pmax) <= 0.1 
                Constraints = [Constraints, -Pd == trace(Y_k{i}*W)];
            else
                Constraints = [Constraints, Pmin-Pd <= trace(Y_k{i}*W)];
                Constraints = [Constraints, Pmax-Pd >= trace(Y_k{i}*W)];
            end
        Constraints = [Constraints, Qmin-Qd <= trace(Y_bar_k{i}*W)];
        Constraints = [Constraints, Qmax-Qd >= trace(Y_bar_k{i}*W)];
        % Schur complements of generator costs
        c = [mpc.gencost(gen_ind,5)*Sbase^2, mpc.gencost(gen_ind,6)*Sbase, mpc.gencost(gen_ind,7)];
        Constraints = [Constraints, [c(2)*trace(Y_k{i}*W) - alpha(i) + c(3) + c(2)*Pd, sqrt(c(1))*trace(Y_k{i}*W) + sqrt(c(1))*Pd ;...
                                    sqrt(c(1))*trace(Y_k{i}*W) + sqrt(c(1))*Pd,-1] <= 0];
    end
end

for i = 1:N_branch
    % BRANCH DATA COLUMNS: RATE_A:6
    Smax = mpc.branch(i,6) / Sbase;
    Constraints = [Constraints, trace(Y_lm{i}*W) <= Smax];
    Constraints = [Constraints, [-Smax^2, trace(Y_lm{i}*W), trace(Y_bar_lm{i}*W); trace(Y_lm{i}*W), -1, 0; trace(Y_bar_lm{i}*W),0,-1] <= 0];
end

% Set options for YALMIP and solver
options =  sdpsettings('solver','mosek');

sol = optimize(Constraints,Objective,options);

W_opt = value(W);
cost_SDP = value(Objective);

end