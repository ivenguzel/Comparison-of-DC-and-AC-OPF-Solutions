%----------------------------------------------------------------------------------------------

% Input:
    % casefile: Network in MATPOWER Casefile Format
% Output:
    % opt_var: Set of optimal variables [PG1*  PG2* ... PG_N*  theta1*  theta2* ... theta_N*]
    % opt_gen: Set of optimal generations at each bus
    % lmp: Set of locational marginal prices (LMPs) for all buses

% Bus Types:
    % PQ bus = 1
    % PV bus = 2
    % reference bus = 3
    % isolated bus = 4

%----------------------------------------------------------------------------------------------

function [opt_var, opt_pgen, lmp] = DC_OPF_Function_Guzel_Suri(casefile)
    
    % Load Case File
    mpc = loadcase(casefile);
    
    % Base MVA
    base_mva = mpc.baseMVA;
    
    % Bus Data
    bus_no = mpc.bus(:,1);
    bus_type = mpc.bus(:,2);
    P_demand = mpc.bus(:,3)/base_mva;
    Q_demand = mpc.bus(:,4)/base_mva;
    G_shunt = mpc.bus(:,5);
    B_shunt = mpc.bus(:,6);
    area_no = mpc.bus(:,7);
    Vm_bus = mpc.bus(:,8);
    theta_bus = mpc.bus(:,9);
    base_voltage = mpc.bus(:,10);
    loss_zone = mpc.bus(:,11);
    Vm_max = mpc.bus(:,12);
    Vm_min = mpc.bus(:,13);
    
    % Generator Data
    gen_bus_no = mpc.gen(:,1);
    gen_P_gen = mpc.gen(:,2)/base_mva;
    gen_Q_gen = mpc.gen(:,3)/base_mva;
    gen_Q_max = mpc.gen(:,4)/base_mva;
    gen_Q_min = mpc.gen(:,5)/base_mva;
    gen_Vg = mpc.gen(:,6);
    gen_mBase = mpc.gen(:,7);
    gen_status = mpc.gen(:,8);
    gen_Pmax = mpc.gen(:,9)/base_mva;
    gen_Pmin = mpc.gen(:,10)/base_mva;
    gen_Pc1 = mpc.gen(:,11);
    gen_Pc2 = mpc.gen(:,12);
    gen_Qc1min = mpc.gen(:,13);
    gen_Qc1max = mpc.gen(:,14);
    gen_Qc2min = mpc.gen(:,15);
    gen_Qc2max = mpc.gen(:,16);
    
    % Branch Data
    branch_from = mpc.branch(:,1);
    branch_to = mpc.branch(:,2);
    branch_res = mpc.branch(:,3);
    branch_x = mpc.branch(:,4);
    branch_b = mpc.branch(:,5);
    branch_rateA = mpc.branch(:,6)/base_mva;
    branch_rateB = mpc.branch(:,7)/base_mva;
    branch_rateC = mpc.branch(:,8)/base_mva;
    branch_transformer_ratio = mpc.branch(:,9);
    branch_transformer_angle = mpc.branch(:,10);
    branch_status_init = mpc.branch(:,11);
    branch_min_angle_diff = mpc.branch(:,12);
    branch_max_angle_diff = mpc.branch(:,13);
    
    % Generator Cost Data
    gencost_model = mpc.gencost(:,1);
    gencost_startup_cost = mpc.gencost(:,2);
    gencost_N = mpc.gencost(:,3);
    gencost_no_of_coeffs = mpc.gencost(:,4);
    gencost_params_c2 = mpc.gencost(:,5);
    gencost_params_c1 = mpc.gencost(:,6);
    gencost_params_c0 = mpc.gencost(:,7);
    
    
    generators = size(gen_bus_no,1);
    number_of_buses = size(bus_no,1);
    variables = generators + number_of_buses;
    
    %% Constructing Y-bus
    
    Y = zeros(number_of_buses);
    i = sqrt(-1);
    
    for j = 1:size(branch_from)
    
        Y(branch_from(j),branch_to(j)) = Y(branch_from(j),branch_to(j)) - 1 / (branch_res(j) + i*branch_x(j));
        Y(branch_to(j),branch_from(j)) = Y(branch_to(j),branch_from(j)) - 1 / (branch_res(j) + i*branch_x(j));
        Y(branch_from(j),branch_from(j)) = Y(branch_from(j),branch_from(j)) + (1 / (branch_res(j) + i*branch_x(j))) + i*branch_b(j)/2;
        Y(branch_to(j),branch_to(j)) = Y(branch_to(j),branch_to(j)) + (1 / (branch_res(j) + i*branch_x(j))) + i*branch_b(j)/2;
    
    end
    b_bus = imag(Y);
    
    %% Generator Constraints (upper lower bounds)
    
    gen_constraint_lower = zeros(generators,variables);
    
    for i = 1:generators
        gen_constraint_lower(i,i) = -1;
    end
    
    gen_constraint_upper = zeros(generators,variables);
    
    for i = 1:generators
        gen_constraint_upper(i,i) = 1;
    end
    
    %% Supply/Demand Constraints (b = 0)
    
    sup_dem_cons = zeros(number_of_buses,variables);
    
    for i = 1:number_of_buses
    
        for j = 1:number_of_buses
    
            % check if branch exists b/w buses
    
            if(isempty(intersect(find(branch_from == i),find(branch_to == j))) == 0) 
    
                sup_dem_cons(i,variables-number_of_buses + i) = sup_dem_cons(i,variables-number_of_buses + i) + b_bus(i,j);
                sup_dem_cons(i,variables-number_of_buses + j) = sup_dem_cons(i,variables-number_of_buses + j) - b_bus(i,j);
    
            elseif(isempty(intersect(find(branch_from == j),find(branch_to == i))) == 0) 
    
                sup_dem_cons(i,variables-number_of_buses + i) = sup_dem_cons(i,variables-number_of_buses + i) + b_bus(i,j);
                sup_dem_cons(i,variables-number_of_buses + j) = sup_dem_cons(i,variables-number_of_buses + j) - b_bus(i,j);
            
            end
        end
    end
    
    for i = 1:generators
        sup_dem_cons(gen_bus_no(i),i) = -1;
    end
    
    %% Line Flow Constraints
    
    line_cons_from = zeros(size(branch_rateA,1),variables);
    line_cons_to = zeros(size(branch_rateA,1),variables);
    
    for i = 1:size(branch_rateA)
        b_branch_from = b_bus(branch_from(i),branch_to(i));
        b_branch_to = b_bus(branch_to(i),branch_from(i));
        
        if branch_rateA(i) ~= 0
            line_cons_from(i, variables - number_of_buses + branch_from(i)) = b_branch_from; 
            line_cons_from(i, variables - number_of_buses + branch_to(i)) = -b_branch_from;
            line_cons_to(i, variables - number_of_buses + branch_to(i)) = b_branch_to; 
            line_cons_to(i, variables - number_of_buses + branch_from(i)) = -b_branch_to;
        end
    
    end
    
    %% Slack Bus Constraints
    
    slack_buses = length(find(bus_type == 3));
    slack_cons = zeros(slack_buses,variables);
    slack_cons(variables-number_of_buses + find(bus_type == 3)) = 1;
    slack_cons_b = zeros(slack_buses,1);
    
    %% Objective Function (optimum Pgen and angles are in p.u. and radians respectively)
    
    alpha = diag(gencost_params_c2);
    beta = gencost_params_c1;
    gamma = gencost_params_c0;
    
    cvx_begin
    
        variable x(variables,1)     % [PG1*  PG2* ... PG_N*  theta1*  theta2* ... theta_N*]
        dual variable y
        minimize( x(1:variables-number_of_buses)' * alpha * x(1:variables-number_of_buses) + beta' * x(1:variables-number_of_buses) +  ones(size(gamma,1),1)' * gamma )
        
        subject to
    
            gen_constraint_upper * x <= gen_Pmax;           % Generation upper bound constraint
            gen_constraint_lower * x <= -gen_Pmin;          % Generation lower bound constraint
            slack_cons * x == slack_cons_b;                 % Slack bus angle constraint
            y: sup_dem_cons * x == -P_demand;               % Supply/Demand constraint
            line_cons_from * x <= branch_rateA;             % Line flow constraint
            line_cons_to * x <= branch_rateA;               % Line flow constraint
    
    cvx_end
    
    %% Calculation of variables to return from function

    opt_var = x;
    opt_pgen = x(1:variables - number_of_buses);
    for i = variables - number_of_buses + 1: variables
        opt_var(i) = rad2deg(opt_var(i));
    end
    lmp = -1*y;

end