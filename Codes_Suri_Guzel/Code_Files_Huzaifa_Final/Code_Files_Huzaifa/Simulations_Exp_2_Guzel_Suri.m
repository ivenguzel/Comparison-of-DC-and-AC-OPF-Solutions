clear
close all

%% Experiment 2 (Active Loads vs Generation for (1) AC OPF and (2) DC OPF)

%----------------------------------------------------------------------------------------------------------

% Output: Plot of Total System Loading vs Total System Generation for both AC OPF output and DC OPF output

%----------------------------------------------------------------------------------------------------------

% Define casefile
casefile = case14_test;

% Get base MVA value
mpc = loadcase(casefile);
base_mva = mpc.baseMVA;

% Define scaling factor range
scale_factor_list = 0.1 : 0.01 : 5;

%% Matpower

% Define empty arrays
P_loads_list = [];
P_gen_list_acopf = [];
P_gen_list_dcopf_matpower = [];

for i = 1 : size(scale_factor_list,2)

    % Load case file
    mpc = loadcase(casefile);
    
    scale_factor = scale_factor_list(i);
    
    % Scale Pgen, Pload in MPC file with scale factor
    mpc.gen(:,2) = mpc.gen(:,2) * scale_factor;
    mpc.bus(:,3) = mpc.bus(:,3) * scale_factor;
    
    P_loads = sum(mpc.bus(:,3));
    
    % Run AC OPF
    ac_opf_matpower = opf(mpc);
    
    % Check if feasible
    ac_opf_success = ac_opf_matpower.success;
    
    % Break loop if unfeasible
    if ac_opf_success == 0
        disp('------------------------------------------------------------')
        disp('FEASIBILITY LIMIT REACHED');
        disp('Number of data points:');
        disp(int2str(size(P_loads_list,2)))
        disp('------------------------------------------------------------')
        break;
    end
    
    % Run DC OPF if AC OPF feasible
    dc_opf_matpower = dcopf(mpc);
    
    % Compute Pgen for both AC OPF and DC OPF
    ac_opf_pgen_matpower = ac_opf_matpower.gen(:,2);
    dc_opf_pgen_matpower = dc_opf_matpower.gen(:,2);
    
    % Store generations and loads in arrays
    P_loads_list = [P_loads_list P_loads];
    P_gen_list_acopf = [P_gen_list_acopf sum(ac_opf_pgen_matpower)];
    P_gen_list_dcopf_matpower = [P_gen_list_dcopf_matpower sum(dc_opf_pgen_matpower)];

end

% Plot generations vs loads
plot(P_loads_list, P_gen_list_acopf,'b','LineWidth',1);
hold on
plot(P_loads_list, P_gen_list_dcopf_matpower,'Color','#F59522','Linewidth',1);
xlabel('Total System Loading (MW)');
ylabel('Total System Generation (MW)');
legend({'AC OPF solution','DC OPF solution'},'Location','northwest');

%% Self

% Define empty arrays
P_loads_list = [];
P_gen_list_acopf = [];
P_gen_list_dcopf_self = [];

for i = 1 : size(scale_factor_list,2)

    % Load case file
    mpc = loadcase(casefile);
    
    scale_factor = scale_factor_list(i);
    
    % Scale Pgen, Pload, Qload in MPC file with scale_factor
    mpc.gen(:,2) = mpc.gen(:,2) * scale_factor;
    mpc.bus(:,3) = mpc.bus(:,3) * scale_factor;
%     mpc.bus(:,4) = mpc.bus(:,4) * scale_factor;
    
    P_loads = sum(mpc.bus(:,3));
    
    % Run AC OPF
    ac_opf_matpower = opf(mpc);
    
    % Check if feasible
    ac_opf_success = ac_opf_matpower.success;
    
    % Break loop if unfeasible
    if ac_opf_success == 0
        disp('------------------------------------------------------------')
        disp('FEASIBILITY LIMIT REACHED');
        disp('Number of data points:');
        disp(int2str(size(P_loads_list,2)))
        disp('------------------------------------------------------------')
        break;
    end
    
    % Run DC OPF if AC OPF feasible
    [~, opt_pgen, ~] = DC_OPF_Function_Guzel_Suri(mpc);
    
    % Compute Pgen for both AC OPF and DC OPF
    ac_opf_pgen_matpower = ac_opf_matpower.gen(:,2);
    dc_opf_pgen_matpower = opt_pgen * base_mva;
    
    % Store generations and loads in arrays
    P_loads_list = [P_loads_list P_loads];
    P_gen_list_acopf = [P_gen_list_acopf sum(ac_opf_pgen_matpower)];
    P_gen_list_dcopf_self = [P_gen_list_dcopf_self sum(dc_opf_pgen_matpower)];

end

% Plot generations vs loads
plot(P_loads_list, P_gen_list_acopf,'b','LineWidth',1);
hold on
plot(P_loads_list, P_gen_list_dcopf_self,'Color','#F59522','Linewidth',1);
xlabel('Total System Loading (MW)');
ylabel('Total System Generation (MW)');
legend({'AC OPF solution','DC OPF solution'},'Location','northwest');