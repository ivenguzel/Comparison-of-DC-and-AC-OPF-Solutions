clear
close all

%% Experiment 1 (Convergence of solution for (1) AC OPF vs (2) DC OPF (Pgen*) -> AC PF)

% Output (written to a .csv file too): An array (exp1) which contains the following columns:
    % Column 1: List of load variations at Load 1
    % Column 2: List of load variations at Load 2
    % Column 3: 1 if AC OPF converged, 0 if it did not
    % Column 4: 1 if AC PF (with Pgen = DC OPF -> Pgen*) converged, 0 if it did not

% Define casefile
casefile = case9_test;

% Get base MVA value
mpc = loadcase(casefile);
base_mva = mpc.baseMVA;

% Define voltage limits for AC PF using DC OPF
voltage_ub = 1.05;
voltage_lb = 0.95;

P1_range = 0.1 : 0.1 : 5;
P2_range = 0.1 : 0.1 : 5;

disp("Number of points: ");
disp(size(P1_range,2)*size(P2_range,2))

P1_list = [];
P2_list = [];
ac_opf_success_list = [];
ac_pf_success_list = [];
 

%% Self implementation

for i = 1:size(P1_range,2)
    for j = 1:size(P2_range,2)

        % Load case file
        mpc = loadcase(casefile);
        
        % Update MPC file using above params
        Pload1_ratio = P1_range(i);
        Pload2_ratio = P2_range(j);
        
        % Scale P1 and P2 loads in MPC file with scaling factor
        mpc.bus(5,3) = mpc.bus(5,3) * Pload1_ratio;     % Load 5
        mpc.bus(7,3) = mpc.bus(7,3) * Pload2_ratio;     % Load 7

        P1_list = [P1_list ; mpc.bus(5,3)];
        P2_list = [P2_list ; mpc.bus(7,3)];
        
        % Run AC OPF
        ac_opf_matpower = opf(mpc);
        
        % Check if feasible
        ac_opf_success = ac_opf_matpower.success;
        
        % Run DC OPF if feasible
        if ac_opf_success == 1
            disp('SUCCESS');
            ac_opf_success_list = [ac_opf_success_list ; 1];
            
            [~, opt_pgen, ~] = DC_OPF_Function_Guzel_Suri(mpc);
        
            % Get Pgen values from DC OPF
            dc_opf_pgen_self = opt_pgen * base_mva;
        
            % Set DC OPF values as new Pgen values of MPC file
            mpc.gen(:,2) = dc_opf_pgen_self;
    
            %Run AC PF
            [v_all, ~, converges] = AC_PF_Function_Guzel_Suri(mpc);
        
            % Check if converges or not
            ac_pf_success = converges;
    
            for k = 1:size(v_all,1)
                if v_all(k) > voltage_ub || v_all(k) < voltage_lb
                    ac_pf_success = 0;
                end
            end
        
            if ac_pf_success == 1
                ac_pf_success_list_self = [ac_pf_success_list_self ; 1];
            else
                ac_pf_success_list_self = [ac_pf_success_list_self ; 0];
            end
        else
            ac_opf_success_list = [ac_opf_success_list ; 0];
            ac_pf_success_list_self = [ac_pf_success_list_self ; 0];
        end
    end
end

exp1_self = [P1_list P2_list ac_opf_success_list ac_pf_success_list_self];
writematrix(exp1_self,'exp1_self.csv')

%% MATPOWER implementation

for i = 1:size(P1_range,2)
    for j = 1:size(P2_range,2)

        % Load case file
        mpc = loadcase(casefile);
        
        % Update MPC file using above params
        Pload1_ratio = P1_range(i);
        Pload2_ratio = P2_range(j);
        
        % Scale P1 and P2 loads in MPC file with scaling factor
        mpc.bus(5,3) = mpc.bus(5,3) * Pload1_ratio;     % Load 5
        mpc.bus(7,3) = mpc.bus(7,3) * Pload2_ratio;     % Load 7

        P1_list = [P1_list ; mpc.bus(5,3)];
        P2_list = [P2_list ; mpc.bus(7,3)];
        
        % Run AC OPF
        ac_opf_matpower = opf(mpc);
        
        % Check if feasible
        ac_opf_success = ac_opf_matpower.success;
        
        % Run DC OPF if feasible
        if ac_opf_success == 1
            disp('SUCCESS');
            ac_opf_success_list = [ac_opf_success_list ; 1];
            
            dc_opf_matpower = dcopf(mpc);
        
            % Get Pgen values from DC OPF
            dc_opf_pgen_matpower = dc_opf_matpower.gen(:,2);
        
            % Set DC OPF values as new Pgen values of MPC file
            mpc.gen(:,2) = dc_opf_pgen_matpower;
    
            %Run AC PF
            ac_pf_matpower = runpf(mpc);
            ac_pf_voltages = ac_pf_matpower.bus(:,8);
        
            % Check if converges or not
            ac_pf_success = ac_pf_matpower.success;
    
            for k = 1:size(ac_pf_voltages,1)
                if ac_pf_voltages(k) > voltage_ub || ac_pf_voltages(k) < voltage_lb
                    ac_pf_success = 0;
                end
            end
        
            if ac_pf_success == 1
                ac_pf_success_list_matpower = [ac_pf_success_list_matpower ; 1];
            else
                ac_pf_success_list_matpower = [ac_pf_success_list_matpower ; 0];
            end
        else
            ac_opf_success_list = [ac_opf_success_list ; 0];
            ac_pf_success_list_matpower = [ac_pf_success_list ; 0];
        end
    end
end

exp1_matpower = [P1_list P2_list ac_opf_success_list ac_pf_success_list_matpower];
writematrix(exp1_matpower,'exp1_matpower.csv')
