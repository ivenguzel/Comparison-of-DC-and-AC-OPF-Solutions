%----------------------------------------------------------------------------

% Input:
    % casefile: Network in MATPOWER Casefile Format
% Output:
    % v_all: Set of all bus voltage magnitudes in ascending order of buses
    % theta_all: Set of all bus voltage angles in ascending order of buses
    % converges: 1 if AC PF converges, 0 if it does not converge

% Bus Types:
    % PQ bus = 1
    % PV bus = 2
    % reference bus = 3
    % isolated bus = 4

%----------------------------------------------------------------------------

function [v_all, theta_all, converges] = AC_PF_Function_Guzel_Suri(casefile)

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
    gencost_params_c2 = 0;
    gencost_params_c1 = mpc.gencost(:,5);
    gencost_params_c0 = mpc.gencost(:,6);
    
    
    generators = size(gen_bus_no,1);
    number_of_buses = size(bus_no,1);
    variables = generators + number_of_buses;
    
    %% Constructing Y-bus
    
    Y = zeros(number_of_buses);
    iota = sqrt(-1);
    
    for s = 1:size(branch_from)
    
        Y(branch_from(s),branch_to(s)) = Y(branch_from(s),branch_to(s)) - 1 / (branch_res(s) + iota*branch_x(s));
        Y(branch_to(s),branch_from(s)) = Y(branch_to(s),branch_from(s)) - 1 / (branch_res(s) + iota*branch_x(s));
        Y(branch_from(s),branch_from(s)) = Y(branch_from(s),branch_from(s)) + (1 / (branch_res(s) + iota*branch_x(s))) + iota*branch_b(s)/2;
        Y(branch_to(s),branch_to(s)) = Y(branch_to(s),branch_to(s)) + (1 / (branch_res(s) + iota*branch_x(s))) + iota*branch_b(s)/2;
    
    end
    
    %% Declare variables and known quantities
    
    pq_buses = sum(bus_type == 1);
    pv_buses = sum(bus_type == 2);
    slack_buses = sum(bus_type == 3);
    
    % Arrays containing all variables (known + unknown)
    % (NO NOT CHANGE)
    p_all = zeros(number_of_buses,1);
    q_all = zeros(number_of_buses,1);
    v_all = ones(number_of_buses,1);
    theta_all = zeros(number_of_buses,1);
    % theta_all = [5; 15; 10; 7; 9];
    
    % Arrays containing only variables that will be computed using NR
    % (NO NOT CHANGE)
    p_var = zeros(number_of_buses - slack_buses,1);
    q_var = zeros(number_of_buses - slack_buses - pv_buses,1);
    v_var = zeros(number_of_buses - slack_buses - pv_buses,1);
    theta_var = zeros(number_of_buses - slack_buses,1);
    
    % Arrays that tell whether certain variables need to be computed (1 = compute)
    % (NO NOT CHANGE)
    p_compute = ones(number_of_buses,1);
    q_compute = ones(number_of_buses,1);
    v_compute = ones(number_of_buses,1);
    theta_compute = ones(number_of_buses,1);
    
    for v = 1:number_of_buses
    
        if bus_type(v) == 2      %PV
            v_compute(v) = 0;
            q_compute(v) = 0;
        elseif bus_type(v) == 3  %slack
            v_compute(v) = 0;
            theta_compute(v) = 0;
            p_compute(v) = 0;
            q_compute(v) = 0;
        end
    
    end
    
    % Arrays with indices of variables that have to be computed
    p_compute_indices = find(p_compute == 1);
    q_compute_indices = find(q_compute == 1);
    v_compute_indices = find(v_compute == 1);
    theta_compute_indices = find(theta_compute == 1);
   
    
    % Initialize variables

    p_gen_opt = gen_P_gen;

    for z = 1:generators
        p_all(gen_bus_no(z)) = p_all(gen_bus_no(z)) + p_gen_opt(z);
        q_all(gen_bus_no(z)) = q_all(gen_bus_no(z)) + gen_Q_gen(z);
    end
    p_all = p_all - P_demand;
    q_all = q_all - Q_demand;     % v and theta already initialized as 1 and 0 respectively
    
    
    % Extract unknown variables from all variables
    p_counter = 1;
    q_counter = 1;
    v_counter = 1;
    theta_counter = 1;
    
    for n = 1:number_of_buses
    
        if p_compute(n) == 1
            p_var(p_counter) = p_all(n);
            p_counter = p_counter + 1;
        end
        if q_compute(n) == 1
            q_var(q_counter) = q_all(n);
            q_counter = q_counter + 1;
        end
        if v_compute(n) == 1
            v_var(v_counter) = v_all(n);
            v_counter = v_counter + 1;
        end
        if theta_compute(n) == 1
            theta_var(theta_counter) = theta_all(n);
            theta_counter = theta_counter + 1;
        end
    
    end
    
    %% Compute initial F = [del_p ; del_q] (p - P, q - Q)
    
    del_p = zeros(size(p_var,1),1);
    p_counter = 1;
    for l = 1:size(p_all,1)
        if(p_compute(l) == 1)
            del_p(p_counter) = acpf_p_eq(l,number_of_buses,v_all,theta_all,Y);
            p_counter = p_counter + 1;
        end
    end
    
    del_p = del_p - p_var;         % del_P = p(v,theta) - P
    
    del_q = zeros(size(q_var,1),1);
    q_counter = 1;
    for l = 1:size(q_all,1)
        if(q_compute(l) == 1)
            del_q(q_counter) = acpf_q_eq(l,number_of_buses,v_all,theta_all,Y);
            q_counter = q_counter + 1;
        end
    end
    
    del_q = del_q - q_var;        % del_Q = q(v,theta) - Q
    
    F = [del_p ; del_q];
    
    
    %% Run NR for (max_iter) iterations
     
    tolerance = 1e-5;
    
    iterate = true;
    
    iter = 0;
    
    max_iter = 1000;
    
    converges = 0;
    
    while (iterate)
        
        % Compute F (P/Q known - p/q predicted)
        
        del_p = zeros(size(p_var,1),1);
        p_counter = 1;
        for m = 1:size(p_all,1)
            if(p_compute(m) == 1)
                del_p(p_counter) = acpf_p_eq(m,number_of_buses,v_all,theta_all,Y);
                p_counter = p_counter + 1;
            end
        end
        del_p = del_p - p_var;         % del_P = p(v,theta) - P
        
        del_q = zeros(size(q_var,1),1);
        q_counter = 1;
        for m = 1:size(q_all,1)
            if(q_compute(m) == 1)
                del_q(q_counter) = acpf_q_eq(m,number_of_buses,v_all,theta_all,Y);
                q_counter = q_counter + 1;
            end
        end
        del_q = del_q - q_var;        % del_Q = q(v,theta) - Q
        
        F = [del_p ; del_q];
    
        % Check if error is within tolerance
        if norm(F,2) > tolerance
    
            % Compute Jacobian
            J = compute_jacobian(Y,number_of_buses,v_all,theta_all,p_compute_indices,q_compute_indices,v_compute_indices,theta_compute_indices);
    
            % Assert Jacobian to be non singular
            assert(abs(det(J)) > 1e-5)
        
            % Compute new variables
            new_var = [v_var; theta_var] - inv(J) * F;
    
            % Add new variables to original variable array
            % v
            for f = 1:size(v_var,1)
                v_all(v_compute_indices(f)) = new_var(f);
                v_var(f) = new_var(f);
            end
            % theta
            for f = 1:size(theta_var,1)
                theta_all(theta_compute_indices(f)) = new_var(size(v_var,1) + f);
                theta_var(f) = new_var(size(v_var,1) + f);
            end
    
            %Increment iterator
            iter = iter + 1;
    
        else
            % Error within tolerance
            iterate = false;
            converges = 1;
        end
    
        if iter == max_iter
            not_converge_disp = ['Did not converge within ', num2str(max_iter), ' iterations'];
            disp(not_converge_disp);
            break
        end
    end
    
    %% Functions to compute ACPF equations
    
    function acpf_p = acpf_p_eq (p_index,n,v_all,theta_all,Y) % returns P equation for (p_index) bus in network
        
        acpf_p = 0;
        p = p_index;
    
        for k = 1:n
            temp = abs(v_all(p)) * abs(v_all(k)) * ...
            (real(Y(p,k)) * cosd(theta_all(p) - theta_all(k)) + imag(Y(p,k)) * sind(theta_all(p) - theta_all(k)));
            acpf_p = acpf_p + temp;
        end
    end
    
    function acpf_q = acpf_q_eq (q_index,n,v_all,theta_all,Y) % returns Q equation for (q_index) bus in network
        
        acpf_q = 0;
        q = q_index;
    
        for k = 1:n
            acpf_q = acpf_q + abs(v_all(q)) * abs(v_all(k)) * ...
            (real(Y(q,k)) * sind(theta_all(q) - theta_all(k)) - imag(Y(q,k)) * cosd(theta_all(q) - theta_all(k)));
        end
    
    end
    
    %% Function to compute Jacobian
    
    function J = compute_jacobian(Y,number_of_buses,v_all,theta_all,p_compute_indices,q_compute_indices,v_compute_indices,theta_compute_indices)
        % [ J11   J12  ; ...
        %   J21   J22  ]
        
        % J11 = del_P / del_v
        % J12 = del_P / del_theta
        % J21 = del_Q / del_v
        % J_22 = del_Q / del_theta
        
        J11 = zeros(size(p_compute_indices,1),size(v_compute_indices,1));
        J12 = zeros(size(p_compute_indices,1),size(theta_compute_indices,1));
        J21 = zeros(size(q_compute_indices,1),size(v_compute_indices,1));
        J22 = zeros(size(q_compute_indices,1),size(theta_compute_indices,1));
        
        % J11 (P , v)
        for i = 1:size(p_compute_indices)
            for j = 1:size(v_compute_indices)
                
                p_i = p_compute_indices(i);
                v_j = v_compute_indices(j);
                
                if p_i == v_j
        
                    temp = 0;

                    for k = 1:number_of_buses
                        v_k = k;

                        if p_i ~= v_k
                            temp = temp + abs(v_all(v_k)) * ( real(Y(p_i,v_k)) * cosd(theta_all(p_i) - theta_all(v_k)) ...
                             + imag(Y(p_i,v_k)) * sind(theta_all(p_i) - theta_all(v_k)) );
                        end
                    end
        
                    J11(i,j) = 2 * abs(v_all(v_j)) * real(Y(v_j,v_j)) + temp;
        
                else
        
                    J11(i,j) = abs(v_all(p_i)) * ( real(Y(p_i,v_j)) * cosd(theta_all(p_i) - theta_all(v_j)) ...
                             + imag(Y(p_i,v_j)) * sind(theta_all(p_i) - theta_all(v_j))  );
        
                end
            end
        end
        
        % J12 (P , theta)
        for i = 1:size(p_compute_indices)
            for j = 1:size(theta_compute_indices)
                
                p_i = p_compute_indices(i);
                t_j = theta_compute_indices(j);
                
                if p_i == t_j
        
                    temp = 0;

                    for k = 1:number_of_buses
                        t_k = k;

                        if p_i ~= t_k
                            temp = temp - abs(v_all(p_i)) * abs(v_all(t_k)) * ( real(Y(p_i,t_k)) * sind(theta_all(p_i) - theta_all(t_k)) ...
                             - imag(Y(p_i,t_k)) * cosd(theta_all(p_i) - theta_all(t_k)) );
                        end
                    end
        
                    J12(i,j) = temp;
        
                else
        
                    J12(i,j) = abs(v_all(p_i)) * abs(v_all(t_j)) * ( real(Y(p_i,t_j)) * sind(theta_all(p_i) - theta_all(t_j)) ...
                             - imag(Y(p_i,t_j)) * cosd(theta_all(p_i) - theta_all(t_j))  );
        
                end
            end
        end
        
        
        % J21 (Q , v)
        for i = 1:size(q_compute_indices)
            for j = 1:size(v_compute_indices)
                
                q_i = q_compute_indices(i);
                v_j = v_compute_indices(j);
                
                if q_i == v_j
        
                    temp = 0;

                    for k = 1:number_of_buses
                        v_k = k;

                        if q_i ~= v_k
                            temp = temp + abs(v_all(v_k)) * ( real(Y(q_i,v_k)) * sind(theta_all(q_i) - theta_all(v_k)) ...
                             - imag(Y(q_i,v_k)) * cosd(theta_all(q_i) - theta_all(v_k)) );
                            
                        end
                    end
                    J21(i,j) = -2 * abs(v_all(v_j)) * imag(Y(v_j,v_j)) + temp;
        
                else
        
                    J21(i,j) = abs(v_all(q_i)) * ( real(Y(q_i,v_j)) * sind(theta_all(q_i) - theta_all(v_j)) ...
                             - imag(Y(q_i,v_j)) * cosd(theta_all(q_i) - theta_all(v_j))  );
        
                end
            end
        end
        
        % J22 (Q , theta)
        for i = 1:size(q_compute_indices)
            for j = 1:size(theta_compute_indices)
                
                q_i = q_compute_indices(i);
                t_j = theta_compute_indices(j);
                
                if q_i == t_j
        
                    temp = 0;

                    for k = 1:number_of_buses
                        t_k = k;

                        if q_i ~= t_k
                            temp = temp + abs(v_all(q_i)) * abs(v_all(t_k)) * ( real(Y(q_i,t_k)) * cosd(theta_all(q_i) - theta_all(t_k)) ...
                             + imag(Y(q_i,t_k)) * sind(theta_all(q_i) - theta_all(t_k)) );
                        end
                    end
                
                    J22(i,j) = temp;
        
                else
        
                    J22(i,j) = -abs(v_all(q_i)) * abs(v_all(t_j)) * ( real(Y(q_i,t_j)) * cosd(theta_all(q_i) - theta_all(t_j)) ...
                             + imag(Y(q_i,t_j)) * sind(theta_all(q_i) - theta_all(t_j))  );
        
                end
            end
        end
        
        J = [J11 J12 ; J21 J22];
    
    end

end