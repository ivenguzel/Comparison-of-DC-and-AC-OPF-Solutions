function [V_opt_abs, V_opt_angle] = decompose_W(mpc, Y_k, Y_bar_k, Y_lm, Y_bar_lm, M_k, W_opt)
    if  rank(W_opt, 1e-5) >= 3
        fprintf('Rank Constraint is not satisfied')
        V_opt_abs = 0;
        V_opt_angle = 0;
    else
        % Eigenvalue decomposition of W_opt
        tol = 1e-5;
        [V,D] = eig(W_opt);
        
        % vectors corresponding to the nonzero eigenvalues of W
        [nz,~] = find(abs(D)>=tol);
        v1 = V(:,nz(1));
        v2 = V(:,nz(2));
        
        % First equation: use the slack bus info to eliminate one of the zeta values
        n_slack = 1; % Slack bus is the first bus in our example
        n = length(v1)/2;
        c = -v1(n+n_slack)/v1(n_slack);
        v_app = v1 + c*v2;
        
        % Second Equation: Sync. Condensor at Bus 3 does not generate power.
        % Then Tr(Y_k{3}.W) ==-Pd =-0.95
        zeta_1 = sqrt(-0.95/trace(Y_k{3}*(v_app*v_app')));
        zeta_2 = zeta_1*c;
        V_opt = (zeta_1 +1i*zeta_2)*( v1(1:n) + 1i*v1(n+1:end));
        
        % Check feasibility 
        X_opt = [real(V_opt);imag(V_opt)];
        W_new = X_opt*X_opt';
        
        % P and Q generation
        Sbase = mpc.baseMVA;
        Pd = mpc.bus(:,3) / Sbase;
        Qd = mpc.bus(:,4) / Sbase; 
        Pinj = [trace(Y_k{1}*W_new); trace(Y_k{2}*W_new);trace(Y_k{3}*W_new)];
        Pg = Pinj + Pd
        
        Qinj = [trace(Y_bar_k{1}*W_new); trace(Y_bar_k{2}*W_new);trace(Y_bar_k{3}*W_new)];
        Qg = Qinj + Qd
        
        % Bus voltage magnitude
        V_bus_mag = sqrt([trace(M_k{1}*W_new); trace(M_k{2}*W_new); trace(M_k{3}*W_new)])

        % Line flow constraint between Line-2-3
        S_23 = sqrt(trace(Y_lm{2}*W_new)^2 + trace(Y_bar_lm{2}*W_new)^2)
        
        % Results
        [V_opt_angle, V_opt_abs]  = cart2pol(real(V_opt),imag(V_opt));
        V_opt_angle = rad2deg(V_opt_angle);
        V_opt_angle = V_opt_angle + (0-V_opt_angle(n_slack));
        for i=1:n
            if V_opt_angle(i) >=180 
                V_opt_angle(i) = V_opt_angle(i) -360;
            end
        end
        Cost_verified = 0.11*(Pg(1)^2)*Sbase^2 + 5*Pg(1)*Sbase + 0.085*(Pg(2)^2)*Sbase^2 + 1.2*Pg(2)*Sbase
    end
end