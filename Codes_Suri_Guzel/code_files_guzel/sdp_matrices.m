function [Y_k, Y_bar_k, Y_lm, Y_bar_lm, M_k] = sdp_relaxation(mpc)

mpopt = mpoption('out.all', 0);
[Ybus, ~,~] = makeYbus(mpc, mpopt);

N_bus = size(mpc.bus,1);
N_branch = size(mpc.branch,1);
e = eye(N_bus);

Y_k = cell(N_bus,1);
Y_bar_k = cell(N_bus,1);
Y_lm = cell(N_branch,1);
Y_bar_lm = cell(N_branch,1);
M_k = cell(N_bus,1);
% M_lm = cell(N_bus,1);   % I am not sure how this constraint is included
% in mpc files. I didn't include it. 

for i = 1:N_bus
    h_k = e(:,i)* e(:,i)'* Ybus;
    H_Re = 0.5*(h_k + h_k');
    H_Im = (h_k'-h_k)/(0.5*i);
    decompose_hermitian = @(H) [ real(H), -imag(H); imag(H), real(H)];
    Y_k{i} = decompose_hermitian(H_Re);
    Y_bar_k{i} = decompose_hermitian(H_Im);
    M_k{i} = [e(:,i)* e(:,i)', zeros(N_bus);  zeros(N_bus), e(:,i)* e(:,i)'];
end

for i = 1:N_branch
    % BR_B:5, BR_X:4, BR_R:3, T_BUS = 2, F_BUS:1
    y = 1/(mpc.branch(i,3) + 1i*mpc.branch(i,4));
    b = 0.5*1i*mpc.branch(i,5);
    l = mpc.branch(i,1);
    m = mpc.branch(i,2);
    h_lm = (y+b)*e(:,l)*e(:,l)' - y*e(:,l)*e(:,m)';
    H_Re = 0.5*(h_lm + h_lm');
    H_Im = 0.5*1i*(h_lm - h_lm');
    Y_lm{i} = decompose_hermitian(H_Re);
    Y_bar_lm{i} = decompose_hermitian(H_Im);
end
end