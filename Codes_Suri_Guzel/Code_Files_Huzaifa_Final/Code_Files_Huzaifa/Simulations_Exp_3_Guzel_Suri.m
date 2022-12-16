clear
close all

%% Experiment 3 (Comparison of LMPs computed using (1) AC OPF and (2) DC OPF)

%----------------------------------------------------------------------------

% Output:
    % lmp_ac_opf: LMPs for all buses computed using AC OPF
    % lmp_dc_opf: LMPs for all buses computed using DC OPF

%----------------------------------------------------------------------------

% Define casefile
casefile = case9;

% Load case file
mpc = loadcase(casefile);

%% Matpower

% Get LMPs from DC OPF
dc_opf_matpower = dcopf(mpc);
lmp_dc_opf_matpower = dc_opf_matpower.bus(:,14);

% Get LMPs from AC OPF
ac_opf_matpower = opf(mpc);
lmp_ac_opf = ac_opf_matpower.bus(:,14);

%% Self

% Get LMPs from DC OPF
[~, ~, lmp_self] = DC_OPF_Function_Guzel_Suri(mpc);
lmp_dc_opf_self = lmp_self;

% Get LMPs from AC OPF
ac_opf_matpower = opf(case5_test);
lmp_ac_opf = ac_opf_matpower.bus(:,14);
