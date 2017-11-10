%% DYNAMIC MODEL OF DDM
const.DDMMode         = 0;                  % MOTOR OR PUMP MODE? 1 is pump else motor
const.Optim           = 0;
const.Temp            = 40+273.15;          % Oil temperature
const.rho             = 872;                % Fluid density [kg/m^3] 870
const.rho_air         = 1.1455;             % Air density [kg/m^3]
const.alpha_mu        = 2.1e-8;             % Fluid coefficient
const.mu              = 0.0000633361*exp(879.7742/((const.Temp)-177.7865))*exp(const.alpha_mu*5e5);   % Dynamic fluid viscosity [kg/(m s)]
const.nu              = const.mu/const.rho; % Fluid viscosity 46e-6 217e-6
const.mu_0            = 0.0221;             % Dynamic fluid viscosity [kg/(m s)]
const.l_stroke        = 2.5e-3;             % Valve stroke length
const.ddotx_0         = 0;                  % Initial acceleration
const.alpha           = 0.005;              % Fluid coefficient from 0.005
const.E_0             = 11500e5;            % Fluid coefficient stiffness from 15500
const.m_coef          = 11.4;               % Fluid coefficient
const.kappa           = 1.4;                % Fluid coefficient
const.p_0             = 101325;             % Atmospheric pressure [Pa]
const.epsilon         = 1e-16;
const.R_air           = 286.9;
const.k_air           = 1.4;
const.c_p             = const.R_air*const.k_air/(const.k_air-1);
% Motor/Pump specific
const.p_H             = 360e5;              % High pressure mfd
const.p_L             = 10e5;                % Low pressure mfd
const.RPM             = 12.1;               % Rated rotations per minute (pump)
const.RPM_M           = 800;                % Rated rotations per minute (motor)
const.dottheta_ini    = const.RPM/60*2*pi;  % Angular velocity (pump)
const.dottheta_ini_M  = const.RPM_M/60*2*pi;% Angular velocity (motor)
const.theta_ini       = 10*pi/180;%pi-thetaOffset;   CHECK THIS!!
% Lobe Pump parameters
const.A_c_p_P         = pi*(45e-3)^2;       % Cross-section area of piston
const.r_P             = 75e-3/2;            % Cylinder stroke [m]
const.lobes           = 16;                 % Amount of pump lobes [-]
const.p_HP            = 360e5;              % High pressure mfd (pump) [Pa]
const.p_LP            = 10e5;               % Low pressure mfd (pump) [Pa]
% Eccentric Motor parameters
const.A_c_p           = pi*(18e-3)^2;       % Piston Area
const.l_strokecyl     = 49.12e-3;           % previous data(4/pi*50e-6)^(1/3);
const.Cylinder_volume = const.l_strokecyl*const.A_c_p;
const.V_0             = 62e-6;              % Cylinder volume at initial position
const.r_e             = const.l_strokecyl/2;% Exentric radius
% Valve initial positions
const.Transient       = 1;                  % Quasi-static orifice equation or with a transient flow development?
const.R_max           = 0.1e-5;
const.h_min           = 0*1e-6;
const.h_minstic       = 0.1e-7;
const.h_max           = 2.5e-3;
const.hdotini         = 0*0.00001;
const.idle            = 0;
const.p_ini1_M        = const.p_L+0.5e5;    % Pressure resulting from position
const.z_Hini1_M       = const.h_min;   % Open or closed?
const.z_Lini1_M       = const.h_max;   % Open or closed?
const.Q_Hini          = 0;
const.Q_Lini          = const.A_c_p*-sin(const.theta_ini)*const.r_e*const.dottheta_ini_M;
% Switching parameters and leakage
const.Q_leakH         = 0/(60000*(const.p_H-const.p_L)); % Leakage coeficient
const.Q_leakL         = 0/(60000*(const.p_H-const.p_L)); % Leakage coeficient
const.theta_HP_M      = (305.5+8)*pi/180;         % Valve switching timing 330.5
const.theta_LP_M      = (142.5+5.5*0.5)*pi/180;         % Valve switching timing 162.5
const.theta_HP1_M     = mod(const.theta_HP_M,2*pi); % Valve switching timing (only one cylinder considered)
const.theta_LP1_M     = mod(const.theta_LP_M,2*pi); % Valve switching timing (only one cylinder considered)
const.theta_HP_P      = pi;                 % Valve switching timing
const.theta_LP_P      = 0;                  % Valve switching timing
% Plunger dimensions
const.m               = 0.017+0*0.002;        % Mass of moving armature ~ from Solidworks [kg] started at 0.07 [kg]
const.L_ring          = 0.7e-3;             % Length of one contact point
const.r_out1          = 18.3e-3;            % Radius of outer ring 1
const.r_in1           = const.r_out1-const.L_ring;%r2+r3-0.50311e-3; % Radius of inner ring 1
const.r_out2          = 11.5e-3+const.L_ring;% Radius of outer ring 2
const.r_in2           = 11.5e-3;            % Radius of inner ring 2
const.L_1             = (const.L_ring)/2;   % Length of outer contact surface
const.L_2             = (const.L_ring)/2;   % Length of inner contact surface
const.L_top           = 0.61e-3/2;          % Length for stiction in full stroke
const.L_inlet         = 15e-3;
const.R1              = 11.5e-3-const.L_1;
const.R2              = const.r_out1-const.L_2;
const.d_mean          = (const.r_in1^2+const.r_out2^2)/2;
const.A_pres          = pi*(const.d_mean);      % [m^2] Effective area for the cylinder pressure on the valves
const.A_s             = 4*pi*(const.L_1*const.R1+const.L_2*const.R2); % From DBR but factor two is included due to difference in definitiion of L
const.A_flow          = pi*(const.r_out1^2-const.r_in2^2);
const.A_flow2         = pi*(const.r_in1^2-const.r_out2^2); % F = pi*((r2-r1)*deltap+p2)-(r4-r3)*p2)
const.D_H             = 2*(const.r_in1-const.r_out2);
const.c_turb          = 0.61;               % Discharge coefficient
const.R_t             = 20;                 % Critical reynolds number
const.p_cr            = const.rho/2*(const.nu*const.R_t/(const.c_turb*const.D_H))^2;
const.fitA_P1         = 2.3870e-04;
const.fitA_P2         = -0.1789;
const.fitA_P3         = -2.7390e-05;
const.r_1             = 11.5e-3;
const.l_1             = 0.35e-3;
const.l_2             = 0.6e-3;
const.W               = 5e-3;
A_H                   = @(l_2) ((const.r_1+4*const.l_1+2.*l_2+const.W).^2-(const.r_1).^2)*pi;
A_L                   = @(l_2) ((const.r_1+2*const.l_1+l_2+const.W).^2-(const.r_1+2*const.l_1+1*l_2).^2)*pi;
A_C                   = @(l_2) ((const.r_1+const.l_1+l_2).^2-(const.r_1+const.l_1).^2+(const.r_1+3*const.l_1+2*l_2+const.W).^2-(const.r_1+3*const.l_1+l_2+const.W).^2)*pi;

% Spring constants
const.F_ini           = 31;                 % Spring preload to compensate for flow force
const.k_s1            = 2.22e3;             % Spring constant 1 [N/m]
const.k_s2            = 31.04e3;            % spring constant 2 [N/m]
const.B_end           = 1e3*1e1;            % End damping upon reaching end stop [Ns/m]
const.K_end           = 1.5e8*1e1;          % End spring [N/m] 1.4e17 from ANSYS
const.tol             = 1e-3;               % Tolerence before damping starts
% Inertia and damp of shaft
const.J_shaft         = 25;                 % Estimated iniertia on the shaft
const.B_shaft         = 0.3*const.J_shaft*900/(const.RPM_M/60*2*pi); % Corresponding damping to ensure constant speed
% Fluid movement coefficients
const.h_b             = 10.4e-3;            %
const.g_1             = 0.05e-3;            % Shearing gap length
const.g_2             = 0.05e-3;            % Shearing gap length
const.g_iarma         = 5e-5;               % Shearing gap length of plunger
const.r_i1            = 14.54133e-3;        % Radius of inner coil
const.r_i2            = const.r_i1+const.g_1; % Radius of inner coil 2
const.r_o1            = 15.6687e-3;         % Radius of outer coil
const.r_o2            = const.r_o1+const.g_2;% Radius of outer coil
const.r_iarma         = 4e-3;               % Radius of inner plunger part
const.h_iarma         = 8.5e-3;             % Height of inner plunger part
const.A_bin           = 2*pi*const.r_i2*const.h_b; % Inner surface area
const.A_bout          = 2*pi*const.r_o1*const.h_b; % Outer surface area
const.A_iarma         = 2*pi*const.r_iarma*const.h_iarma/const.g_iarma;
const.L_pipe          = 12.5e-3;            % Length of the venting bores
const.D_pipe          = 2e-3;               % Diameter of ventin bores
const.B_bcoef         = const.mu*(const.A_bin/const.g_1+const.A_bout/const.g_2);
const.A_top           = (const.r_o1^2-const.r_i2^2)*pi;
const.V_top           = 2.5e-3*pi*(const.r_o2^2-const.r_i1^2)+(19.2365e-3^2-14.63411e-3^2)*pi*1.1e-3;
const.K_hagen1        = pi*const.r_i2*const.g_1^3/(6*const.mu*const.h_b);
const.K_hagen2        = pi*const.r_o2*const.g_2^3/(6*const.mu*const.h_b);
const.K_hagen3        = 128*const.mu*const.L_pipe/(pi*const.D_pipe^4);
const.d_mean          = (17.6e-3^2+12.4e-3^2)/2; % Mean diameter of
const.A_c_valves      = pi*(const.d_mean);    % [m^2] Effective area for the cylinder pressure on the valves
const.K_a             = 8/3*const.r_out1^3*const.rho; % FROM DBR paper DFP draft
const.K_v             = 16*const.mu*const.r_out1;
const.K_d             = 1.1*const.rho*0.5*pi*const.r_out1^2; % 1/2*c_d*rho*A   prøvet 1.17
const.K_h             = 128/(3*pi)*15e-3^2*const.rho*sqrt((const.mu/const.rho)/pi);
const.Flow_gain       = 1e3; %0.5e3
load('K_a1.mat');                           % LPM of added mass
load('K_v1.mat');                           % LPM of viscous damping
load('K_d1.mat');                           % LPM of damping coef
%Stribeck friction on plunger
const.F_cou           = 0.4*7.5;            % Coulomb friction 7.0034, 4.7456
const.z_br            = 0.05;               % Break away velocity 0.013, 0.4799
const.z_th            = 1e-4;               % "zero" velocity to avoid a singularity
const.mu_dry          = 0.4*7.5;            % Dry friction, tried: 4.0366, 3.8350
const.F_pin           = 8.5;                % From the spring with wire connection 11.1374  4.9959  8.5;
const.K_fac           = 1;                  % Scaling factor for testing purpose
const.Q_br            = 10/60000;
const.F_str           = 0;
%% Flow area functions
directory           = strcat('C:\MAIN\CFD\CFD_NCB25_big_motor\Variations\areatest\File_output_area'); % define your directory   Pressure_tests_002\area
file                = sprintf('%s%s',directory,strcat('\Analysis_data.out')); % _updown0
Analysis_data3      = importdata(file);
x_bot               = find(Analysis_data3.data(:,7) >= 2.485);
x_topend            = find(Analysis_data3.data(:,7) < 0);
A_c_valves_pfit     = fit(Analysis_data3.data(8:x_bot(1),7)*1e-3,-Analysis_data3.data(8:x_bot(1),15)./Analysis_data3.data(8:x_bot(1),8),'poly4');
A_c_valves_pospfit  = fit(Analysis_data3.data(x_bot(end):x_topend(1),7)*1e-3,-Analysis_data3.data(x_bot(end):x_topend(1),15)./Analysis_data3.data(x_bot(end):x_topend(1),8),'poly3');

directory           = strcat('C:\MAIN\CFD\CFD_NCB25_big_motor\Variations\Test_1\File_output'); % define your directory
file                = sprintf('%s%s',directory,strcat('\Analysis_data_updown04.out'));
Analysis_data3      = importdata(file);

directory           = strcat('C:\MAIN\');   % define your directory - is needed again because of loaded data from another folder

%% Actuator coefficients
F_act                 = -100;               % Actuator force
const.Voltage         = 70;                 % Possible applied voltage from -80
const.L               = 0.77e-3;            % Estimated inductance of the coil
const.N               = 66;                 % Number of coil windings
const.r_coil          = 15e-3;              % Radius of the coil
const.B_gap           = 0.8;                % Mean flux density in the air gap
const.Actuator_gain   = -2*pi*const.r_coil*const.B_gap*const.N;  %-2*pi*15e-3*66*0.8;% Linear actuator gain [N/A] L*B B~[0.7-0.8]
const.R_coil          = 42.18/15.88;        % Resitance in the coil
const.back_gain       = 1;                  % Factor of back emf gain
const.eta_act         = 1;

%% Stiction specific
const.A_stop          = 4*pi*const.L_top*(6.3e-3+5.7e-3)/2;
const.delta1          = (const.r_out1 - const.r_in1)/const.r_in1;   % Stiction factor from Scheidl et al.
const.delta2          = (const.r_out2 - const.r_in2)/const.r_in2;   % Stiction factor from Scheidl et al.
const.delta3          = (6.3e-3 - 5.7e-3)/5.7e-3;                   % Stiction factor from Scheidl et al.
const.K_stic1         = const.mu*const.r_in1^4*const.delta1^3;      % Stiction factor for lower outer part from Scheidl et al.
const.K_stic2         = const.mu*const.r_in2^4*const.delta2^3;      % Stiction factor for lower inner part from Scheidl et al.
const.K_stic3         = const.mu*5.7e-3^4*const.delta3^3;           % Stiction factor for top part of plunger from Scheidl et al.

%% Hertizian pressure stuff
const.E1 = 96e9;
const.E2 = 200e9;
const.nu_1 = 0.36;
const.nu_2 = 0.3;

% Common simulation settings
Decimation      = 2;
Sample_time     = -1;
factor          = 5e-3;             % reduction of the INSANE sampletime
sample          = 1/8e-7;
fixed_time_step = 1/sample*1e-1;         % Time step 0.5e-6 by DBR

%% C_D
for i = 1:60
    directory2 = sprintf('C:/MAIN/CFD/CFD_static/Variations/Areaprestest_%s/File_output',strcat(num2str(i))); % define your directory
    file = sprintf('%s%s',directory2,strcat('\Analysis_data.out'));
    Analysis_data = importdata(file);
    field1 = 'F_fluid'; value1 = {Analysis_data.data(:,1)};
    field2 = 'Flow1'; value2 = {Analysis_data.data(:,2)};
    field3 = 'Flow2'; value3 = {Analysis_data.data(:,3)};
    staticCFDstruct3(i) = struct(field1,value1,field2,value2,field3,value3);
end
Openings = [2.5 2.0 1.5 1.0 0.5 0.15]*1e-3;
pres_seq = (linspace(4.32e5,5.68e5,10)-5e5);
Q = @(deltap,z) (2*pi*z*(const.r_in1+const.r_out2)+const.epsilon)*sqrt(2/const.rho*abs(deltap))*sign(deltap);
counter = 0;
for k = 1:length(pres_seq)
    for i = 1:length(Openings)
        counter = counter+1;
        C_d(k,i) = staticCFDstruct3(counter).Flow1(end)/const.rho/Q(pres_seq(k),Openings(i));
    end
end