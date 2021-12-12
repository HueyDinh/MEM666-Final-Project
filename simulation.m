clc;clear; % Clear Console and Work Space
step_size = 0.0001; % Define Time Step Size
t_series = 0:step_size:10; % Generate Time Series of Constant Step Size
num_snapshot = length(t_series); % Number of Entry in Time Series
% States always go in order [x1 y1 z1 x2 y2 z2 x3 y3 z3 X Y Z]
all_q = zeros(12,num_snapshot); % Placeholder Matrix - Store All 0th Order States of All Time
all_dq = zeros(12,num_snapshot); % Placeholder Matrix - Store 1st Order States
all_d2q = zeros(12,num_snapshot); % Placeholder Matrix - Store 2nd Order States
all_lambda = zeros(12,num_snapshot); % Placeholder Matrix - Store Lagrange Multiples

L_big = 0.3; % Side Length of the Base (Big) Triangle
L_small= 0.035; % Side Length of the Effector (Small) Triangle

% Constant Definition (all in fundamental SI units):
% mc: mass of cartridge; M: mass of effector body, l: linkage length;
const.mc = 0.1; const.M = 0.4; const.l = 0.25; const.g = 9.81;
% freq,R,pitch: Parameters for Spiral Sketching (see: vel_control_input.m)
const.freq = 2*pi; const.R = 0.05; const.pitch = -0.005; const.z_start = 0.30;
% B = Frame Pillars Coordinate in Global Frame, E = Linkage Connection
% Points in Effector Frame. 1,2,3 = Cartridge number
const.E1x = L_small/sqrt(3); const.E2x = -L_small/sqrt(3)/2; const.E3x = -L_small/sqrt(3)/2;
const.E1y = 0              ; const.E2y = L_small/2         ; const.E3y = -L_small/2;
const.B1x = L_big/sqrt(3); const.B2x = -L_big/sqrt(3)/2; const.B3x = -L_big/sqrt(3)/2;
const.B1y = 0            ; const.B2y = L_big/2         ; const.B3y = -L_big/2;

% Performing Inverse Kinemtic Analysis (Given Desired X,Y,Z movment) to:
% a) Determine neccesary z1,z2,z3 movement
% b) Determine the Coorect Initial Condition.
[all_d2z1, all_d2z2, all_d2z3,... % Neccessary cartridge acceleration
 all_dz1, all_dz2, all_dz3,... % Neccessary cartridge velocity 
 all_z1,all_z2,all_z3,... % Neccessary cartridge position
 all_X,all_Y,all_Z,... % (Designed) Effector Coordinates
 all_dX,all_dY,all_dZ] = vel_control_input(t_series,step_size,const); % Effector Velocity

init_q = [const.B1x,const.B1y,all_z1(1),...
          const.B2x,const.B2y,all_z2(1),...
          const.B3x,const.B3y,all_z3(1),...
          all_X(1),all_Y(1),all_Z(1)]';
init_dq = [0 0 all_dz1(1) 0 0 all_dz2(1) 0 0 all_dz3(1)...
           all_dX(1) all_dY(1) all_dZ(1)]';
init_lambda = zeros(12,1);

all_q(:,1) = init_q;
all_dq(:,1) = init_dq;

for i=2:num_snapshot % Simulation Loop
    all_dq(:,i) = all_d2q(:,i-1)*step_size + all_dq(:,i-1); % Euler Method to Calculate Current 1st Order States
    dq = num2cell(all_dq(:,i)); % Package for Input

    all_q(:,i) = (all_dq(:,i-1)+all_dq(:,i))*step_size/2 + all_q(:,i-1); % Trapezoidal Method to Calculate Current 0th Order States
    q = num2cell(all_q(:,i)); % Package for Input

    [d2,lambda] = dynamics(t_series(i),q{:},dq{:},...
                           all_d2z1(i),all_d2z2(i),all_d2z3(i),const); % Iterate through Dynamics Matrix (see: dynamics.m)
    all_d2q(:,i) = d2; % Logging Values
    all_lambda(:,i) = lambda;
end

x1_series = all_q(1,:); x2_series = all_q(4,:); x3_series = all_q(7,:);
y1_series = all_q(2,:); y2_series = all_q(5,:); y3_series = all_q(8,:);
z1_series = all_q(3,:); z2_series = all_q(6,:); z3_series = all_q(9,:);
X_series = all_q(10,:); Y_series = all_q(11,:); Z_series = all_q(12,:);

f1 = all_lambda(3,:) ; f2 = all_lambda(6,:) ; f3 = all_lambda(9,:);
s1xy = all_lambda(1:2,:); s2xy = all_lambda(4:5,:); s3xy = all_lambda(7:8,:);

rot_120deg = [cosd(120)  sind(120);
              -sind(120) cosd(120)];

s1xy_T = s1xy;
s2xy_T = rot_120deg*s2xy;
s3xy_T = rot_120deg*rot_120deg*s3xy;

% 3D Plots of the Trajectory (Both Effector and Cartridges)
figure;
plot3(X_series,Y_series,Z_series,"DisplayName","Effector") % Note: Positive Z Points Down. The Picture is Up-Side-Down
title("Cartridge and Effector Global Coordinate Trace Over Time");
hold on;
plot3(x1_series,y1_series,z1_series,"DisplayName","Cartridge 1");
plot3(x2_series,y2_series,z2_series,"DisplayName","Cartridge 2");
plot3(x3_series,y3_series,z3_series,"DisplayName","Cartridge 3");
legend("Location","best");
xlabel("x (m)"); ylabel("y (m)"); zlabel("z (m)");
hold off;
grid on; grid minor;

% Comparison Between Commanded z Values and Actual z Values of Cartridge
figure;
plot(t_series,z1_series,"DisplayName","z1 Actual")
title("Comparison Between Commanded z Values and Actual z Values of Cartridge")
hold on;
plot(t_series,all_z1,"-.","DisplayName","z1 Command")
plot(t_series,z2_series,"DisplayName","z2 Actual")
plot(t_series,all_z2,"-.","DisplayName","z2 Command")
plot(t_series,z3_series,"DisplayName","z3 Actual")
plot(t_series,all_z3,"-.","DisplayName","z3 Command")
legend("Location","best");
xlabel("Time (s)");
ylabel("Z Displacement (m)")

% Control Force
figure;
plot(t_series,f1,"DisplayName","F1")
title("Cartridge Control Forces Over Time")
hold on;
plot(t_series,f2,"DisplayName","F2")
plot(t_series,f3,"DisplayName","F3")
xlabel("Time (s)")
ylabel("Force (N)")
legend;
% Frame Shaking
figure;
plot(s1xy_T(1,:),s1xy_T(2,:),"DisplayName","Rail 1")
title("Cartridge 1 Frame Shake")
xlabel("Radial (N)"); ylabel("Tangential (N)")

figure;
plot(s2xy_T(1,:),s2xy_T(2,:),"DisplayName","Rail 2")
title("Cartridge 2 Frame Shake")
xlabel("Radial (N)"); ylabel("Tangential (N)")

figure;
plot(s3xy_T(1,:),s3xy_T(2,:),"DisplayName","Rail 3")
title("Cartridge 3 Frame Shake")
xlabel("Radial (N)"); ylabel("Tangential (N)")
