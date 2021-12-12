clc;clear;
step_size = 0.0001;
t_series = 0:step_size:60;
num_snapshot = length(t_series);
all_q = zeros(12,num_snapshot);
all_dq = zeros(12,num_snapshot);
all_d2q = zeros(12,num_snapshot);
all_lambda = zeros(12,num_snapshot);

const.mc = 0.1; const.M = 0.1; const.l = 0.25369; const.freq = 0.5*pi; const.R = 0.05;
const.E1x = 0.02; const.E2x = -0.01  ; const.E3x = -0.01;
const.E1y = 0   ; const.E2y = 0.01732; const.E3y = -0.01732;
const.B1x = 0.16615; const.B2x = -0.08307; const.B3x = -0.08307;
const.B1y = 0;       const.B2y =  0.14389; const.B3y = -0.14389;
const.g = 9.81;

[all_dzc1, all_dzc2, all_dzc3,...
    all_z1,all_z2,all_z3] = vel_control_input(t_series,step_size,const);

init_q = [const.B1x,const.B1y,all_z1(1),...
          const.B2x,const.B2y,all_z2(1),...
          const.B3x,const.B3y,all_z3(1),const.R,0,0.2]';
init_dq = zeros(12,1);
init_lambda = zeros(12,1);

all_q(:,1) = init_q;
all_dq(:,1) = init_dq;

for i=2:num_snapshot
    all_dq(:,i) = all_d2q(:,i-1)*step_size + all_dq(:,i-1);
    dq = num2cell(all_dq(:,i));

    all_q(:,i) = all_dq(:,i-1)*step_size + all_q(:,i-1);
    q = num2cell(all_q(:,i));

    [d2,lambda] = dynamics(t_series(i),q{:},dq{:},...
                           all_dzc1(i),all_dzc2(i),all_dzc3(i),const);
    all_d2q(:,i) = d2;
    all_lambda(:,i) = lambda;
end

x1_series = all_q(1,:); x2_series = all_q(4,:); x3_series = all_q(7,:);
y1_series = all_q(2,:); y2_series = all_q(5,:); y3_series = all_q(8,:);
z1_series = all_q(3,:); z2_series = all_q(6,:); z3_series = all_q(9,:);
X_series = all_q(10,:); Y_series = all_q(11,:); Z_series = all_q(12,:);
figure
plot3(X_series,Y_series,-Z_series)
figure;
plot(t_series,z1_series)
hold on;
plot(t_series,all_z1,"--")
plot(t_series,z2_series)
plot(t_series,z3_series)
legend
