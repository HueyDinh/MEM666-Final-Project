function [d2z1,d2z2,d2z3,z1,z2,z3] = vel_control_input(t_series,t_step,const)
     R = const.R;
     freq = const.freq;
     X_series = R*cos(freq*t_series);
     Y_series = R*sin(freq*t_series);
     Z_series = ones(1, length(t_series))*0.2;
     inv_kin_const = @(X,Y,Z) inverse_kinematic(X,Y,Z,const);
     [z1,z2,z3] = arrayfun(inv_kin_const,X_series,Y_series,Z_series);
     pos_control = [z1;z2;z3];
     vel_control = ([pos_control(:,2:end) pos_control(:,end)] - pos_control)/t_step;
     accel_control = ([vel_control(:,2:end) vel_control(:,end)] - vel_control)/t_step;
     d2z1 = accel_control(1,:);
     d2z2 = accel_control(2,:);
     d2z3 = accel_control(3,:);
     z1 = pos_control(1,:);
     z2 = pos_control(2,:);
     z3 = pos_control(3,:);
end

function [z1,z2,z3] = inverse_kinematic(X,Y,Z,const)
 l = const.l;
 E1x = const.E1x; E2x = const.E2x; E3x = const.E3x;
 E1y = const.E1y; E2y = const.E2y; E3y = const.E3y;
 B1x = const.B1x; B2x = const.B2x; B3x = const.B3x;
 B1y = const.B1y; B2y = const.B2y; B3y = const.B3y;
 z1= Z - sqrt(l^2 - (X+E1x-B1x)^2 - (Y+E1y-B1y)^2);
 z2= Z - sqrt(l^2 - (X+E2x-B2x)^2 - (Y+E2y-B2y)^2);
 z3= Z - sqrt(l^2 - (X+E3x-B3x)^2 - (Y+E3y-B3y)^2);
end