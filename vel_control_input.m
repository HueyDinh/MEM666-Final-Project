function [d2z1,d2z2,d2z3,...% Return the required cartridge acceleration
          dz1,dz2,dz3,... % Return the required cartridge velocity
          z1,z2,z3,... % Return the required cartridge position
          X_series,Y_series,Z_series,... % Desired effector path
          dX_series,dY_series,dZ_series] = vel_control_input(t_series,t_step,const) % Desired effector velocity
     R = const.R;
     freq = const.freq;
     pitch = const.pitch; 
     z_start = const.z_start; % Spiral sketching parameters
     X_series = R*cos(freq*t_series);
     Y_series = R*sin(freq*t_series);
     Z_series = pitch*t_series + z_start;
     dX_series = -R*freq*sin(freq*t_series);
     dY_series = R*freq*cos(freq*t_series);
     dZ_series = pitch*ones(1, length(t_series));
     inv_kin_const = @(X,Y,Z,dX,dY,dZ) inverse_kinematic(X,Y,Z,dX,dY,dZ,const);
     [z1,z2,z3,dz1,dz2,dz3] = arrayfun(inv_kin_const,X_series,Y_series,Z_series,...
                                       dX_series,dY_series,dZ_series);
     vel_control = [dz1;dz2;dz3];
     accel_control = ([vel_control(:,2:end) vel_control(:,end)] - vel_control)/t_step;
     d2z1 = accel_control(1,:);
     d2z2 = accel_control(2,:);
     d2z3 = accel_control(3,:);
end

function [z1,z2,z3,dz1,dz2,dz3] = inverse_kinematic(X,Y,Z,dX,dY,dZ,const) % See report for derivation
 l = const.l;
 E1x = const.E1x; E2x = const.E2x; E3x = const.E3x;
 E1y = const.E1y; E2y = const.E2y; E3y = const.E3y;
 B1x = const.B1x; B2x = const.B2x; B3x = const.B3x;
 B1y = const.B1y; B2y = const.B2y; B3y = const.B3y;
 z1= Z - sqrt(l^2 - (X+E1x-B1x)^2 - (Y+E1y-B1y)^2);
 z2= Z - sqrt(l^2 - (X+E2x-B2x)^2 - (Y+E2y-B2y)^2);
 z3= Z - sqrt(l^2 - (X+E3x-B3x)^2 - (Y+E3y-B3y)^2);
 dz1 = dZ +(l^2 - (X+E1x-B1x)^2 - (Y+E1y-B1y)^2)^(-1/2)*((X+E1x-B1x)*dX+(Y+E1y-B1y)*dY);
 dz2 = dZ +(l^2 - (X+E2x-B2x)^2 - (Y+E2y-B2y)^2)^(-1/2)*((X+E2x-B2x)*dX+(Y+E2y-B2y)*dY);
 dz3 = dZ +(l^2 - (X+E3x-B3x)^2 - (Y+E3y-B3y)^2)^(-1/2)*((X+E3x-B3x)*dX+(Y+E3y-B3y)*dY);
end