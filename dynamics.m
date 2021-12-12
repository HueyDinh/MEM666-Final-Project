% x1 y1 z1 x2 y2 z2 x3 y3 z3 X Y Z is the default order of states
function [d2,lambda] = dynamics(t,x1,y1,z1,x2,y2,z2,x3,y3,z3,X,Y,Z,...
                                     dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3,dX,dY,dZ,...
                                     d2z1,d2z2,d2z3,const)
mc = const.mc; M = const.M; g = const.g; % Unpacking const object
E1x = const.E1x; E2x = const.E2x; E3x = const.E3x;
E1y = const.E1y; E2y = const.E2y; E3y = const.E3y;
a_matrix_trans = [1 0 0 0 0 0 0 0 0 -X+x1-E1x 0         0; % A matrix !transposed!
                  0 1 0 0 0 0 0 0 0 -Y+y1-E1y 0         0;
                  0 0 1 0 0 0 0 0 0 z1-Z      0         0;
                  0 0 0 1 0 0 0 0 0 0         -X+x2-E2x 0;
                  0 0 0 0 1 0 0 0 0 0         -Y+y2-E2y 0;
                  0 0 0 0 0 1 0 0 0 0         z2-Z      0;
                  0 0 0 0 0 0 1 0 0 0         0         -X+x3-E3x;
                  0 0 0 0 0 0 0 1 0 0         0         -Y+y3-E3y;
                  0 0 0 0 0 0 0 0 1 0         0         z3-Z;
                  0 0 0 0 0 0 0 0 0 X-x1+E1x  X-x2+E2x  X-x3+E3x;
                  0 0 0 0 0 0 0 0 0 Y-y1+E1y  Y-y2+E2y  Y-y3+E3y;
                  0 0 0 0 0 0 0 0 0 Z-z1      Z-z2      Z-z3     ];
dbdt = [0 0 -d2z1 0 0 -d2z2 0 0 -d2z3 0 0 0]';  % b stem from terms that are 
M_matrix = [mc 0  0  0  0  0  0  0  0  0  0  0; % not coefficients of generalized velocities
            0  mc 0  0  0  0  0  0  0  0  0  0; % In this case, b contains the velocity control terms.
            0  0  mc 0  0  0  0  0  0  0  0  0;
            0  0  0  mc 0  0  0  0  0  0  0  0; % Mass matrix
            0  0  0  0  mc 0  0  0  0  0  0  0;
            0  0  0  0  0  mc 0  0  0  0  0  0;
            0  0  0  0  0  0  mc 0  0  0  0  0;
            0  0  0  0  0  0  0  mc 0  0  0  0;
            0  0  0  0  0  0  0  0  mc 0  0  0;
            0  0  0  0  0  0  0  0  0  M  0  0;
            0  0  0  0  0  0  0  0  0  0  M  0;
            0  0  0  0  0  0  0  0  0  0  0  M];
F = [0;
     0;
     mc*g; % Other Lower Terms Consist only of Gravitational Forces
     0;
     0;
     mc*g;
     0;
     0;
     mc*g;
     0;
     0;
     M*g];
dadt_trans = [0 0 0 0 0 0 0 0 0 -dX+dx1  0       0;
              0 0 0 0 0 0 0 0 0 -dY+dy1  0       0;
              0 0 0 0 0 0 0 0 0 dz1-dZ   0       0;
              0 0 0 0 0 0 0 0 0 0        -dX+dx2 0;
              0 0 0 0 0 0 0 0 0 0        -dY+dy2 0;
              0 0 0 0 0 0 0 0 0 0        dz2-dZ  0;
              0 0 0 0 0 0 0 0 0 0        0       -dX+dx3;
              0 0 0 0 0 0 0 0 0 0        0       -dY+dy3;
              0 0 0 0 0 0 0 0 0 0        0       dz3-dZ;
              0 0 0 0 0 0 0 0 0 dX-dx1   dX-dx2  dX-dx3;
              0 0 0 0 0 0 0 0 0 dY-dy1   dY-dy2  dY-dy3;
              0 0 0 0 0 0 0 0 0 dZ-dz1   dZ-dz2  dZ-dz3];
dyn_matrix = [M_matrix            -a_matrix_trans; % Dynamic Matrix
              -a_matrix_trans'    zeros(12,12)];
q_dot = [dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3,dX,dY,dZ]';
advance = dyn_matrix^-1*[F;dadt_trans'*q_dot+dbdt]; % Calculate accelerations and lambdas
d2 = advance(1:12); lambda = advance(13:24);
end