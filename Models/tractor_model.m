function [zdot] = Tractor_model(z,u,parameters)
%TRACTOR_MODEL Summary of this function goes here
%% Read parameters, states and inputs

% u=[u1;
%    u2];


%parameters
Lt      =   parameters(1,1);                 % Wheelbase (m)
Hi      =   parameters(2,1);                 % Initial heading of the tractor (rad)
Hf      =   parameters(3,1);                 % Final heading of the tractor (rad)
d       =   parameters(4,1);                 % Row width (m)
Li      =   parameters(5,1);

%states

xt       =   z(1,1);            % inertial X position (m)
yt       =   z(2,1);            % inertial Y position (m)
psit     =   z(3,1);            % yaw angle (rad)
vt       =   z(4,1);            % body x velocity (m/s) 

%inputs

deltat   =   u(1,1);              %steering angle
at       =   u(2,1);              %acceleration


%% Model equations

zdot(1,1)=vt*cos(psit);
zdot(2,1)=vt*sin(psit);
zdot(3,1)=vt*tan(deltat)/Lt;
zdot(4,1)=at;


end

