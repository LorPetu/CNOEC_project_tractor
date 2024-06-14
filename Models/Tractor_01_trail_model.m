function [zdot] = Tractor_01_trail_model(z,u,parameters)
%TRACTOR_MODEL Summary of this function goes here
%% Read parameters, states and inputs

%parameters
Lt      =   parameters(1,1);                % wheelbase (m)
Li      =   parameters(2,1);                % wheelbase of implement (m)
Lh      =   parameters(3,1);                   % hook length
d       =   parameters(4,1);                % row width (m)

%states

xt       =   z(1,1);                % inertial X position (m)
yt       =   z(2,1);                % inertial Y position (m)
psit     =   z(3,1);                % yaw angle (rad)
vt       =   z(4,1);                % body x velocity (m/s) 

xi       =   z(5,1);               % implemen inertial X position (m)
yi       =   z(6,1);               % implement inertial Y position (m)
psii     =   z(7,1);               % implement yaw angle (rad)
vi       =   z(8,1);               % implement body x velocity (m/s) 


%inputs

deltat   =   u(1,1);                %steering angle
at       =   u(2,1);                %acceleration

%% Model equations

deltai=psit-psii;
vi=vt*cos(deltai);


zdot(1,1)=vt*cos(psit);             % xtdot -->  derivative of tractor x position
zdot(2,1)=vt*sin(psit);             % ytdot -->  derivative of tractor y position
zdot(3,1)=vt*tan(deltat)/Lt;        % psitdot -->  derivative of tractor yaw angle
zdot(4,1)=at;                       % vtdot --> derivative of tractor speed

zdot(5,1)=vi*cos(psii);           % derivative of tractor x position
zdot(6,1)=vi*sin(psii);           % derivative of tractor y position
zdot(7,1)=vi*tan(deltai)/Li;      % derivative of tractor yaw angle
zdot(8,1)=at*cos(deltai)-((vt^2*sin(deltai)*tan(deltat))/Lt)+((vt^2*sin(deltai)^2)/Li); % derivative of tractor speed


end

