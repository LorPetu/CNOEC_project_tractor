function [zdot] = Tractor_02_trail_model(z,u,th)
%TRACTOR_MODEL Summary of this function goes here
%% Read parameters, states and inputs

%parameters
Lt      =   th(1,1);                % wheelbase (m)
d       =   th(2,1);                % row width (m)
asat    =   th(3,1);                % maximum accelleration (m/s^2)   
deltamax=   th(4,1);                % maximum steering angle (rad)
deltamin=   th(5,1);                % minimum steering angle (rad)
vsat    =   th(6,1);                % maximum velocity (m/s)
Li      =   th(7,1);                % wheelbase of implement (m)

%states

xt       =   z(1,1);                % inertial X position (m)
yt       =   z(2,1);                % inertial Y position (m)
psit     =   z(3,1);                % yaw angle (rad)
vt       =   z(4,1);                % body x velocity (m/s) 

xi1      =   z(5,1);                % implemen inertial X position (m)
yi1      =   z(6,1);                % implement inertial Y position (m)
psii1    =   z(7,1);                % implement yaw angle (rad)
vi1      =   z(8,1);                % implement body x velocity (m/s) 

xi2      =   z(9,1);                % implemen inertial X position (m)
yi2      =   z(10,1);               % implement inertial Y position (m)
psii2    =   z(11,1);               % implement yaw angle (rad)
vi2      =   z(12,1);               % implement body x velocity (m/s) 


%inputs

deltat   =   u(1,1);                %steering angle
at       =   u(1,2);                %acceleration


%% Model equations

deltai1=psit-psii1;
vi1=vt*cos(deltai1);
deltai2=psii1-psii2;
vi2=vi1*cos(deltai2);


zdot(1,1)=vt*cos(psit);             % xtdot -->  derivative of tractor x position
zdot(2,1)=vt*sin(psit);             % ytdot -->  derivative of tractor y position
zdot(3,1)=vt*tan(deltat)/Lt;        % psitdot -->  derivative of tractor yaw angle
zdot(4,1)=at;                       % vtdot --> derivative of tractor speed


zdot(5,1)=vi1*cos(psii1);           % xi1dot -->  derivative of tractor x position
zdot(6,1)=vi1*sin(psii1);           % yi1dot --> derivative of tractor y position
zdot(7,1)=vi1*tan(deltai1)/Li;      % psii1dot --> derivative of tractor yaw angle
zdot(8,1)=at*cos(deltai1)-vt*sin(deltai1)*(zdot(3,1)-zdot(7,1)); % vi1dot --> derivative of tractor speed

zdot(9,1)=vi2*cos(psii2);           % xi2dot --> derivative of tractor x position
zdot(10,1)=vi2*sin(psii2);          % yi2dot --> derivative of tractor y position
zdot(11,1)=vi2*tan(deltai2)/Li;     % psii2dot --> derivative of tractor yaw angle
zdot(12,1)=zdot(8,1)*cos(deltai2)-vi1*sin(deltai1)*(zdot(7,1)-zdot(11,1)); % vi2dot --> derivative of tractor speed


%   Detailed explanation goes here

end

