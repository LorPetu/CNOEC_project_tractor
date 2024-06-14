function [z_sim] = Tractor_traj(U,z0,zf,Nu,Ns,parameters,Optimization_opt)
%% Build vector of inputs

Np=ceil((Ns+1)/Nu);

u_in        =   [U(1:Np,1)';
                U(Np+1:2*Np,1)'];

s =    U(2*Np+1:end-1,1);
Ts=     U(end,1);

disp(s)
%% Run simulation with FFD
zdot=zeros(8,1);
z_sim      =   zeros(8,Ns+1);
z_sim(:,1) =   z0;
f_1 = 0;


  for ind=2:Ns+1
    u               =  u_in(:,ceil(ind/Nu));
    zdot               =   Tractor_01_trail_model(z_sim(:,ind-1),u,parameters);
    z_sim(:,ind)       =   z_sim(:,ind-1)+Ts*zdot;

    % Cost function evaluation
     f_1=f_1+0*((z_sim(1:2, ind)-z_sim(1:2, ind-1))'*(z_sim(1:2, ind)-z_sim(1:2, ind-1)));

  end

delta_delta=0;%u_in(1,2:end)-u_in(1,1:end-1);
delta_acc=0;%u_in(2,2:end)-u_in(2,1:end-1);
% 
% f = f +0*(delta_acc*delta_acc')+ 1*(delta_delta*delta_delta')...
%     + 5e3*(s'*s)+1e4*Ts;

fprintf("fcn1: %.2f  delta_acc: %.2f   delta_delta: %.2f  slack_violation: %.2f   T_s: %.2f \r\r",f_1, 0.5*(delta_acc*delta_acc'),1*(delta_delta*delta_delta'),2e2*(s'*s),4e2*Ts)

%% Post-processing the results
Tend=Ns*Ts;

plx     = z_sim(1,:);
ply     = z_sim(2,:);
ang     = z_sim(3,:);
vel     = z_sim(4,:);

vsat = Optimization_opt.vsat;
maxvsat = vsat*ones(Ns+1);
minvsat = -vsat*ones(Ns+1);
deltasat = Optimization_opt.deltasat;
maxdeltasat = deltasat*ones(Np);
mindeltasat = -deltasat*ones(Np);
asat = Optimization_opt.asat;
maxasat = asat*ones(Np);
minasat = -asat*ones(Np);

del     = u_in(1,:);
acc     = u_in(2,:);


plxf = zf(1,1);
plyf = zf(2,1);

time_s=linspace(0,Tend,Ns+1);
time_p=linspace(0,Tend,Np);

% figure(1);
% subplot(5,2,1:6),plot(plx,ply,plxf,plyf,"xr"),daspect([1,1,1]),grid on,axis([-10 10 -10 10])
% xlabel('x'), ylabel('y'),title('traiettoria')
% subplot(5,2,7),plot(time_s,ang),xlabel('Time (s)'),ylabel('psi'),
% subplot(5,2,8),plot(time_s,vel*3.6,time_s,maxvsat*3.6,time_s,minvsat*3.6),xlabel('Time (s)'),ylabel('velocità [km/h]')
% subplot(5,2,9);plot(time_p,del,time_p,maxdeltasat,time_p,mindeltasat),xlabel('Time (s)'),ylabel('delta ')
% subplot(5,2,10);plot(time_p,acc,time_p,maxasat,time_p,minasat),xlabel('Time (s)'),ylabel('acc [m/s]');


end

