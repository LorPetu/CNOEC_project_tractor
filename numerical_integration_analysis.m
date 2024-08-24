clear all
close all
clc

warning('off', 'all');
addpath(genpath('Models/'));
addpath(genpath('Optimization_functions/'));
% errore calcolato sulla distanza x e y rispetto a ode45 (ground truth)
%% Model Parameters

Lt        =   3;                  % Wheelbase (m)
Li        =   2;                  % Wheelbase of implements
d         =   4.0;                  % Row width (m)

parameters=[Lt;Li;d];

%% initial states
psit      =    pi/2;              % yaw angle (rad)
psii     =     psit;               % implement yaw angle (rad)
vt        =    20/3.6;             % body x velocity (m/s) 
vi       =   vt;               % implement body x velocity (m/s)
xi       =   0;               % implemen inertial X position (m)
yi       =   0;               % implement inertial Y position (m)
xt        =  xi +Li*cos(psii);                 % inertial X position (m)
yt        =  yi +Li*sin(psii);                 % inertial Y position (m)


z0=[xt;yt;psit;vt;xi;yi;psii;vi];



%% Input variables

j= linspace(0.01, 0.02, 10);
j=[0.01,0.02,0.05,0.1,0.2,0.5];

k=length(j);

tempo=zeros(k,5);

Tend=10;


for k=1:length(j)
       
    Ts=j(k); 
    disp( num2str(Ts))
    Ns=ceil(Tend/Ts);

    delta=[-0.5*ones(1,Ns/2),0.5*ones(1,Ns/2)];
    acc=zeros(1,Ns);
    u             = [delta; acc];
    
    
    errore_ffd=zeros(1,Ns-1);
    errore_rk2=zeros(1,Ns-1);
    errore_rk3=zeros(1,Ns-1);
    errore_rk4=zeros(1,Ns-1);

    zdot        =   zeros(8,1);
    z_prime_rk2 = zeros(8,1);
    z_sim_ffd      =   zeros(8,Ns);
    z_sim_ffd(:,1) =   z0;
    z_sim_rk2     =   zeros(8,Ns);
    z_sim_rk2(:,1) =   z0;
    z_sim_ode     =   zeros(8,Ns);
    z_sim_ode(:,1) =   z0;
 %% simulation using ode45

    tic
    for ind=2:Ns
        
        z_temp= ode45(@(t,z)Tractor_01_trail_model(z, u(:, ind-1), parameters),[0,Ts], z_sim_ode(:, ind-1));
        z_sim_ode(:, ind) = z_temp.y(:,end);

    end 
    tempo(k,1)=toc;

    % figure()
    % plot(z_sim_ode(1,:), z_sim_ode(2,:)), hold on;
    % plot(z_sim_ode(5,:), z_sim_ode(6,:)), daspect([1,1,1]),xlabel('X [m]'),ylabel('Y [m]'),grid on;
    % legend('tractor', 'implement');
    % title('ODE45')
    
    
    %% simulation using rk2
    tic 
    for ind=2:Ns
    
        z_prime_rk2       =   z_sim_rk2(:,ind-1)+Ts/2*Tractor_01_trail_model(z_sim_rk2(:,ind-1),u(:,ind-1),parameters);
        z_sim_rk2(:,ind)  =   z_sim_rk2(:,ind-1)+Ts*Tractor_01_trail_model(z_prime_rk2,u(:,ind-1),parameters);
     
         errore_rk2(k,ind-1)=sqrt((z_sim_rk2(1,ind)-z_sim_ode(1,ind))^2+(z_sim_rk2(2,ind)-z_sim_ode(2,ind))^2);
    end 
    tempo(k, 2)=toc;
    
    
    % figure()
    % plot(z_sim_rk2(1,:), z_sim_rk2(2,:)), hold on;
    % plot(z_sim_rk2(5,:), z_sim_rk2(6,:)), daspect([1,1,1]),xlabel('X [m]'),ylabel('Y [m]'),grid on;
    % legend('tractor', 'implement');
    % title('RK2')
    
    %% RK3
    z_sim_rk3    =   zeros(8,Ns);
    z_sim_rk3(:,1) =   z0;
    tic
    for ind = 2:Ns
        k1 = Tractor_01_trail_model(z_sim_rk3(:,ind-1), u(:,ind-1), parameters);
        z_prime_rk3_1 = z_sim_rk3(:,ind-1) + Ts/2 * k1;
        k2 = Tractor_01_trail_model(z_prime_rk3_1, u(:,ind-1), parameters);
        z_prime_rk3_2 = z_sim_rk3(:,ind-1) + Ts * k2;
        k3 = Tractor_01_trail_model(z_prime_rk3_2, u(:,ind-1), parameters);
        z_sim_rk3(:,ind) = z_sim_rk3(:,ind-1) + Ts/6 * (k1 + 4*k2 + k3);

        errore_rk3(k,ind-1)=sqrt((z_sim_rk3(1,ind)-z_sim_ode(1,ind))^2+(z_sim_rk3(2,ind)-z_sim_ode(2,ind))^2);
    end
    tempo(k,3)=toc;
    
    % figure()
    % plot(z_sim_rk3(1,:), z_sim_rk3(2,:)), hold on;
    % plot(z_sim_rk3(5,:), z_sim_rk3(6,:)), daspect([1,1,1]),xlabel('X [m]'),ylabel('Y [m]'),grid on;
    % legend('tractor', 'implement');
    % title('RK3')
    % 
    %% Rk4
    
    z_sim_rk4    =   zeros(8,Ns);
    z_sim_rk4(:,1) =   z0;
    
    tic 
    for ind = 2:Ns
        k1 = Tractor_01_trail_model(z_sim_rk4(:,ind-1), u(:,ind-1), parameters);
        k2 = Tractor_01_trail_model(z_sim_rk4(:,ind-1) + Ts/2 * k1, u(:,ind-1), parameters);
        k3 = Tractor_01_trail_model(z_sim_rk4(:,ind-1) + Ts/2 * k2, u(:,ind-1), parameters);
        k4 = Tractor_01_trail_model(z_sim_rk4(:,ind-1) + Ts * k3, u(:,ind-1), parameters);
        z_sim_rk4(:,ind) = z_sim_rk4(:,ind-1) + Ts/6 * (k1 + 2*k2 + 2*k3 + k4);
        
        errore_rk4(k,ind-1)=sqrt((z_sim_rk4(1,ind)-z_sim_ode(1,ind))^2+(z_sim_rk4(2,ind)-z_sim_ode(2,ind))^2);
    end
    tempo(k,4)=toc;
    
    % figure()
    % plot(z_sim_rk4(1,:), z_sim_rk4(2,:)), hold on;
    % plot(z_sim_rk4(5,:), z_sim_rk4(6,:)), daspect([1,1,1]),xlabel('X [m]'),ylabel('Y [m]'),grid on;
    % legend('tractor', 'implement');
    % title('RK4')
    
    %% simulation using FFD
    tic
    for ind=2:Ns
        
        zdot               =   Tractor_01_trail_model(z_sim_ffd(:,ind-1),u(:,ind-1),parameters);
        z_sim_ffd(:,ind)       =   z_sim_ffd(:,ind-1)+Ts*zdot;
        
        errore_ffd(k,ind-1)=sqrt((z_sim_ffd(1,ind)-z_sim_ode(1,ind))^2+(z_sim_ffd(2,ind)-z_sim_ode(2,ind))^2);
    end 
    tempo(k,5)=toc;
    
    
    % figure()
    % plot(z_sim_ffd(1,:), z_sim_ffd(2,:)), hold on;
    % plot(z_sim_ffd(5,:), z_sim_ffd(6,:)), daspect([1,1,1]),xlabel('X [m]'),ylabel('Y [m]'),grid on;
    % legend('tractor', 'implement');
    % title('FFD')


    figure()
    plot(linspace(0, Tend, Ns-1), errore_rk2(k,:)); hold on;
    plot(linspace(0, Tend, Ns-1), errore_rk3(k,:))
    plot(linspace(0, Tend, Ns-1), errore_rk4(k,:))
    plot(linspace(0, Tend, Ns-1), errore_ffd(k,:))
    xlabel('tempo'),ylabel('errore'),grid on;
    legend('RK2', 'RK3', 'RK4','FFD');
    title('errors with respect to ODE45', ['Ts=' num2str(j(k))])


    figure()
    plot(z_sim_rk2(1,:), z_sim_rk2(2,:));hold on;
    plot(z_sim_rk3(1,:), z_sim_rk3(2,:));hold on;
    plot(z_sim_rk4(1,:), z_sim_rk4(2,:));hold on;
    plot(z_sim_ffd(1,:), z_sim_ffd(2,:));hold on;
    plot(z_sim_ode(1,:), z_sim_ode(2,:));hold on;
    daspect([1,1,1]),xlabel('X [m]'),ylabel('Y [m]'),grid on;
    legend('RK2', 'RK3', 'RK4','FFD', 'ODE45');
    title('comparison between trajectories', ['Ts=' num2str(j(k))])
end

%%


colnames = {'ODE45','RK2', 'RK3', 'RK4','FFD'};
rownames= string(j);

for i = 1:size(tempo, 5)
    fprintf('\nTable of computation time for each integration method and each sampling time tested:\n', i);
    T = array2table(tempo(:,:,i), 'VariableNames', colnames,'RowNames', rownames); 
    disp(T); 
end




