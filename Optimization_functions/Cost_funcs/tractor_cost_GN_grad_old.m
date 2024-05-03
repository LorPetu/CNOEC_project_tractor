function [F,zsim] = tractor_cost_GN_grad_old(U,z0,zf,Ts,Q,Qf,Qdot,R,vsat,deltasat,asat,alpha,beta,parameters)


N=length(U)/2;


zsim        =   zeros(4*(N+1),1);
zsim(1:4,1) =   z0;
ztemp       =   z0;



%inizializzazione funzione di costo gauss- newton

F           =   [zeros(4*N,1);
                zeros(4*N,1);
                alpha*exp(-beta*(vsat+ztemp(4,1)));
                alpha*exp(-beta*(vsat+ztemp(4,1)));
                alpha*exp(-beta*(deltasat+U(1,1)));
                alpha*exp(-beta*(deltasat-U(1,1)));
                alpha*exp(-beta*(asat+U(2,1)));
                alpha*exp(-beta*(asat+U(2,1)));
                zeros(2,1);
                zeros(4,1);];


for ind=2:N+1
    % Update the state
    zdot = Tractor_00_trail_model (ztemp,[U((ind-1)*2-1,1);U((ind-1)*2,1)],parameters);
    ztemp                               =   ztemp+Ts*zdot;
    zsim((ind-1)*4+1:ind*4,1)           =   ztemp;
    
   e=zf-ztemp;
    u=[U((ind-1)*2-1,1);
        U((ind-1)*2,1)];
   
   % Update the cost function 
    F((ind-2)*4+1:(ind-1)*4,1)          =   Q*[e(:,1)];
    F((ind-2)*4+4*N+1:(ind-1)*4+4*N,1)          =   Qdot*zdot;
    F((ind-2)*2+8*N+1:(ind-1)*2+8*N,1)              =   R*u;
    F(8*N+2*N+1,1)                    =   alpha*exp(-beta*(vsat+ztemp(4,1)));
    F(8*N+2*N+2,1)                    =   alpha*exp(-beta*(vsat-ztemp(4,1)));
    F(8*N+2*N+3,1)                    =   alpha*exp(-beta*(deltasat+U((ind-1)*2-1,1)));
    F(8*N+2*N+4,1)                    =   alpha*exp(-beta*(deltasat-U((ind-1)*2-1,1)));
    F(8*N+2*N+5,1)                    =   alpha*exp(-beta*(asat+U((ind-1)*2,1)));
    F(8*N+2*N+6,1)                    =   alpha*exp(-beta*(asat-U((ind-1)*2,1)));
    

    F((ind-1)*4,1)                =   i*F((ind-1)*4,1);
   
end

% Update the terminal cost and its Jacobian
F(8*N+2*N+7:8*N+2*N+10,1)                =   Qf*[e(:,1)];




