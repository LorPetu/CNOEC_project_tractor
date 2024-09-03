function stop = myOutputFunction(U, optimValues, state)
    persistent numViolations Ns constr_param;  
   
    constr_param= evalin('base', 'constr_param');
    Ns = evalin('base', 'Ns');
    stop = false; 

    switch state
        case 'init'
            numViolations = 0;  
        case 'iter'
            if any(optimValues.constrviolation > 1e-6) 
                numViolations = numViolations + 1;
            else
                numViolations = 0; 
            end

            if numViolations > 50  %&& constr_param.c_vel==0
                disp('Ottimizzazione fermata: i vincoli non sono stati soddisfatti per troppe iterazioni consecutive.');
                stop = true;
            end
        case 'done'
           
    end
end