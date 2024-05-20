%% ##########  Start gridding from here ######################
% ##################################

m_iter   =[0.5 0.4 0.3 0.2 0]';
q_iter   =[10 8 7 6 5.95 5.5 5]'; 

for m_index=1:size(m_iter)

    for q_index=1:size(q_iter)
        clearvars -except m_iter m_index q_iter q_index ;
        try
            [out] = gridding_opt_routines(m_iter(m_index),q_iter(q_index));
        catch E
            warning(E.identifier,'Error occurred: %s', E.message);
            % nothing
            continue;
        end

    
    end
    
end


