%% EM Contamination Function

function gk = EM_contamination(Yl,Fi,pu,tau,K,itnum1,M)

    % initialization-------------------------------------------------------------------------------------------------
    
    for k = 1:K
        gk(k,:,1) = randn(1,M)+1i*randn(1,M); 
    end

    %beta = rand(1,K);
    %beta/sum(beta); 
    beta = 1/(K*2);       

    for i = 1:itnum1
        
        % E-Step-----------------------------------------------------------------------------------------------------
        
        for k = 1:K
            zk(:,:,k) = sqrt(pu)*Fi(:,k)*gk(k,:,i);
        end
        for k = 1:K
            Yk(:,:,k) = zk(:,:,k) + beta*(Yl - sum(zk,3));
        end
        
        
        % M-Step------------------------------------------------------------------------------------------------------
        
        for k = 1:K
            gk(k,:,i+1) = (1/(tau*sqrt(pu)))*Fi(:,k)'*Yk(:,:,k);
        end
        
        
    end
end


 