function d = estdiagpost(output,Q,sigma2,HybridOptions)
% This function is to approximate the diagonal entries of the posterior 
% covariance matrix
% 
%   output        - output file from hybrid methods that contains 
%                   matrices B and V generated from Golub-Kahan
%                   bidiagonalization
%   Q             - prior covariance matrix
%   sigma2        - prior noise variance
%   HybridOptions - {'genhybr','meanest'}
%
% Latest Update by T.Cho, Nov. 28, 2021   
    
    switch HybridOptions
        case 'genhybr'
            n_beta = 0;
            n_s = size(Q,1);
            
            Bk = output.B; k = size(Bk,2);
            QV = zeros(n_s + n_beta, k);
            for i = 1:k
                QV(:,i) = Q*output.V(:,i);
            end
            
            I = speye(size(QV,1));
            Q_ds = full(I(:,1))'*(Q*full(I(:,1))); 
            
            lambda = output.alpha;
            
            [~,s,v] = svd(Bk,0);
            s = diag(s);
            d12 = s./sqrt(s.^2 + lambda.^2);
            invL = (QV*v)*diag(d12);

            d = Q_ds*ones(n_s,1);
            d = d - diaglowrank(invL);
            d = d./lambda.^2*sigma2;
            
        case 'meanest'
            [n_s,n_beta] = size(Q.W);
            
            Bk = output.B; k = size(Bk,2);
            QV = zeros(n_s + n_beta, k);
            for i = 1:k
                QV(:,i) = Q*output.V(:,i);
            end
            
            I = speye(size(QV,1));
            Q_ds = full(I(:,1))'*(Q*full(I(:,1))); 
            Q_dbeta = full(I(:,end))'*(Q*full(I(:,end))); 
            
            lambda = output.alpha;
            
            [~,s,v] = svd(Bk,0);
            s = diag(s);
            d12 = s./sqrt(s.^2 + lambda.^2);
            invL = (QV*v)*diag(d12);

            if length(sigma2) == 1
                d = sigma2*[Q_ds*ones(n_s,1); Q_dbeta* ones(n_beta,1)];
            else
                d = [sigma2.*(Q_ds*ones(n_s,1)); Q_dbeta* ones(n_beta,1)];
            end

            d = d - diaglowrank(invL);
            d = d./lambda.^2;
    end
       
end

    
function d = diaglowrank(V)

    n = size(V,1);
    d = zeros(n,1);
    for i = 1:n
        d(i) = V(i,:)*V(i,:)';
    end
    
end