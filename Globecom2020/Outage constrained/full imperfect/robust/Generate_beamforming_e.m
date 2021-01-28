function [e_opt,flag] = Generate_beamforming_e(N, M, K, H, G, H_error, G_error,...
                F_ini, e_ini, x_ini, prob, noise_maxpower, trans_maxpower, rate_min)
            
%%%%%  Paramenter generation  %%%%%
for k=1:K
    gamma(k)=2^(rate_min)/(2^(rate_min)-1);
   

    PHI(:,:,k)=2^(rate_min)/(2^(rate_min)-1)*F_ini(:,k)*F_ini(:,k)'-F_ini*F_ini';
     [t1,t2]=eig(G(:,:,k)*PHI(:,:,k)*G(:,:,k)');
     constE(k)=0*max(diag(t2));
%      constE(k)=gamma(k)*F_ini(:,k)'*G(:,:,k)'*G(:,:,k)*F_ini(:,k);
    C(:,:,k)=[constE(k)*eye(M)-G(:,:,k)*PHI(:,:,k)*G(:,:,k)'    -G(:,:,k)*PHI(:,:,k)*H(:,k);...
              -(G(:,:,k)*PHI(:,:,k)*H(:,k))'             0 ];
    R(:,:,k)=[G(:,:,k)*PHI(:,:,k)*PHI(:,:,k)*G(:,:,k)'   G(:,:,k)*PHI(:,:,k)*PHI(:,:,k)*H(:,k);...
              (G(:,:,k)*PHI(:,:,k)*PHI(:,:,k)*H(:,k))'   0 ];
    [T1,T2]=eig(-(H_error(k)^2+G_error(k)^2*M)*PHI(:,:,k));
%     y(k)=max([max(diag(T2)) 0]);
    if max(real(diag(T2)))>=0
        y(k)=max(real(diag(T2)));
    else
        y(k)=0;
    end
end

cvx_solver mosek
cvx_save_prefs

cvx_begin quiet
    variable E_hat(M+1,M+1) hermitian
    variable x(K) 
    variable SINR_residual(K) 

    
    expressions    LMI_S(M+1,M+1,K) LMI_IN(M+K,M+K,K);
    
 
 
    for k=1:K
        constraint_1(k)=(H_error(k)^2+G_error(k)^2*M)*trace(PHI(:,:,k))-sqrt(2*log(1/prob))*x(k)...
                        +log(prob)*y(k)-trace(C(:,:,k)*E_hat) -SINR_residual(k) -noise_maxpower...
                        +M*constE(k)+H(:,k)'*PHI(:,:,k)*H(:,k);
        constraint_2(k)=(H_error(k)^2+G_error(k)^2*M)^2*norm(PHI(:,:,k),'fro')^2 ...
                        +2*(H_error(k)^2+G_error(k)^2*M)*(trace(R(:,:,k)*E_hat)+H(:,k)'*PHI(:,:,k)*PHI(:,:,k)*H(:,k))...
                        -(2*real(x_ini(k)*x(k))-x_ini(k)^2);
      
    end 
    
    maximize 0*sum(real(constraint_1))+sum(SINR_residual)
  
    subject to
    
         real(constraint_1)>=0;
         real(constraint_2)<=0;
         E_hat  == hermitian_semidefinite(M+1);
         diag(E_hat) ==1;  
         SINR_residual>=0;
         x>=0;
     
     
cvx_end
if cvx_status(1)=='S'  ||  cvx_status(3)=='a'
    flag=1;

[t1,t2]=eig(E_hat);
location=find( abs(diag(t2))>10^(-6));
if size(location,1)==1
    e_hat=t1(:,location)*t2(location,location)^(1/2);
else
    for i=1:1000
        flag_2=1;
        b1(:,i)=t1*t2^(1/2)*sqrt(1/2)*(randn(M+1,1) + sqrt(-1)*  randn(M+1,1));
        b2(:,i)=exp(1j*angle(b1(:,i)/b1(M+1,i)));
        %%%%%  Obj value  %%%%%
        for k=1:K
            Obj(k,i)=b2(:,i)'*Matrix(:,:,k)*b2(:,i)+const(k);
            if Obj(k,i)<0
                flag_2=0;
                break;
            end
        end
        if flag_2==0
           Obj(:,i)=0;
        end    
    end
    Obj_sum=sum(Obj,1);
    [X,location]=max(Obj_sum);
%     [value location]=min(Obj);
    e_hat=b2(:,location);
end
e_opt=exp(1j*angle(e_hat(1:M)/e_hat(M+1)));

else
    flag=0;
    e_opt=ones(M,1);
end

end