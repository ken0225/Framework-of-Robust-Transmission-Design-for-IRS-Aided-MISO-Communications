function [F_opt,power_opt,x_opt,flag] = Generate_beamforming_F(N, M, K, H, G, G_error,...
              F_ini, e_ini, prob, noise_maxpower, trans_maxpower, rate_min)


cvx_solver mosek
cvx_save_prefs

cvx_begin quiet
    variable Gamma(N,N,K) hermitian
    variable x(K) 
    variable y(K)
    expressions    PHI(N,N,K)  constant(K) constraint_1(K)...
                   constraint_2(K)  constraint_3(N,N,K);
    
    for k=1:K
        
        PHI(:,:,k)=2^(rate_min)/(2^(rate_min)-1)*Gamma(:,:,k)-sum(Gamma,3);
        constant(k)=(H(:,k)+G(:,:,k)'*e_ini)' * PHI(:,:,k) * (H(:,k)+G(:,:,k)'*e_ini)...
                    - noise_maxpower;
        
        constraint_1(k)=G_error(k)^2*M*trace(PHI(:,:,k))-sqrt(2*log(1/prob))*x(k)...
                        +log(prob)*y(k)+constant(k);
        element_1=G_error(k)^2*M*vec(PHI(:,:,k));
        element_2=sqrt(2*M)*G_error(k)*PHI(:,:,k)*(H(:,k)+G(:,:,k)'*e_ini);
        constraint_2(k)=norm([element_1;element_2],2);
        constraint_3(:,:,k)=y(k)*eye(N,N)+G_error(k)^2*M*PHI(:,:,k);
    end
    
    minimize trace(sum(Gamma,3))
  
    subject to
    
         real(constraint_1)>=0;
         real(constraint_2-x)<=0;
         y>=0;
    
         for k=1:K
             constraint_3(:,:,k)  == hermitian_semidefinite(N);
             Gamma(:,:,k) == hermitian_semidefinite(N);
         end
     
cvx_end

if cvx_status(1)=='S' || cvx_status(3)=='a' 
    flag=1;

   
   %%%%%  test rate constraints  %%%%%
for k=1:K
    Matrix(:,:,k)=(H(:,k)+G(:,:,k)'*e_ini)*(H(:,k)+G(:,:,k)'*e_ini)';
    PHI(:,:,k)=2^(rate_min)/(2^(rate_min)-1)*Gamma(:,:,k)-sum(Gamma,3);
    Obj_original(k)=trace(PHI(:,:,k)*Matrix(:,:,k))-noise_maxpower;
end
%%%%%          end            %%%%%

for k=1:K
    [t1,t2]=eig(Gamma(:,:,k));
    [X,location]=max( abs(diag(t2)) );
    if length(location)==1
        F_opt(:,k)=t1(:,location)*t2(location,location)^(1/2);
    else
        for i=1:500
            b2(:,i)=t1*t2^(1/2)*sqrt(1/2)*(randn(N,1) + sqrt(-1)*  randn(N,1));
            Obj(i)=norm(b2(:,i),2);
        end
        [value locat]=min(Obj);
        F_opt(:,k)=b2(:,locat);
    end
end
power_opt=trace(F_opt*F_opt');


%%%%%  test rate constraints  %%%%%
for k=1:K
    PHI(:,:,k)=2^(rate_min)/(2^(rate_min)-1)*F_opt(:,k)*F_opt(:,k)'-F_opt*F_opt';
    Obj_new(k)=trace(PHI(:,:,k)*Matrix(:,:,k))-noise_maxpower;

    element_1=G_error(k)^2*M*vec(PHI(:,:,k));
    element_2=sqrt(2*M)*G_error(k)*PHI(:,:,k)*(H(:,k)+G(:,:,k)'*e_ini);
    x_opt(k)=norm([element_1;element_2],2);
end
%%%%%          end            %%%%%
A=1;
else
    flag=0;
    F_opt=ones(N,K);
    power_opt=0;
    x_opt=0;
end

end

