function [e_opt,flag] = Generate_beamforming_e(N, M, K, H, G, H_error, G_error,...
                F_ini, e_ini, noise_maxpower, trans_maxpower, rate_min)

            
%%  test
%     e=randn(M,1);
%     F=ones(N,K)*sqrt(trans_maxpower/(N*K));
%     for k=1:K
%         delta_G(:,:,k)=randn(M,N) + sqrt(-1)*  randn(M,N);
%         delta_H(:,k)=randn(N,1) + sqrt(-1)*  randn(N,1);
%         %%%%%  Generate  LMI_S  %%%%%
%         temp_D_ini = [kron(F_ini(:,k),1);  kron(F_ini(:,k),conj(e_ini))];
%         temp_D     = [kron(F(:,k),1);      kron(F(:,k),conj(e))];
%         D(:,:,k)   = temp_D_ini*temp_D';
%         Z(:,:,k)   = temp_D_ini*temp_D_ini';
%         
%         temp_1     = F(:,k)*F_ini(:,k)'*(H(:,k)+G(:,:,k)'*e_ini);
%         temp_2     = F_ini(:,k)*F(:,k)'*(H(:,k)+G(:,:,k)'*e);
%         temp_3     = F_ini(:,k)*F_ini(:,k)'*(H(:,k)+G(:,:,k)'*e_ini); 
%         d1(:,k)    = [temp_1; conj( vec(e*temp_1') )];
%         d2(:,k)    = [temp_2; conj( vec(e_ini*temp_2') )];
%         zz(:,k)     = [temp_3; conj( vec(e_ini*temp_3') )];
%                   
%         temp_c_ini = F_ini(:,k)'*(H(:,k)+G(:,:,k)'*e_ini);        
%         temp_c     = F(:,k)'*(H(:,k)+G(:,:,k)'*e);  
%         const_d(k) = temp_c_ini'*temp_c;
%         const_z(k) = temp_c_ini'*temp_c_ini;
%         
%         constant(k)= 2*real(const_d(k))-const_z(k);
%         
%         X(:,:,k)    = D(:,:,k)+D(:,:,k)'-Z(:,:,k);
%         x(:,k)      = d1(:,k)+d2(:,k)-zz(:,k);
%            
%     
%         vari(:,k)=[delta_H(:,k);conj(vec(delta_G(:,:,k)))];
%         rate_new(k)=vari(:,k)'*X(:,:,k)*vari(:,k)...
%                     +2*real( x(:,k)'*vari(:,k) )...
%                     +constant(k);
%                 
%         rate_old(k)=2*real( ((H(:,k)+delta_H(:,k))'+e_ini'*(G(:,:,k)+delta_G(:,:,k)))*F_ini(:,k)*F(:,k)'*((H(:,k)+delta_H(:,k))+(G(:,:,k)+delta_G(:,:,k))'*e) )...
%                     -((H(:,k)+delta_H(:,k))'+e_ini'*(G(:,:,k)+delta_G(:,:,k)))*F_ini(:,k)*F_ini(:,k)'*((H(:,k)+delta_H(:,k))+(G(:,:,k)+delta_G(:,:,k))'*e_ini);
%                 
%      
%     end

%%           
            
            
cvx_solver mosek
cvx_save_prefs

cvx_begin quiet
    variable e(M,1) complex
    variable z(K) 
    variable epselo(M)
    variable SINR_residual(K)
    variable relax_S_H(K)
    variable relax_S_G(K)
    variable relax_IN_H(K)
    variable relax_IN_G(K)
    
%     expressions    LMI_S(M*N+1,M*N+1,K) LMI_IN(N*M+K,N*M+K,K);
    expressions    LMI_S(N+M*N+1,N+M*N+1,K) LMI_IN(M+K,M+K,K);
    
    F=F_ini;
    for k=1:K
        %%%%%  Generate  LMI_S  %%%%%
        temp_D_ini = [kron(F_ini(:,k),1);  kron(F_ini(:,k),conj(e_ini))];
        temp_D     = [kron(F(:,k),1);      kron(F(:,k),conj(e))];
        D(:,:,k)   = temp_D_ini*temp_D';
        Z(:,:,k)   = temp_D_ini*temp_D_ini';
        
        temp_1     = F(:,k)*F_ini(:,k)'*(H(:,k)+G(:,:,k)'*e_ini);
        temp_2     = F_ini(:,k)*F(:,k)'*(H(:,k)+G(:,:,k)'*e);
        temp_3     = F_ini(:,k)*F_ini(:,k)'*(H(:,k)+G(:,:,k)'*e_ini); 
        d1(:,k)    = [temp_1; conj( vec(e*temp_1') )];
        d2(:,k)    = [temp_2; conj( vec(e_ini*temp_2') )];
        zz(:,k)    = [temp_3; conj( vec(e_ini*temp_3') )];
                  
        temp_c_ini = F_ini(:,k)'*(H(:,k)+G(:,:,k)'*e_ini);        
        temp_c     = F(:,k)'*(H(:,k)+G(:,:,k)'*e);  
        const_d(k) = temp_c_ini'*temp_c;
        const_z(k) = temp_c_ini'*temp_c_ini;
        
        constant(k)= 2*real(const_d(k))-const_z(k)-z(k)*(2^(rate_min)-1)-SINR_residual(k)...
                     -relax_S_H(k)*H_error(k)^2-relax_S_G(k)*G_error(k)^2;
        
        X(:,:,k)    = D(:,:,k)+D(:,:,k)'-Z(:,:,k)...
                      +[relax_S_H(k)*eye(N) zeros(N,N*M);...
                        zeros(N*M,N) relax_S_G(k)*eye(N*M)];
        x(:,k)      = d1(:,k)+d2(:,k)-zz(:,k);
        LMI_S(:,:,k)=[X(:,:,k)    x(:,k);...
                      x(:,k)'     constant(k)];
         

                   
         %%%%%  Generate  LMI_IN 1 %%%%%
%          LMI_IN(2*N+K,2*N+K,K);
%          F_k=[F(:,1:k-1)  F(:,k+1:K) ];
%          temp_0=z(k)-noise_maxpower-relax_IN_G(k)*M-relax_IN_H(k);
%          temp_1=(H(:,k)'+e'*G(:,:,k))*F_k;        
%          LMI_IN(:,:,k)=[temp_0       temp_1           zeros(1,N)            zeros(1,N);...
%                         temp_1'      eye(K-1)         G_error(k)*F_k'       H_error(k)*F_k';...
%                         zeros(N,1)   G_error(k)*F_k   relax_IN_G(k)*eye(N)  zeros(N,N);...
%                         zeros(N,1)   H_error(k)*F_k   zeros(N,N)            relax_IN_H(k)*eye(N)];
                    
                    

         LMI_IN(M+K,M+K,K);
         F_k=[F(:,1:k-1)  F(:,k+1:K) ];
         temp_0=z(k)-noise_maxpower;
         temp_1=(H(:,k)'+e'*G(:,:,k))*F_k;
         A(:,:,k)=[temp_0       temp_1;...       
                   temp_1'      eye(K-1)]; 
         B_G(:,:,k)=[e zeros(M,K-1)];
         B_H(:,:,k)=zeros(K,K);
         B_H(1,1,k)=1;
         C_GH(:,:,k)=[zeros(N,1)  F_k];
         temp_2=A(:,:,k)-(relax_IN_G(k)+relax_IN_H(k))*C_GH(:,:,k)'*C_GH(:,:,k);           
%          LMI_IN(:,:,k)=[temp_2                   G_error(k)*B_G(:,:,k)'   H_error(k)*B_H(:,:,k)';...
%                         G_error(k)*B_G(:,:,k)    relax_IN_G(k)*eye(M)     zeros(M,K);...
%                         H_error(k)*B_H(:,:,k)    zeros(K,M)               relax_IN_H(k)*eye(K)];
                    
         LMI_IN(:,:,k)=[temp_2                   G_error(k)*B_G(:,:,k)';...  
                        G_error(k)*B_G(:,:,k)    relax_IN_G(k)*eye(M)];...   
%                       
         %%%%%  Generate  LMI_IN 2 %%%%%
%        LMI_IN(2*K,2*K,K);
%          F_k=[F(:,1:k-1)  F(:,k+1:K) ];
%          temp_0=z(k)-noise_maxpower;
%          temp_1=(H(:,k)'+e'*G(:,:,k))*F_k;
%          A(:,:,k)=[temp_0       temp_1;...       
%                    temp_1'      eye(K-1)]; 
%          B(:,:,k)=[zeros(M*N,1)  kron(conj(e),F_k)];
%          C(:,:,k)=zeros(K,K);
%          C(1,1,k)=1;
%                        
%          LMI_IN(:,:,k)=[A(:,:,k)-relax_IN(k)*C(:,:,k)'*C(:,:,k)   G_error(k)*B(:,:,k)';...
%                         G_error(k)*B(:,:,k)                              relax_IN(k)*eye(M*N)];
    end 
    
    maximize -0*sum(epselo)+sum(SINR_residual)
  
    subject to
         for k=1:K
             LMI_S(:,:,k)  == hermitian_semidefinite(N+M*N+1);
%              LMI_IN(:,:,k) == hermitian_semidefinite(N*M+K);
             LMI_IN(:,:,k) == hermitian_semidefinite(M+K);
             relax_IN_H(k)>0;
             relax_IN_G(k)>0;
             relax_S_H(k)>0;
             relax_S_G(k)>0;
             z(k)>=noise_maxpower;             
         end

         for m=1:M
             norm(e(m),2)<=sqrt(1);

         end
         
%          for m=1:M
%              norm(e(m),2)<=sqrt(1+epselo(m));
%              e_ini(m)'*e_ini(m)-2*real(e_ini(m)'*e(m))<=epselo(m)-1;
%          end
%          epselo>=0;
         
        SINR_residual>=0;
             
     
cvx_end

    e_opt=exp(1i*angle(e));

if cvx_status(1)=='S'  ||  cvx_status(3)=='a'
    flag=1;

else
    flag=0;
end

%% TEST
%%%%%  before
   for k=1:K
      y(k,1)=norm((H(:,k)'+e'*G(:,:,k))*F(:,k),2)^2;
      
      z_temp(k,1)=norm((H(:,k)'+e'*G(:,:,k))*F,2)^2 ...
             -y(k,1)+noise_maxpower;
      Rate(k)=log2(1+y(k,1)/z_temp(k,1));
   end
  obj=sum(Rate);
%   OBJ=sum(log2(1+SINR));

%%%%%  after
   for k=1:K
      y(k,1)=norm((H(:,k)'+e_opt'*G(:,:,k))*F(:,k),2)^2;
      
      z_temp(k,1)=norm((H(:,k)'+e_opt'*G(:,:,k))*F,2)^2 ...
             -y(k,1)+noise_maxpower;
      Rate_opt(k)=log2(1+y(k,1)/z_temp(k,1));
   end
  obj_opt=sum(Rate_opt);
%   OBJ=sum(log2(1+SINR));
end