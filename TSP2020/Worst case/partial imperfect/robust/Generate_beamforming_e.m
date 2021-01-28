function [e_opt,flag] = Generate_beamforming_e(N, M, K, H, G, G_error,...
                F_ini, e_ini, noise_maxpower, trans_maxpower, rate_min)

cvx_solver mosek
cvx_save_prefs

cvx_begin quiet
    variable e(M,1) complex
    variable z(K) 
    variable epselo(M)
    variable SINR_residual(K)
    variable relax_scaler_S(K)
    variable relax_scaler_IN(K)
    
    expressions    LMI_S(M*N+1,M*N+1,K) LMI_IN(M+K,M+K,K);
    
    E=diag(e);  E_ini=diag(e_ini);
    F=F_ini;
    for k=1:K
        %%%%%  Generate  LMI_S  %%%%%
        X(:,:,k)=(kron(F_ini(:,k)*F(:,k)',conj(e_ini)*e.'))...
                 +(kron(F_ini(:,k)*F(:,k)',conj(e_ini)*e.'))'...
                 -kron(F_ini(:,k)*F_ini(:,k)',conj(e_ini)*e_ini.');       
             
        x(:,k)=vec(e*(H(:,k)'+e_ini'*G(:,:,k))*F_ini(:,k)*F(:,k)')...
               +vec(e_ini*(H(:,k)'+e'*G(:,:,k))*F(:,k)*F_ini(:,k)')...
               -vec(e_ini*(H(:,k)'+e_ini'*G(:,:,k))*F_ini(:,k)*F_ini(:,k)');
           
        c(k)=2*real(  (H(:,k)'+e_ini'*G(:,:,k))*F_ini(:,k)*F(:,k)'*(H(:,k)+G(:,:,k)'*e)   )...
             -(H(:,k)'+e_ini'*G(:,:,k))*F_ini(:,k)*F_ini(:,k)'*(H(:,k)+G(:,:,k)'*e_ini);  
    
        constant(k)=c(k)-z(k)*(2^(rate_min)-1)-SINR_residual(k)-relax_scaler_S(k)*G_error(k)^2;

        LMI_S(:,:,k)=[relax_scaler_S(k)*eye(M*N)+conj(X(:,:,k))    x(:,k);...
                      x(:,k)'                                      constant(k)];
         

                   
         %%%%%  Generate  LMI_IN 1 %%%%%
         LMI_IN(M+K,M+K,K);
         F_k=[F(:,1:k-1)  F(:,k+1:K) ];
         temp_0=z(k)-noise_maxpower;
         temp_1=(H(:,k)'+e'*G(:,:,k))*F_k;
         A(:,:,k)=[temp_0       temp_1;...       
                   temp_1'      eye(K-1)]; 
         B(:,:,k)=[e zeros(M,K-1)];
         C(:,:,k)=[zeros(N,1)  F_k];
                       
         LMI_IN(:,:,k)=[A(:,:,k)-relax_scaler_IN(k)*C(:,:,k)'*C(:,:,k)   G_error(k)*B(:,:,k)';...
                        G_error(k)*B(:,:,k)                              relax_scaler_IN(k)*eye(M)];
                    
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
%          LMI_IN(:,:,k)=[A(:,:,k)-relax_scaler_IN(k)*C(:,:,k)'*C(:,:,k)   G_error(k)*B(:,:,k)';...
%                         G_error(k)*B(:,:,k)                              relax_scaler_IN(k)*eye(M*N)];
    end 
    
    maximize -0*sum(epselo)+sum(SINR_residual)
  
    subject to
         for k=1:K
             LMI_S(:,:,k)  == hermitian_semidefinite(M*N+1);
%              LMI_IN(:,:,k) == hermitian_semidefinite(N*M+K);
             LMI_IN(:,:,k) == hermitian_semidefinite(M+K);
             relax_scaler_IN(k)>0;
             relax_scaler_S(k)>0;
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