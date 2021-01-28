function [F,power_opt,flag] = Generate_beamforming_F(N, M,  K,...
          H_dr,H_d,H_r,H_error, F_ini, e_ini,noise_maxpower,...
          trans_maxpower,rate_min)


cvx_solver mosek
cvx_save_prefs

cvx_begin quiet
    variable F(N,K) complex
    variables z(K) 
    variable relax_scaler_S(K)
    variable relax_scaler_IN(K)
    expressions    LMI_S(M+1,M+1,K) LMI_IN(M+K,M+K,K);
    
    E=diag(e_ini);  E_ini=diag(e_ini);
    for k=1:K
        %%%%%  Generate  LMI_S  %%%%%
        X(:,:,k)=E*H_dr*F(:,k)*F_ini(:,k)'*H_dr'*E_ini' ...
                 +(E*H_dr*F(:,k)*F_ini(:,k)'*H_dr'*E_ini')'...
                 -E_ini*H_dr*F_ini(:,k)*F_ini(:,k)'*H_dr'*E_ini';
        x(:,k)=E*H_dr*F(:,k)*F_ini(:,k)'*H_d(:,k)...
             +E_ini*H_dr*F_ini(:,k)*F(:,k)'*H_d(:,k)...
               -E_ini*H_dr*F_ini(:,k)*F_ini(:,k)'*H_d(:,k);
        c(k)=H_d(:,k)'*(F(:,k)*F_ini(:,k)'+F_ini(:,k)*F(:,k)'...
             -F_ini(:,k)*F_ini(:,k)')*H_d(:,k);
        constrant(k)=real(H_r(:,k)'*X(:,:,k)*H_r(:,k)+x(:,k)'*H_r(:,k)+H_r(:,k)'*x(:,k)+c(k))...
                     -z(k)*rate_min-relax_scaler_S(k)*H_error(k)^2;
         LMI_S(:,:,k)=[relax_scaler_S(k)*eye(M)+X(:,:,k)    x(:,k)+X(:,:,k)'*H_r(:,k);...
                       x(:,k)'+H_r(:,k)'*X(:,:,k)         constrant(k)];
         
         %%%%%  Generate  LMI_IN  %%%%%
         F_k=[F(:,1:k-1)  F(:,k+1:K) ];
         temp_0=z(k)-noise_maxpower-relax_scaler_IN(k);
         temp_1=(H_d(:,k)'+H_r(:,k)'*E*H_dr)*F_k;
         temp_2=E*H_dr*F_k;
         LMI_IN(:,:,k)=[temp_0       temp_1           zeros(1,M);...
                        temp_1'      eye(K-1)         H_error(k)*temp_2';...
                        zeros(M,1)   H_error(k)*temp_2   relax_scaler_IN(k)*eye(M)];
    end 
    
    minimize norm(F,'fro')
  
    subject to
         for k=1:K
             LMI_S(:,:,k)  == hermitian_semidefinite(M+1);
             LMI_IN(:,:,k) == hermitian_semidefinite(M+K);
             relax_scaler_IN(k)>0;
             relax_scaler_S(k)>0;
             z(k)>=noise_maxpower;
         end
     
cvx_end

   power_opt=norm(F,'fro')^2;
   for k=1:K
      y(k,1)=norm((H_d(:,k)'+H_r(:,k)'*E_ini*H_dr)*F(:,k),2)^2;
      z_temp(k,1)=norm((H_d(:,k)'+H_r(:,k)'*E_ini*H_dr)*F,2)^2 ...
             -y(k,1)+noise_maxpower;
      sinr(k)=y(k,1)/z_temp(k,1);
   end
  obj=sum(log2(1+sinr));
%   OBJ=sum(log2(1+SINR));

if cvx_status(1)=='S' || cvx_status(3)=='a' 
    flag=1;

else
    flag=0;
end
end

