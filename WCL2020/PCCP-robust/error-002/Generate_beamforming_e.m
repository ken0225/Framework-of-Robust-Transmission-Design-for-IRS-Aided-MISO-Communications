function [e_opt,flag] = Generate_beamforming_e(N, M,  K, H_dr,H_d,H_r,H_error,...
                   F_ini, e_ini,  noise_maxpower, trans_maxpower,rate_min,n)
e_total(:,1)=e_ini;
scaler=10;
step=5;
scaler_max=max(200-n*20,30);
for n_outter=1:10
    scaler=10;
e_inner= sqrt(1/2)*(randn(M,1)+sqrt(-1)*randn(M,1));
e_inner=exp(1j.*angle(e_inner));
for temp=1:50
    
    cvx_solver mosek
    cvx_save_prefs

    cvx_begin quiet
        variable e(M,1) complex
        variable z(K) 
        variable epselo(2*M)
        variable SINR_residual(K)
        variable relax_scaler_S(K)
        variable relax_scaler_IN(K)

    expressions    LMI_S(M+1,M+1,K) LMI_IN(M+K,M+K,K)...
                   X(M,M,K) x(M,K) C(K) constrant(K)...
                   temp_0(K) temp_1(1,K) temp_2(M,K);

        E=diag(e);  E_ini=diag(e_ini); E_inner=diag(e_inner);
        F=F_ini;
        for k=1:K
            %%%%%  Generate  LMI_S  %%%%%
            X(:,:,k)=E*H_dr*F(:,k)*F_ini(:,k)'*H_dr'*E_ini'...
                     +(E*H_dr*F(:,k)*F_ini(:,k)'*H_dr'*E_ini')'...
                     -E_ini*H_dr*F_ini(:,k)*F_ini(:,k)'*H_dr'*E_ini';
            x(:,k)=E*H_dr*F(:,k)*F_ini(:,k)'*H_d(:,k)...
                 +E_ini*H_dr*F_ini(:,k)*F(:,k)'*H_d(:,k)...
                   -E_ini*H_dr*F_ini(:,k)*F_ini(:,k)'*H_d(:,k);
            c(k)=H_d(:,k)'*(F(:,k)*F_ini(:,k)'+F_ini(:,k)*F(:,k)'...
                 -F_ini(:,k)*F_ini(:,k)')*H_d(:,k);
            constrant(k)=real(H_r(:,k)'*X(:,:,k)*H_r(:,k)+x(:,k)'*H_r(:,k)+H_r(:,k)'*x(:,k)+c(k))...
                         -z(k)*rate_min-SINR_residual(k)-relax_scaler_S(k)*H_error(k)^2;
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

        maximize -scaler*sum(epselo)+sum(SINR_residual)

        subject to
             for k=1:K
                 LMI_S(:,:,k)  == hermitian_semidefinite(M+1);
                 LMI_IN(:,:,k) == hermitian_semidefinite(M+K);
                 relax_scaler_IN(k)>0;
                 relax_scaler_S(k)>0;
                 z(k)>=noise_maxpower;             
             end
             for m=1:M
                 norm(e(m),2)<=sqrt(1+epselo(m));
                 e_inner(m)'*e_inner(m)-2*real(e_inner(m)'*e(m))<=epselo(M+m)-1;
             end
             epselo>=0;

            SINR_residual>=0;

    cvx_end
    
    if cvx_status(1)=='S' || cvx_status(3)=='a'
        e_inner=e;
        a=[scaler+step,scaler_max];
        obj(temp+1)=-scaler*sum(epselo)+sum(SINR_residual);
        scaler=min(a);
        e_total(:,temp+1)=e;
        test_1(temp)=norm(e_total(:,temp+1)-e_total(:,temp),2)^2;
        test_2(temp)=sum(epselo);
    
        if test_1(temp)<=0.001 && test_2(temp)<=0.001
            break;
        end
    end
end
if temp<50
    break;
end

end


e_opt=e;

if cvx_status(1)=='S' ||  cvx_status(3)=='a'
    flag=1;

else
    flag=0;
end
end