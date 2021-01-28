%By Zhou Gui
%From 2019-10-19 to 
close all;
clear all;clc;
warning('off');
rand('twister',mod(floor(now*8640000),2^31-1));
%% Parameters Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% '1' stands for source-relays; '2' stands for relays-destination
N          = 8;            % array number of BS
M          = 8;            % array number of IRS
K          = 4;            % number of users 

SNR_dB     = 10;     % dBW
%%%%% noise
N0=10^((-174-30) / 10); %-174dBm
B=10^7; %10MHz
noise_maxpower_original   = N0*B;            % % W
noise_maxpower_original   = 10^((-90) / 10);            % % W
%%%%% end
% noise_maxpower   = 1;            % % W
trans_maxpower_all =0; % 

error         = 0.02; 
rate_min_dB   = [3:1:6]  ;   %bit/s/hz


%% Simulation loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_loop = 100; 
load('H_dr_all');
load('H_d_all');
load('H_r_all');
% for i_p     = 1 : length(rate_min_dB)
Power=zeros(100,length(rate_min_dB));
Rate=zeros(K,100,length(rate_min_dB));
for loop = 1 : num_loop
    outerflag=1;   
    T1=cputime;          
    H_dr=H_BI_ALL(1:M,1:N,loop);
    H_d=H_BU_ALL(1:N,1:K,loop)/sqrt(noise_maxpower_original);
    H_r=H_IU_ALL(1:M,1:K,loop)/sqrt(noise_maxpower_original);
    noise_maxpower=1;
    for k=1:K
        H_error(k)=error*norm(H_r(:,k),2);
 
    end
%%  For different SNR  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('  loop |  num_J  |  SNR  |  i  |  trans_SNR | relay_SNR \n');
for i_p     = 1 : length(rate_min_dB)
 
    t0=cputime;
    trans_maxpower=100 ;%trans_maxpower_all(i_p);
    sinr_min   = 2^(rate_min_dB (i_p))-1;%10^(rate_min_dB (i_p)  / 10); 
    %%%%%  Initialize F and e beamforming  %%%%%
     e_ini=ones(M,1);

     F_ini=ones(N,K)*sqrt(trans_maxpower/(N*K));
    F(:,:,1)=full(F_ini);
    e(:,1)=e_ini;

    num_iterative = 15000;
    for n  = 1 : num_iterative
       %%%%%  Optimize F  %%%%%
        [F_1,power_sub,innerflag] = Generate_beamforming_F(N, M,  K,...
                          H_dr,H_d,H_r,H_error, F(:,:,n), e(:,n),...
                          noise_maxpower, trans_maxpower,sinr_min);
        power(n+1)=power_sub;
        F(:,:,n+1)=F_1;
        if innerflag==0
            outerflag=0;
            break;
        end     

        %%%%%  Optimize e  %%%%%
        [e_1,flag_e] = Generate_beamforming_e(N, M,  K,...
                        H_dr,H_d,H_r,H_error, F(:,:,n+1), e(:,n),...
                        noise_maxpower, trans_maxpower,sinr_min,n);
        e(:,n+1)=e_1;
        if flag_e==0
            outerflag=0;
            break;
        end   
       
        %%%%%  stop criterion  %%%%%
       
    
        fprintf('   %g  |  %g  |  %g  \n',loop, rate_min_dB(i_p), n);
        if abs(power(n+1)-power(n))<10^(-3) 
            break;
        end
        
    end
    if outerflag==0
        break;
    end
    %%%%%  Generate the achievable rate of each user  %%%%%
      F_temp=F_1;
    E=diag(e_1);
    flag=ones(1000,1);
    for i_loop=1:1000
      for k=1:K
          temp=sqrt(1/2)*(randn(M,1) + sqrt(-1)* randn(M,1));
          H_r_error_channel(:,k)=H_error(k)*temp/norm(temp,2)^2*rand(1);
          H_prac=H_r(:,k)+H_r_error_channel(:,k);
          y(k,1)=norm((H_d(:,k)'+H_prac'*E*H_dr)*F_temp(:,k),2)^2;
          z_ini(k,1)=norm((H_d(:,k)'+H_prac'*E*H_dr)*F_temp,2)^2 ...
                 -y(k,1)+noise_maxpower;
          RATE(k,i_loop)=log2(1+y(k,1)/z_ini(k,1));
          if  RATE(k,i_loop) < rate_min_dB 
              flag((i_loop-1)*K+k)=0;
          end
      end
    end
    Rate(:,i_p,loop)=sum(RATE,2)/i_loop;
    Rate_ratio(i_p,loop)=sum(flag~=0)/length(flag);
     %%%%%  End  %%%%%
   
    Power(loop,i_p)=real(power(n+1));

    t2=cputime;
    CPU_Time(loop,i_p)=t2-t0;
    iteration(loop,i_p)=n;

end
    save('Power_e1','Power');
    save('Rate_e1','Rate');
%     save('iteration','iteration');
    T2=cputime;
    E1=T2-T1;
end
a=1;
    
