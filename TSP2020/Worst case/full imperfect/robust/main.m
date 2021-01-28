%By Zhou Gui
%From 2019-10-19 to 
close all;
clear all;clc;
warning('off');
rand('twister',mod(floor(now*8640000),2^31-1));
%% Parameters Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% '1' stands for source-relays; '2' stands for relays-destination
N          = 6;            % array number of BS
M          = 6;            % array number of IRS
K          = 2;            % number of users 

SNR_dB     = 10;     % dBW
%%%%% noise
N0=10^((-174-30) / 10); %-174dBm
B=10^7; %10MHz
% noise_maxpower_original   = N0*B;            % % W
noise_maxpower_original   = 10^((-80-30) / 10);             % % W
%%%%% end
% noise_maxpower   = 1;            % % W
trans_maxpower_all =0; % trans_power=1 适用速率(0.5-6.5);trans_power=2 适用速率(7.5);trans_power=6 适用速率(8.5);

error_outage_G       = 0.01; 
error_outage_H       = 0.02; 
rate_min_dB   = [1:0.5:4]  ;   %bit/s/hz
prob=0.05;
error_G        = sqrt( error_outage_G^2*chi2inv(1-prob,2*N*M)/2 );
error_H        = sqrt( error_outage_H^2*chi2inv(1-prob,2*N)/2 );
%% Simulation loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_loop = 200; 
load('H_d_all');
load('G_r_all');
load('Channel_error');
% for i_p     = 1 : length(rate_min_dB)
Power=zeros(num_loop,length(rate_min_dB));
Rate=zeros(K,length(rate_min_dB),num_loop);
for loop = 65 : num_loop
    outerflag=1;   
    T1=cputime;          
    H=H_d_all(1:N,1:K,loop)/sqrt(noise_maxpower_original);
    G=G_r_all(1:M,1:N,1:K,loop)/sqrt(noise_maxpower_original);
    noise_maxpower=1;
%     H=H_d_all(1:N,1:K,loop);
%     G=G_r_all(1:M,1:N,1:K,loop);
%     noise_maxpower=noise_maxpower_original;
    for k=1:K
        G_error(k)=error_G*norm(G(:,:,k),'fro');
        H_error(k)=error_H*norm(H(:,k),'fro');
    end
%%  For different SNR  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('  loop |  num_J  |  SNR  |  i  |  trans_SNR | relay_SNR \n');
for i_p     = 1 : length(rate_min_dB)
 
    t0=cputime;
    trans_maxpower=100 ;%trans_maxpower_all(i_p);
    rate_min   = rate_min_dB (i_p);%10^(rate_min_dB (i_p)  / 10); 
    %%%%%  Initialize F and e beamforming  %%%%%
    e_temp=randn(M,1) + sqrt(-1)*  randn(M,1);
    e_ini=exp(1j*angle(e_temp));
    e_ini=ones(M,1);

%     F_ini=ones(N,K)*sqrt(trans_maxpower/(N*K));
    F_ini=randn(N,K)*sqrt(trans_maxpower/(N*K));
    F(:,:,1)=full(F_ini);
    e(:,1)=e_ini;

    num_iterative = 15000;
    for n  = 1 : num_iterative
       %%%%%  Optimize F  %%%%%
        [F_1,power_sub,innerflag] = Generate_beamforming_F(N, M, K, H, G, H_error,...
            G_error, F(:,:,n), e(:,n), noise_maxpower, trans_maxpower, rate_min);
        power(n+1)=power_sub;
        F(:,:,n+1)=F_1;
        if innerflag==0
            outerflag=0;
            break;
        end     
       
        
        %%%%%  Optimize e  %%%%%
        [e_1,flag_e] = Generate_beamforming_e(N, M, K, H, G, H_error, G_error,...
                F(:,:,n+1), e(:,n), noise_maxpower, trans_maxpower, rate_min);
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
    flag=ones(1000,1);
    for i_loop=1:1000
      for k=1:K
          G_error_channel(:,:,k)=sqrt(G_error(k)^2/2)*(unifrnd(-1,1,M,N) + sqrt(-1)* unifrnd(-1,1,M,N));
          G_prac(:,:,k)=G(:,:,k)+G_error_channel(:,:,k);
          H_error_channel(:,k)=sqrt(H_error(k)^2/2)*(unifrnd(-1,1,N,1) + sqrt(-1)* unifrnd(-1,1,N,1));
          H_prac(:,k)=H(:,k)+H_error_channel(:,k);
          y(k,1)=norm((H_prac(:,k)'+e_1'*G_prac(:,:,k))*F_temp(:,k),2)^2;
          z_ini(k,1)=norm((H_prac(:,k)'+e_1'*G_prac(:,:,k))*F_temp,2)^2 ...
                 -y(k,1)+noise_maxpower;
          RATE(k,i_loop)=log2(1+y(k,1)/z_ini(k,1));
          if  RATE(k,i_loop) < rate_min
              flag(i_loop)=0;
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
    save('Power','Power');
    save('Rate','Rate');
    save('Rate_ratio','Rate_ratio');
    T2=cputime;
    E1=T2-T1;
end
a=1;
    
