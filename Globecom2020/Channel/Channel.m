
close all;
clear all;clc;
N          = 32;            % array number of BS
M          = 32;            % array number of IRS
K          = 4;            % number of users in each group


%%%%% Large scale path loss
PL_0=10^(-30/10); %dB the channel gain at the reference distance
x_bs=0;
y_bs=0;
x_irs=50;
y_irs=10;

d_BI=sqrt((x_irs-x_bs)^2+(y_irs-y_bs)^2); %m distance from the BS to IRS
pathloss_BI=sqrt(PL_0*(d_BI)^(-2.2));    % Large-scale pass loss from the BS to the IRS


x_user(1)=70; y_user(1)=0;
x_user(2)=65; y_user(2)=-5;
x_user(3)=70; y_user(3)=5;
x_user(4)=73; y_user(4)=1;
for k=1:K
    d_BU(k)=sqrt((x_bs-x_user(k))^2+(y_bs-y_user(k))^2);%sqrt(d^2+d_v^2);  %m distance from the BS to the users
    d_IU(k)=sqrt((x_irs-x_user(k))^2+(y_irs-y_user(k))^2);  %m distance from the IRS to the users
    pathloss_BU(k)=sqrt(PL_0*(d_BU(k))^(-4));  % Large-scale pass loss from the BS to the users
    pathloss_IU(k)=sqrt(PL_0*(d_IU(k))^(-2));  % Large-scale pass loss from the IRS to the users
end

%% Simulation loop %%%%%%%%%%%%%%%%%%%%00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000%%%%%%%%%%%%%%%%%%%%%%%%%%
num_loop = 500; 

for loop = 1 : num_loop
    T1=cputime;
    H_dr_temp=sqrt(1/2)*(randn(M,N) + sqrt(-1)*  randn(M,N)); % small scale pass loss from the BS to the IRS
    H_dr_all(:,:,loop)=pathloss_BI*H_dr_temp;

    for k=1:K
        H_d_temp=sqrt(1/2)*(randn(N,1) + sqrt(-1)*  randn(N,1)); % small scale pass loss from the BS to the user
        H_d_all(:,k,loop)=pathloss_BU(k)*H_d_temp;
        H_r_temp=sqrt(1/2)*(randn(M,1) + sqrt(-1)*  randn(M,1)); % small scale pass loss from the IRS to the users
        H_r_all(:,k,loop)=pathloss_IU(k)*H_r_temp;
        G_r_all(:,:,k,loop)=diag(H_r_all(:,k,loop)')*H_dr_all(:,:,loop);
    end 
    
end




save('H_d_all','H_d_all');
save('G_r_all','G_r_all');