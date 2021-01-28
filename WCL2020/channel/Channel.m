
close all;
clear all;clc;
N          = 6;            % array number of BS
M          = 16;            % array number of IRS
K          = 4;            % number of users in each group
Rican      = 5;
Rican_BU   = 5;
Rican_BI   = 5;
Rican_IU   = 10;

%%%%% Large scale path loss
PL_0=10^(-30/10); %dB the path loss at the reference distance, it can be also choosen 
x_bs=0;
y_bs=0;
x_irs=50;
y_irs=10;

d_BI=sqrt((x_irs-x_bs)^2+(y_irs-y_bs)^2); %m distance from the BS to IRS
pathloss_BI=sqrt(PL_0*(d_BI)^(-2.2));    % Large-scale pass loss from the BS to the IRS




%% Simulation loop %%%%%%%%%%%%%%%%%%%%00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000%%%%%%%%%%%%%%%%%%%%%%%%%%
num_loop = 4; 

for i_outer = 1 : num_loop
    %%%%%% location of users  %%%%%
    radius=5;
    x0=70;
    y0=0;
    user_radius=radius*sqrt(rand(K,1));
    seta=2*pi*rand(K,1);
    x_user=x0+user_radius.*cos(seta);
    y_user=y0+user_radius.*sin(seta);
    %%%%%  end   %%%%%
    %%%%%% LOS of BS-user  %%%%%
    N_t = linspace(0,N-1,N).';
    theta_AOD_BU=pi/3;
    theta_AOA_BU=0;
    steering_AOD_BU   = exp(1j*pi*N_t.*sin(theta_AOD_BU));
    steering_AOA_BU   = exp(1j*pi*1*sin(theta_AOA_BU));
    H_d_LOS=(steering_AOA_BU*steering_AOD_BU')';
    %%%%%  end   %%%%%
    %%%%%% LOS of BS-IRS  %%%%%
    N_M = linspace(0,M-1,M).';
    theta_AOA_BI=atan((x_irs-x_bs)/(y_irs-y_bs));
    theta_AOD_BI=pi/2-theta_AOA_BI;
    pha_AOA_BI=pi/6;
    steering_AOD_BI   = exp(1j*pi*N_t.*sin(theta_AOD_BI));
    for m=1:M
        steering_AOA_BI(m,1)   = exp(1j*pi* (floor(m/4)*sin(pha_AOA_BI)*sin(theta_AOA_BI)...
                                 +(m-floor(m/4)*4)*sin(pha_AOA_BI)*cos(theta_AOA_BI)) );
    end
    H_dr_LOS=steering_AOA_BI*steering_AOD_BI';
    %%%%%  end   %%%%%
    %%%%%% LOS of IRS-user  %%%%%
    theta_AOA_IU=atan((x_irs-x_user)./(y_irs-y_user));
    theta_AOD_IU=pi/2-theta_AOA_IU;
    pha_AOD_IU=atan(-10./sqrt((x_irs-x_user).^2+(y_irs-y_user).^2)); %10 means height
    

    %%%%%  end   %%%%%
    d_BU=sqrt((x_bs-x_user).^2+(y_bs-y_user).^2);%sqrt(d^2+d_v^2);  %m distance from the BS to the users
    d_IU=sqrt((x_irs-x_user).^2+(y_irs-y_user).^2);  %m distance from the IRS to the users
    pathloss_BU=sqrt(PL_0*(d_BU).^(-4));  % Large-scale pass loss from the BS to the users
    pathloss_IU=sqrt(PL_0*(d_IU).^(-2));  % Large-scale pass loss from the IRS to the users
    
    for i_inner=1:50
        loop=(i_outer-1)*50+i_inner;
        T1=cputime;
        H_dr_NLOS=sqrt(1/2)*(randn(M,N) + sqrt(-1)*  randn(M,N)); % small scale pass loss from the BS to the IRS
        H_BI_ALL(:,:,loop)=pathloss_BI*(sqrt(Rican_BI/(1+Rican_BI))*H_dr_LOS+sqrt(1/(1+Rican_BI))*H_dr_NLOS);

        for k=1:K
            H_d_NLOS=sqrt(1/2)*(randn(N,1) + sqrt(-1)*  randn(N,1)); % small scale pass loss from the BS to the user
            H_BU_ALL(:,k,loop)=pathloss_BU(k,1)*(sqrt(Rican_BU/(1+Rican_BU))*H_d_LOS+sqrt(1/(1+Rican_BU))*H_d_NLOS);

            for m=1:M
                steering_AOD_IU(m,1)   = exp(1j*pi* (floor(m/4)*sin(pha_AOD_IU(k,1))*sin(theta_AOD_IU(k,1))...
                    +(m-floor(m/4)*4)*sin(pha_AOD_IU(k,1))*cos(theta_AOD_IU(k,1))) );
            end
            steering_AOA_IU   = exp(1j*pi*1.*sin(theta_AOA_IU(k,1)));
            H_r_LOS(:,k)=steering_AOD_IU*steering_AOA_IU';
            H_r_NLOS=sqrt(1/2)*(randn(M,1) + sqrt(-1)*  randn(M,1)); % small scale pass loss from the IRS to the users
            H_IU_ALL(:,k,loop)=pathloss_IU(k,1)*(sqrt(Rican_IU/(1+Rican_IU))*H_r_LOS(:,k)+sqrt(1/(1+Rican_IU))*H_r_NLOS);
            H_total(:,:,k,loop)=[diag(H_IU_ALL(:,k,loop))'*H_BI_ALL(:,:,loop); H_BU_ALL(:,k,loop)'];

        end 
   
    end
end


save('H_BU_ALL','H_BU_ALL');
save('H_IU_ALL','H_IU_ALL');
save('H_BI_ALL','H_BI_ALL');
save('H_total','H_total');