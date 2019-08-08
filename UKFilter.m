%
clc;
clear all;
close all force;
addpath('./utility');
addpath(genpath('./utility/UKF/'));
addpath('/home/valerio/Programs/Matlab_utils/altmany-export_fig');

%% SETTINGS

% 1) movimento lungo asse y
% 2) movimento lungo asse x
% 3) movimento circolare
caso = 3;

num_filters = 500;

% input std
std_v     = 0.1; % [m/s]
std_omega = 0.1745; % [rad/s]

% Initial position covariance 
cov_x     = 1^2; % x,y position [m]
cov_theta = (pi/180)^2; % theta position [rad]

% Measure std
std_measure = 0.01; % [m]

% Sampling times
dt = 0.1;               % simulation time step
rfid_measure_dt = 0.1;  % measure frequency

l = 5;   % rfid semidistance

plot_setting.line_width = 7;
plot_setting.style = {'k','--g','-.r',':b'};
plot_setting.font_size = 25;
%% Parameters
% num_rfid = 2;
output.rms_x = [];
output.rms_y = [];
output.rms_th = [];

for num_rfid = 1:4
    clearvars -except std_measure cov_x cov_theta std_v std_omega caso plot_setting output num_rfid l rfid_measure_dt dt t_stop num_filters legenda
    % Initial position and covariance of the robot
    initialPose = [10   -10  pi/4]';

    
    initialCov = diag([cov_x,cov_x,cov_theta]);
    
    % Odometry noise
    add_odom_noise = true;
    
    % Plo params
    plot_dt  = 1/5;          % seconds between frame
    L = 0.4; % wheelbase
    
    
    RFID = [0, 0,  l, -l
        l, -l, 0,  0]; % Position of RFID
    RFID = RFID(:,1:num_rfid);
    
    % Measure noise
    add_measure_noise = true;
    R = (std_measure)^2*eye(size(RFID,2),size(RFID,2)); % measure covariance
    
    %% Define the path
    fixed_v = 1;
    switch caso
        case 1 % movimento lungo asse y
            fixed_omega = 0.0;
            initialPose = [0; 20; -pi/2];
            t_stop = 40; %300;  % seconds to simulate

        case 2
            fixed_omega = 0.0;
            initialPose = [-50; 2; 0];
            t_stop = 100; %300;  % seconds to simulate

        case 3
            fixed_omega = 0.1667;
            initialPose = [-6; 1; -pi/2];
            t_stop = 40; %300;  % seconds to simulate
        
        otherwise
            error("caso non implementato");
    end
    
    %% Initialize a unicycle-like Robot
    
    close all force
    tartufino = UnicycleBot(initialPose, dt);
    
    filters = {};
    for i=1:num_filters
        noisy_initial_pose = mvnrnd(initialPose,initialCov)';
        
        m0 = [noisy_initial_pose(:); 0 ; 0];
        P0 = zeros(5,5);
        P0(1:3,1:3) = initialCov;
        filters{i,1} = UKF2(m0,P0); %EKF(m0,P0);
        filters{i,2} = ['ukf',num2str(i)];
    end
    
    %% Init for loop
    simulationTime = 0;
    
    real_trace = [];   % traccia seguita veramente dal robot
    filters_trace = {};  % traccia data dai filtri
    
    % Initialize
    for i=1:size(filters,1)
        filters_trace{i} = [];
    end
    
    h_waitbar = waitbar(0, 'Simulating...');
    
    % Store info
    past_state = tartufino.Pose;
    real_trace(:,end+1) = tartufino.Pose;
    for i=1:size(filters,1)
        filters_trace{i}(1:5,end+1) = filters{i,1}.x;
        filters_trace{i}(6,end) = sqrt(filters{i,1}.P(1,1));
        filters_trace{i}(7,end) = sqrt(filters{i,1}.P(2,2));
        filters_trace{i}(8,end) = sqrt(filters{i,1}.P(3,3));
    end
    
    %% Simulation loop
    last_measure = simulationTime;
    while simulationTime < t_stop % if time is not up
        
        % Compute the command to send to the robot
        uCmd(1) = fixed_v;  % linear velocity
        % Feed forward
        uCmd(2) = fixed_omega;  % linear velocity

        past_state = tartufino.Pose;
        
        % Drive the real robot
        drive(tartufino, uCmd);
        
        % Compute the measure and predict
        cov_u = [std_v^2, 0; 0, std_omega^2];
        for i=1:size(filters,1)
            % Flag per aggiungere o meno rumore di input
            if add_odom_noise
                measure_uCmd = mvnrnd(uCmd,cov_u)';
            else
                measure_uCmd = uCmd;
            end
            filters{i,1}.Prediction(measure_uCmd,dt,cov_u);
        end
        
        % Measure correction
        if (simulationTime-last_measure)>=rfid_measure_dt
            
            % Get the measure
            ideal_measurement = rfidReadings2(tartufino.Pose, past_state,RFID);
            
            if ~isnan(ideal_measurement)
                for i=1:size(filters,1)
                    % Add measure noise flag!
                    if add_measure_noise
                        measurement = mvnrnd(ideal_measurement, R)';
                        filters{i,1}.UpdateRFID(measurement,R,RFID);
                    else
                        filters{i,1}.UpdateRFID(ideal_measurement,R,RFID);
                    end
                end
            end
            last_measure = simulationTime;
        end
        
        % keep track of the path
        real_trace(:,end+1) = tartufino.Pose;
        for i=1:size(filters,1)
            filters_trace{i}(1:5,end+1) = filters{i,1}.x;
            filters_trace{i}(6,end) = sqrt(filters{i,1}.P(1,1));
            filters_trace{i}(7,end) = sqrt(filters{i,1}.P(2,2));
            filters_trace{i}(8,end) = sqrt(filters{i,1}.P(3,3));
%             filters_trace{i}(8 + 1:num_rfid,end) = filters{i,1}.innovation(:);
        end
        % Update simulation time
        simulationTime = simulationTime + dt;
        
        h_waitbar = waitbar(simulationTime/t_stop);
    end
    close(h_waitbar);
    
    
    %% Plot delle tracce per avere idea qualitativa di quello che e' successo
    if false
        figure();
        ax = gca;
        hold on; axis equal; box; grid on;
        
        legenda = {};
        for i=1:size(filters,1)
            
            plot(ax,filters_trace{i}(1,1),filters_trace{i}(2,1),'.')
            plot(ax,filters_trace{i}(1,1:end),filters_trace{i}(2,1:end))%,'linewidth', 4);
            legenda{end+1} = 'estimated trajectory';
        end
        plot(ax,real_trace(1,:),real_trace(2,:),'r', 'linewidth', 4);
        legenda{end+1} = 'actual trajectory';
        plot(RFID(1,:),RFID(2,:),'ks','MarkerFaceColor','k','markersize',17);
        legenda{end+1} = 'RFID location';
        %legend(legenda(:));
        xlabel('x [m]')
        ylabel('y [m]')
        set(gca,'FontSize', 35)
        hold off
    end
    
    %% Compute RMSE
    n_filt = size(filters,1);
    n_stamp = numel(real_trace(1,:));
    error_x = zeros(n_filt, n_stamp);
    error_y = zeros(n_filt, n_stamp);
    error_th = zeros(n_filt, n_stamp);
    
    for i=1:n_filt
        error_x(i,:) = filters_trace{i}(1,:) - real_trace(1,:);
        error_y(i,:) = filters_trace{i}(2,:) - real_trace(2,:);
        delta_ang = (filters_trace{i}(3,:) - real_trace(3,:))*180/pi;
        while(sum(delta_ang<0)>0)
            delta_ang = delta_ang + 360;
        end
        error_th(i,:) = rem(delta_ang + 180, 360) - 180;
    end
    
    rms_x = sqrt(sum(error_x.^2,1)/n_stamp);
    rms_y = sqrt(sum(error_y.^2,1)/n_stamp);
    rms_th = sqrt(sum(error_th.^2,1)/n_stamp);
    
    times = dt*(0:(n_stamp-1));
    
    output.rms_x = [output.rms_x; rms_x];
    output.rms_y = [output.rms_y; rms_y];
    output.rms_th = [output.rms_th; rms_th];
end

%%
figure(1); clf;
subplot(3,1,1); hold on; box;
for i=1:4
    plot(times, output.rms_x(i,:), plot_setting.style{i}, 'linewidth', plot_setting.line_width);
end
legend('1 RFID','2 RFID','3 RFID','4 RFID','Orientation','horizontal');
ylabel('r_x [m]');
xlim([0,t_stop]);
set(gca,'FontSize', plot_setting.font_size)
subplot(3,1,2); hold on; box;
for i=1:4
plot(times, output.rms_y(i,:), plot_setting.style{i}, 'linewidth', plot_setting.line_width);
end
ylabel('r_y [m]');
xlim([0,t_stop]);
set(gca,'FontSize', plot_setting.font_size)
subplot(3,1,3); hold on; box;
for i=1:4
plot(times, output.rms_th(i,:), plot_setting.style{i}, 'linewidth', plot_setting.line_width);
end
ylabel('r_\theta [deg]');
xlabel('time [s]');
xlim([0,t_stop]);
set(gca,'FontSize', plot_setting.font_size)
%%

% %%
% times = dt*(0:(numel(filters_trace{i}(9,:))-1));
% figure(); hold on; box;
% plot(times, filters_trace{i}(9,:), 'linewidth', 4);
% plot(times, filters_trace{i}(10,:), 'linewidth', 4);
% xlabel("time [s]");
% ylabel("innovation [m]");
% legend('RFID 1','RFID 2');
% set(gca,'FontSize', plot_setting.font_size)
%
% %%
% font_size = 25;
% figure();
% for i=1:size(filters,1)
%     subplot(3,1,1);
%     hold on; box; axis tight;
%     plot(times, filters_trace{i}(1,:) - real_trace(1,:),'linewidth', 2)
%     plot(times, filters_trace{i}(6,:) ,'r','linewidth', 2)
%     plot(times, -filters_trace{i}(6,:) ,'r','linewidth', 2)
%     ylabel('e_x [m]');
% %     ylim([-4,4]);
%     %xlabel("time [s]");
%     set(gca,'FontSize', font_size);
%     subplot(3,1,2);
%     hold on; box; axis tight;
%     plot(times, filters_trace{i}(2,:) - real_trace(2,:),'linewidth', 2)
%     plot(times, filters_trace{i}(7,:) ,'r','linewidth', 2)
%     plot(times, -filters_trace{i}(7,:) ,'r','linewidth', 2)
%     ylabel('e_y [m]');
% %     ylim([-4,4]);
%     %xlabel("time [s]");
%     set(gca,'FontSize', font_size);
%     subplot(3,1,3);
%     hold on; box; axis tight;
%     delta_ang = (filters_trace{i}(3,:) - real_trace(3,:))*180/pi;
%     while(sum(delta_ang<0)>0)
%         delta_ang = delta_ang + 2*pi;
%     end
%     delta_ang = rem(delta_ang + 180, 360) - 180;
%     plot(times, delta_ang,'linewidth', 2)
%     plot(times, filters_trace{i}(8,:)*180/pi ,'r','linewidth', 2)
%     plot(times, -filters_trace{i}(8,:)*180/pi ,'r','linewidth', 2)
%     ylabel('e_\theta [deg]');
%     xlabel("time [s]");
% %     ylim([-60,60]);
%     set(gca,'FontSize', font_size);
% end



%% TRASH
% figure(2); clf; hold on; axis equal;
% quiver(initialPositionList(1,:),initialPositionList(2,:),cos(initialPositionList(3,:)),sin(initialPositionList(3,:)),0)
% quiver(failedPositionList(1,:),failedPositionList(2,:),cos(failedPositionList(3,:)),sin(failedPositionList(3,:)),0)
% plot(initialPositionList(1,:), initialPositionList(2,:), '.g','markersize',6)
% plot(failedPositionList(1,:), failedPositionList(2,:), '.r','markersize',6)
% plot(path(1,:),path(2,:),'--k','linewidth',4);
% plot(RFID(1,:),RFID(2,:),'s','MarkerFaceColor','k','markersize',14);
% plot(initialPose(1),initialPose(2),'gx','linewidth',4,'markersize',14)
