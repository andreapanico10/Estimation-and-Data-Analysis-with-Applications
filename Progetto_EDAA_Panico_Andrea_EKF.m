%extended Kalman Filter EKF
clear all;
close all;
clc;

% first value = distance measeure
% second value = orientation measures
% 1 1 means AND
% 0 1 means ONLY orientation measures
% 1 0 means ONLY orientation measures
truth_table = [1 1];

% number of beacons evenly distributed around the work area 
n_beacons = 4;
beacons_positions = zeros(2,n_beacons);

% find n_beacons equidistant points on a circumference
radius = 17; % [m] circumference of radius = radius
circumference_lenght = 2*pi*radius;
distance_between_beacons = circumference_lenght/n_beacons; %[m] distance between two consecutives beacons 
alfa = (2*pi)/n_beacons; 

circle_center = [5,5]'; % center of the circumference formed by beacons
% Beacons coordinates in the cartesian plane needed by EKF approach
for i = 1: n_beacons
    beacons_positions(:,i) = circle_center + [radius*cos((i-1)*alfa); radius*sin((i-1)*alfa)]; %  [ [m] [m] ]
end
for mc_i = 1:1000
% Computation parameters
T = 100;          % Running time [sec]
Tc = 0.15;        % Sampling time [sec]
N = round(T/Tc);  % Sampling steps

sound_speed = 343; %[m/s] % propagation speed of an acoustic wave
beacon_range = sound_speed*Tc/2; %[m] max range reachable by a beacon (by round trip time)

% System inputs for both estimation approaches
% Robot trajectory

v = 0.75 * ones(N,1)';          % [m/s] velocity of the robot

omega = zeros(N,1)';

for i = 1:N
    omega(i) =  3.5*pi/180;  % angular speed  [rad/s]
end


% Process noise covariance matrix for EKF approach
sigma_x = 0.001; %[m]
sigma_y = 0.001; %[m]
sigma_theta = 0.1 * pi/180; %[rad]
Q = diag([sigma_x^2 sigma_y^2 sigma_theta^2]); % Covariance matrix 

% Measurements noise covariance matrix definition for EKF approach
% for hypothesys all the devices have got same standard deviation for
% distance measurement error and they have got same standard deviation for
% orientation measurement

sigma_distance = 0.3; %[m]
sigma_theta = 10*pi/180; %[rad]

% EKF initialization
X(:,1) = [5 -7.5 1 * pi/180]';            % Real initial position of the robot
x_pred = [5 -2.5 1 * pi/180]';               % Initial position of the estimation
x_hat = [5 -2.5 1 * pi/180]';                % Initial position of the correction
P = diag([0.5 0.5 (5*pi/180)^2]);
valid_measures= zeros(N,n_beacons*2)';    
x =             zeros(N,1)';
y =             zeros(N,1)';
theta =         zeros(N,1)';
zn =            zeros(N,n_beacons*2)';
z =             zeros(N,n_beacons*2)';
z_subtracted_full = zeros(N,n_beacons*2)';
est_error =     zeros(N,1)';
P1 = [];
P2 = [];
P3 = [];

for k = 1:N
    
    %real model
    
    X(:,k+1) = [X(1,k) + Tc*v(k)*cos(X(3,k));
                X(2,k) + Tc*v(k)*sin(X(3,k));
                X(3,k) + Tc*omega(k);] + sqrt(Q) * randn(3,1);
    
    % 3 components extraction
    
    x(:,k+1) = X(1,k+1);
    y(:,k+1) = X(2,k+1);
    theta(:,k+1) = X(3,k+1);

    % Measurement equation
    
    for beacon = 1:n_beacons
        % z = real output
        % distance measurement
        z(beacon,k+1) = sqrt((beacons_positions(1, beacon)-x(:,k))^2 + (beacons_positions(2, beacon)-y(:,k))^2);                
        % orientation measurement
        z(n_beacons + beacon,k+1) = atan2(beacons_positions(2, beacon) - y(:,k), beacons_positions(1, beacon) - x(:,k)) - theta(:,k);
    end
    
    
    % zn = Real measurement contaminated by noise
    if truth_table(1)  % Only distance measurements available
        R_distance = sigma_distance^2 * eye(n_beacons);
        R = R_distance;
    end
    if truth_table(2) % Only orientation measurements available
        R_theta = sigma_theta^2 * eye(n_beacons);
        R = R_theta;
    end
    
    % if both measures type are available, we have a block diagonal matrix
    if (truth_table(1) && truth_table(2))     
        R = blkdiag(R_distance,R_theta);
    end
    
    if (truth_table(1) && truth_table(2))
        zn(:,k+1) = z(:,k+1) + sqrt(R) * randn(n_beacons*2,1);
    elseif truth_table(1) == 0  
        zn((n_beacons+1):(2*n_beacons),k+1) = z((n_beacons+1):2*n_beacons,k+1) + sqrt(R) * randn(n_beacons,1);
    elseif truth_table(2) == 0  
        zn(1:n_beacons,k+1) = z(1:n_beacons,k+1) + sqrt(R) * randn(n_beacons,1);
    end

    
    % not all beacons are able to reach the vehicle, 
    % so set the measurements of the beacons out of range to 0
    
    for measure = 1:n_beacons
        flight_time = 2 * (z(measure,k+1)/sound_speed); % round trip time (2x)
        % hypothesis of non-synchronization of the beacons
        if flight_time < Tc 
            if truth_table(1)
                valid_measures(measure, k) = 1;
            end
            if truth_table(2)
                valid_measures(measure+n_beacons, k) = 1;
            end
        else % if beacon is not reachable, the measure is lost
            if truth_table(1)
                valid_measures(measure, k) = 0;
                zn(measure,k+1) = 0;
                z(measure,k+1) = 0;
            end

            if truth_table(2)
                valid_measures(measure+n_beacons, k) = 0;
                zn(measure+n_beacons,k+1) = 0;
                z(measure+n_beacons,k+1) = 0;
            end
        end    
    end
    
    % # of valid measures in order do initialize matrix R 
    n_valid_measures = 0;
    for i = 1:2*n_beacons
        if valid_measures(i,k) == 1
                % each sensor gives 1 correct measurement
                n_valid_measures = n_valid_measures+1;   
        end
    end
    
    % Redefine R in order to delete zero-row refer to the measures not available
    R = eye(n_valid_measures); 
    if (~truth_table(1) || ~truth_table(2))
        for i = 1:(n_valid_measures)
            if (truth_table(1) && ~truth_table(1))
                R(i,i) = sigma_distance^2;
            end
            if (truth_table(2) && ~truth_table(1))
                R(i,i) = sigma_theta^2;
            end
        end
    else
        for i = 1:(n_valid_measures/2)
            
                R(i,i) = sigma_distance^2;
                R(n_valid_measures/2 + i,n_valid_measures/2 + i) = sigma_theta^2;
            
        end
    end
    % EKF Prediction step

    F = [ 1 0 -Tc * v(k) * sin(x_hat(3,k));
          0 1  Tc * v(k) * cos(x_hat(3,k));
          0 0           1                 ];  
      
    L = eye(3);
    
    
    P_pred = F * P * F' + L * Q * L';                      % P(k+1|k)
    x_pred(:,k+1) = [x_hat(1,k) + Tc*v(k)*cos(x_hat(3,k)); % x(k+1|k)
                     x_hat(2,k) + Tc*v(k)*sin(x_hat(3,k));
                     x_hat(3,k) + Tc*omega(k);          ];
            
            
    % EKF Correction step
    
    H_full = zeros(n_beacons*2,3);
    % setting to 0 the components of H for which measures are not available
    for id_beacon = 1:n_beacons
        if (zn(id_beacon,k+1) == 0 && zn(n_beacons + id_beacon, k+1) == 0)
            if truth_table(1)
                H_full(id_beacon,1) = 0;             % x coordinate of Jacobian relative to distance
                H_full(id_beacon,2) = 0;             % y coordinate of Jacobian relative to distance
                z_subtracted_full(id_beacon, k+1) = 0;

            end
            if truth_table(2)
                H_full(n_beacons + id_beacon,1) = 0; % x coordinate of Jacobian relative to orientation
                H_full(n_beacons + id_beacon,2) = 0; % y coordinate of Jacobian relative to orientation
                H_full(n_beacons + id_beacon,3) = 0; % theta coordinate of Jacobian relative to orientation
                z_subtracted_full(n_beacons + id_beacon,k+1) = 0;

            end
            
        else
            if truth_table(1)
                H_full(id_beacon,1) = (x_pred(1,k+1) - beacons_positions(1, id_beacon)) / sqrt( (beacons_positions(1, id_beacon) - x_pred(1,k+1))^2 + (beacons_positions(2, id_beacon) - x_pred(2,k+1))^2  );
                H_full(id_beacon,2) = (x_pred(2,k+1) - beacons_positions(2, id_beacon)) / sqrt( (beacons_positions(1, id_beacon) - x_pred(1,k+1))^2 + (beacons_positions(2, id_beacon) - x_pred(2,k+1))^2  );
                z_subtracted_full(id_beacon, k+1) = sqrt((beacons_positions(1, id_beacon)-x_pred(1,k+1))^2 + (beacons_positions(2, id_beacon)-x_pred(2,k+1))^2);

            end
            if truth_table(2)
                H_full(n_beacons + id_beacon,1) = (-(x_pred(2,k+1) - beacons_positions(2, id_beacon)))/((x_pred(1,k+1) - beacons_positions(1, id_beacon))^2 + (x_pred(2,k+1) - beacons_positions(2, id_beacon))^2);
                H_full(n_beacons + id_beacon,2) = (x_pred(1,k+1) - beacons_positions(1, id_beacon))/((x_pred(1,k+1) - beacons_positions(1, id_beacon))^2 + (x_pred(2,k+1) - beacons_positions(2, id_beacon))^2);
                H_full(n_beacons + id_beacon,3) = -1;
               
                z_subtracted_full(n_beacons + id_beacon,k+1) = atan2(beacons_positions(2, id_beacon) - x_pred(2,k+1), beacons_positions(1, id_beacon) - x_pred(1,k+1)) - x_pred(3,k+1);
            end
        end 
    end
    
     H = zeros(n_valid_measures,3);
     counter = 1;
     for i = 1:n_beacons*2
         if valid_measures(i,k) == 1
            H(counter,:) = H_full(i,:); 
            counter = counter + 1;
         end
     end

    
    M = eye(n_valid_measures);
  
    z_subtracted = zeros(n_valid_measures,1);
    counter = 1;
    for i = 1:n_beacons*2
         if valid_measures(i,k) == 1
            z_subtracted(counter) = z_subtracted_full(i,k+1);  
            counter = counter +1;
         end
    end
     
    %z_n is zn with zero-row suppressed
    z_n = zeros(n_valid_measures,1);
    counter = 1;
    for i = 1:n_beacons*2
         if valid_measures(i,k) == 1
            z_n(counter) = zn(i,k+1);  
            counter = counter +1;
         end
    end
    
    W = P_pred * H'* inv(H * P_pred * H'+ M * R * M');                  % Gain W
    P = (eye(3) - W * H) * P_pred;                                      % P(k+1|k+1)
    x_hat(:,k+1) = x_pred(:,k+1) + W*(z_n - z_subtracted);              % x(k+1|k+1)
    
    % estimation errors 
    est_error(1,k+1) = X(1,k+1) - x_hat(1,k+1);
    est_error(2,k+1) = X(2,k+1) - x_hat(2,k+1);
    est_error(3,k+1) = X(3,k+1) - x_hat(3,k+1);
    
    tot_meas = size(est_error);

        
    %x position estimation error covariance
    P1 = [P1, P_pred(1,1) P(1,1)];
    
    %y position estimation error covariance
    P2 = [P2, P_pred(2,2) P(2,2)];
    
    %theta estimation error covariance
    P3 = [P3, P_pred(3,3) P(3,3)];
end

t = Tc * (1:N+1); % Time [s]

    est_err_sample_1(1,mc_i) = est_error(1, randi([1,tot_meas(2)],1));
    est_err_sample_2(2,mc_i) = est_error(2, randi([1,tot_meas(2)],1));
    est_err_sample_3(3,mc_i) = est_error(3, randi([1,tot_meas(2)],1));

%%%%%%% plot official %%%%%
% 
% figure;
% scatter(X(1,:),X(2,:),'b', '.');hold on;
% scatter(x_hat(1,:), x_hat(2,:),'g', '.');
% scatter( beacons_positions(1,:) ,beacons_positions(2,:), 60, 'r','filled');
% scatter(X(1,1),X(2,1),60,'r');
% scatter(X(1,k),X(2,k),60,'b', 'filled');
% 
% for beacon = 1:n_beacons 
%     viscircles([beacons_positions(1, beacon) beacons_positions(2, beacon)],beacon_range,'Color',[0.3010 0.7450 0.9330],'LineWidth',0.3 );
% end
% grid on; grid minor;
% xlabel('x[m]');
% ylabel('y[m]');
% legend('real position','estimate', 'beacons', 'Start', 'Robot')
% 
% 
% figure;
% grid on; grid minor;
% plot(t((N/8):(N-1)), est_error(1,(N/8):(N-1)),'b'); hold on;
% plot(t((N/8):(N-1)), est_error(2,(N/8):(N-1)),'r');
% plot(t((N/8):(N-1)), est_error(3,(N/8):(N-1)),'g');
% xlabel('t[s]')
% ylabel('error');
% legend('estimation error on x', 'estimation error on y', 'estimation error on theta');
% 
% 
% % Plot the real vehicle trajectory and the beacons
% figure('name','Real trajectory of vehicle')
% scatter( beacons_positions(1,:) ,beacons_positions(2,:), 60, 'r','filled'); 
% hold on;
% for beacon = 1:n_beacons 
%     viscircles([beacons_positions(1, beacon) beacons_positions(2, beacon)],beacon_range,'Color',[0.3010 0.7450 0.9330],'LineWidth',0.3 );
% end
% plot(X(1,:), X(2,:), 'b'); hold on;
% grid on; grid minor;
% xlabel('x [m]');
% ylabel('y [m]');
% legend('Beacon','Real trajectory');
% 
% % Plot the range measurements of distances
% if n_beacons < 10 && truth_table(1)
%     figure('name','Range measurements');
%     for beacon = 1:n_beacons
%         subplot_number = n_beacons*100+10+beacon;
%         subplot(subplot_number), plot(t, zn(beacon,:),'m'); 
%         hold on;
%         subplot(subplot_number), plot(t, z(beacon,:),'b');
%         title(['\fontsize{14} Distance Measurements beacon ',num2str(beacon)]);
%         xlabel('Time [sec]');
%         ylabel('d [m]');
%         legend('Measured','Real');
%     end
% end
% 
% % Plot the range measurements of orientations
% if n_beacons < 10 && truth_table(2)
%     figure('name','Range measurements');
%     for beacon = 1:n_beacons
%         subplot_number = n_beacons*100+10+beacon;
%         subplot(subplot_number), plot(t, zn(beacon,:),'m'); 
%         hold on;
%         subplot(subplot_number), plot(t, z(beacon,:),'b');
%         title(['\fontsize{14} Orientation Measurements beacon ',num2str(beacon)]);
%         xlabel('Time [sec]');
%         ylabel('rad [rad]');
%         legend('Measured','Real');
%     end
% end
% 
% 
% 
% 
% figure('name', 'Comparison between real and estimated pose (EKF approach)');
% subplot(311), plot(t, X(1,:),'b', t, x_hat(1,:),'r--');
% title('\fontsize{14} x(k): real and estimate');
% legend('Real','EKF estimate');
% ylabel('x [m]');
% xlabel('Time [s]');
% subplot(312),plot(t, X(2,:),'b', t, x_hat(2,:),'r--');
% title('\fontsize{14} y(k): real and estimate');
% legend('Real','EKF estimate');
% ylabel('y [m]');
% xlabel('Time [s]');
% subplot(313),plot(t, X(3,:),'b', t, x_hat(3,:),'r--');
% title('\fontsize{14} theta(k): real and estimate');
% legend('Real','EKF estimate');
% ylabel('[rad]');
% xlabel('Time [s]');
% 
% 
% 
% figure;
% subplot(311),plot(P1,'LineWidth',5); hold on;
% title('x position estimation error covariance [m^2] trend');
% xlabel('iteration');
% ylabel('x estimation error covariance [m^2]');
% legend('P(1,1)')
%  
% subplot(312),plot(P2,'LineWidth',5); hold on;
% title('y position estimation error covariance [m^2] trend');
% xlabel('iteration');
% ylabel('y estimation error covariance [m^2]');
% legend('P(2,2)')
% 
% subplot(313), plot(P3,'LineWidth',5); hold on;
% title('theta position estimation error covariance [rad^2] trend');
% xlabel('iteration');
% ylabel('theta error covariance [rad^2]');
% legend('P(3,3)')
end

    average_1 = sum(est_err_sample_1(1,:))/tot_meas(2);
    average_2 = sum(est_err_sample_2(2,:))/tot_meas(2);
    average_3 = sum(est_err_sample_3(3,:))/tot_meas(2);
