classdef UnicycleBot  < handle
    %UnicycleBot Create a simulated unycicle-like robot object
    
    properties(SetAccess = public)
        %Dt Time step for simulation of the robot
        Dt
        
        %Trace Pose history of the robot
        Trace = []
        
        %BGs Best estimated (guessed) robot trace
        BGs = []
        
        %Pose Current pose of the robot
        Pose
        
        %L wheelbase
        L = 0.2
    end
    
    properties(Access = private)
        %HCenter Graphics handle for center of rear axle of the car robot
        HCenter
        
        %HTrajectory Graphics handle for the trajectory of the car
        HTrajectory
        
        %HChassis Graphics handle for chassis
        HChassis
        
        %HAngleCov Graphics handle for uncertanty in heading direction
        HAngleCov
        
        %HEllipse Graphics handle for uncertanty elipse in xy plance
        HPosCov
        
        %HParticles Graphics handle for particles
        %HParticles
        
        %HRw1 Graphics handle for rear wheel 1
        HRw1
        
        %HRw2 Graphics handle for rear wheel 2
        HRw2
        
        %HBestGuesses Graphics handle for best current pose estimation
        HBestGuesses
        
        %Cnt1 Internal counter 1
        Cnt1
        
        %FigureHandle the handle of the figure
        FigureHandle
        
        %AxesHandle the handle of the plot axes
        AxesHandle
    end
    
    
    methods
        
        
        function obj = UnicycleBot(currentPose, dt, varargin)
            %ExampleHelperCarBot Constructor
            
            
            obj.Dt = dt;
            obj.Pose = currentPose;
            obj.Trace = obj.Pose;
            obj.Cnt1 = 0;
            
            if ~isempty(varargin)
                obj.L = varargin{1};
                obj.FigureHandle = figure('Name', 'CarBot');
                % clear the figure
                obj.AxesHandle = axes(obj.FigureHandle);
                cla(obj.AxesHandle)

                % customize the figure
                obj.FigureHandle.Position = [100 100 1000 500];
                axis(obj.AxesHandle, 'equal');
                xlim(obj.AxesHandle, [-0.1,12]);
                ylim(obj.AxesHandle, [-1,4]);
                grid(obj.AxesHandle, 'on');
                box(obj.AxesHandle, 'on');

                hold(obj.AxesHandle, 'on')

                %obj.HParticles = scatter(ax, 0,0,'MarkerEdgeColor','g', 'Marker', '.');
                obj.HPosCov     = plot([1],[1]);
                obj.HAngleCov   = fill([1,1,1],[1,1,1],'y');
                obj.HCenter     = plot(obj.AxesHandle, obj.Pose(1), obj.Pose(2),'o'); % center of the rear axle of the car robot
                obj.HTrajectory = plot(obj.AxesHandle, obj.Trace(1,:), obj.Trace(2,:),'r'); % car trace (of best guess)

                obj.HChassis = plot(obj.AxesHandle, 0,0,'b'); % chassis
                obj.HRw1 = plot(obj.AxesHandle, 0,0,'linewidth',3, 'color','b');         % rear wheel 1
                obj.HRw2 = plot(obj.AxesHandle, 0,0,'linewidth',3, 'color','b');         % rear wheel 2


                obj.HBestGuesses = plot(obj.AxesHandle, 0,0,'ks-'); % best guess of pose

                legend(obj.AxesHandle, [obj.HTrajectory, obj.HBestGuesses],'actual pose','estimated pose','Location','northwest');

                title(obj.AxesHandle, 't = 0');
                xlabel(obj.AxesHandle, 'x (m)');
                ylabel(obj.AxesHandle, 'y (m)');
                hold(obj.AxesHandle, 'off');

                obj.drawBot()
            end
        end
        
        function [angle] = normalize_rad(~,angle)
            %normalize_rad normalize angle in radiant between -pi and pi
            while angle>pi
                angle = angle - 2*pi;
            end
            while angle<-pi
                angle = angle + 2*pi;
            end
        end
        
        function rho = path_folowing(obj,P_frenet)
            %Comput the curvature to follow to reach the frenet point
            
            k = 10; % control constant
            
            % Estract useful values
            dP_frenet(1:3,1) = P_frenet(1:3) - obj.Pose; % delta pose poin_pose-bot_pose
            theta_f =  obj.normalize_rad(P_frenet(3));    % absolute angle of frenet point
            thetaDes =  obj.normalize_rad(dP_frenet(3));
            l = -[cos(theta_f+pi/2),sin(theta_f+pi/2)]*dP_frenet(1:2);
            c = P_frenet(4); % curvature in the point
            
            thetaTilde = -thetaDes;
            delta = -pi/2*tanh(l/2);
            
            Ddelta_dl = -pi/2*(1-tanh(l/2)^2)/2;
            gamma = c*cos(thetaTilde)/(1-c*l) + sin(thetaTilde)*Ddelta_dl;
            
            tmp = obj.normalize_rad(thetaTilde-delta);
            rho = gamma - k*tmp;
        end
        
        
        function drive(obj, uCmd)
            %Drive Move the robot forward
            % contaminate commanded motion with noise
            u = uCmd + 0.001*randn(1,2)*0;
            
            poseDot = 0*obj.Pose;
            poseDot(1) = u(1)*cos(obj.Pose(3));
            poseDot(2) = u(1)*sin(obj.Pose(3));
            poseDot(3) = u(2);
            
            obj.Pose = obj.Pose + obj.Dt*poseDot;
            
            obj.Trace = [obj.Trace, obj.Pose];
            
        end
        
        function [h,h_dot] = rfidReadings(obj,uCmd,RFID)
            %rfidReadings simulate the reading of RFIDs
            %   Detailed explanation goes here
            v_norm = uCmd(1);
            omega  = uCmd(2);
            
            measure_num = size(RFID,2);
            
            % SLOW SOLUTION
            h = zeros(measure_num,1);
            h_dot = zeros(measure_num,1);
             
            x_p = obj.Pose(1,:);
            y_p = obj.Pose(2,:);
            theta_p = obj.Pose(3,:);
            % v_norm = 0.8; %state(4,:);
            
%             v = v_norm.*[cos(theta_p); sin(theta_p)];
            
            for i=1:measure_num
                x_rfid = RFID(1,i);
                y_rfid = RFID(2,i);
                
                dx = x_rfid - x_p;
                dy = y_rfid - y_p;
                
                alpha = atan2(dy,dx) - theta_p;
                
                normaliz = sqrt(dx.*dx + dy.*dy);
                normaliz(normaliz<1e-4) = Inf;
                
                % %         % SLOW VERSION
                % %         versore = [dx;dy]/normaliz;
                % %         d_dot_norm = v'*versore;
                % %         h(i) = d_dot_norm ;
                
                % FAST VERSION
                % h(i,:) = (v(1,:).*dx + v(2,:).*dy)./normaliz;                
                h(i,:) = v_norm*cos(alpha);
%                 if norm(h(i,:)-(v(1,:).*dx + v(2,:).*dy)./normaliz) > 1e-5
%                    error('funzione implementata male'); 
%                 end
                h_dot(i,:) = v_norm*cos(alpha)/normaliz - sin(alpha)*v_norm*omega;
            end
            
            %if(numel(nargout))>1
            %varargout{1} = h_dot;
            %end
        end
        
        
        function drawBot(obj)
            %drawCarBot Routine to draw the car-like robot
            r = obj.L/6;
            p = obj.Pose;
            x = p(1);
            y = p(2);
            theta = p(3);
            
            
            chassis = [x  + obj.L/2*cos(theta-pi/2), y + obj.L/2*sin(theta-pi/2);
                x + obj.L/2*cos(theta+pi/2), y + obj.L/2*sin(theta+pi/2);
                x, y;
                x + obj.L*cos(theta), y + obj.L*sin(theta);
                ];
            
            rw1 = [ x  + obj.L/2*cos(theta-pi/2) + r*cos(theta), y  + obj.L/2*sin(theta-pi/2) + r*sin(theta);
                x  + obj.L/2*cos(theta-pi/2) - r*cos(theta), y  + obj.L/2*sin(theta-pi/2) - r*sin(theta)];
            rw2 = [ x  + obj.L/2*cos(theta+pi/2) + r*cos(theta), y  + obj.L/2*sin(theta+pi/2) + r*sin(theta);
                x  + obj.L/2*cos(theta+pi/2) - r*cos(theta), y  + obj.L/2*sin(theta+pi/2) - r*sin(theta)];
            
            
            obj.HChassis.XData = chassis(:,1);
            obj.HChassis.YData = chassis(:,2);  % chassis
            
            obj.HRw1.XData = rw1(:,1);
            obj.HRw1.YData = rw1(:,2);  %rear wheel 1
            
            obj.HRw2.XData = rw2(:,1);
            obj.HRw2.YData = rw2(:,2);  %rear wheel 2
        end
        
        function updatePlot(obj, currentBestGuess, currentBestGuess_cov, t)
            % updatePlot
            
            obj.Cnt1 = obj.Cnt1 + 1;
            
            % render particles
            %particles = particleFilter.Particles;
            %obj.HParticles.XData = particles(1:end,1);
            %obj.HParticles.YData = particles(1:end,2);
            
            if obj.Cnt1 == 8
                obj.BGs = [obj.BGs; currentBestGuess(:)'];
                obj.Cnt1 = 0;
                % draw best estimated robot trace (by observer)
                obj.HBestGuesses.XData = obj.BGs(:,1);
                obj.HBestGuesses.YData = obj.BGs(:,2);
            end
            
            
            % draw car rear axle center
            obj.HCenter.XData = obj.Pose(1);
            obj.HCenter.YData = obj.Pose(2);
            
            % Plot xy covariance
            [x_el,y_el] = error_ellipse(currentBestGuess_cov(1:2,1:2),currentBestGuess(1:2),0.95);
            obj.HPosCov.XData = x_el;
            obj.HPosCov.YData = y_el;
            
            % Plot angle covariance
            ang_cov = 3*sqrt(currentBestGuess_cov(3,3));
            n = ceil(ang_cov/0.05);
            angles = currentBestGuess(3) + ang_cov*linspace(-1,1,2*n);
            obj.HAngleCov.XData = currentBestGuess(1)+[0,cos(angles),0]*obj.L/3;
            obj.HAngleCov.YData = currentBestGuess(2)+[0,sin(angles),0]*obj.L/3;
            
            % draw car trajectory
            obj.HTrajectory.XData = obj.Trace(1,:);
            obj.HTrajectory.YData = obj.Trace(2,:);
            
            % draw car-like robot
            obj.drawBot();
            
            ax = get(obj.FigureHandle, 'currentaxes');
            title(ax, ['t = ', num2str(t)]);
            
            % capture snapshots for time frame 101 and 265
            if floor(t/obj.Dt) == 101
                snapnow
            end
            
            if floor(t/obj.Dt) == 265
                snapnow
            end
        end        
    end
end

