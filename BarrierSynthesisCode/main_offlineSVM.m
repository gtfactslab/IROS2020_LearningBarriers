close all;
clc;

rosshutdown;
rosinit('http://localhost:11311')   % Initialize ROS master and global node.
sub = rossubscriber('/robot0/odom'); % Create a subscriber.
sub1 = rossubscriber('/robot0/laser_0'); % Subsribe to Lidar Data
pub = rospublisher('/robot0/cmd_vel','geometry_msgs/Twist');  % Create a publisher

bar = barrier.svm;                          % Create SVM object

spec(:,1) = [-0.46,-0.47,0.2741,0.3741]; % Specs of first obstacle
spec(:,2) = [0.166,0.50,0.2741,0.2741];     % Specs of second obstacle
spec(:,3) = [0.415,-0.64,0.15,0.15];     % Specs of second obstacle
spec(:,4) = [ -1.025,0.49,0.2741,0.3741];     % Specs of second obstacle
spec(:,5) = [1.135,-0.06,0.2741,0.3741];     % Specs of second obstacle

dg = 0.01;                                  % 1cm grid size,
                                            % presuming units in meters
                                            
[bwd, grid] = barrier.interpolant.buildRobotarium(dg); % Build the interpolant object
bwd.specifyByRadii(spec);                   % Generate the interpolants.

[x_state, ~, ~, ~, ~] = bar.ROSinterface(sub,sub1);
msg = rosmessage(pub); % Empty message determined by the topic published by pub.
msg.Linear.X = 0.0;    % Initialize Robot with 0 velocities.
msg.Linear.Y = 0.0;
msg.Linear.Z = 0.0;
msg.Angular.X = 0.0;
msg.Angular.Y = 0.0;
msg.Angular.Z = 0.0;
send(pub,msg);

% x_goal = [2.7 - 1.6; 1.57 - 1];
% plot(x_goal(1), x_goal(2),'rX')                    % Plot goal location
% hold on
% plot(x_state(1), x_state(2), 'gX')                 % PLot initial condition
% axis([-1.6 1.6 -1 1])
% hold on    
% PlotGoalsObstacles();
% hold on
iter = [];
h_eval = [];
i = 0;

% Plt_data1 = [];        % Plot the path of Robot
% Plt_data1 = [Plt_data1; x_state(1,1); x_state(2,1)];
% p1 = plot(Plt_data1(1), Plt_data1(2), 'k-.', 'LineWidth', 3);
% drawnow

xx = load('Offline5Obs_BiasedSVM_2.mat', 'xx');
xx = xx.xx;
yy = load('Offline5Obs_BiasedSVM_2.mat', 'yy');
yy = yy.yy;

gdata_sdist = load('Offline5Obs_BiasedSVM_2.mat', 'gdata_sdist');
h = gdata_sdist.gdata_sdist;
[px, py] = gradient(h);                                        % Gradient of barrier function
h = imgaussfilt(h, 2);
h = h*0.01;
h = -h;
% figure(5);
% surf(xx, yy, gdata_sdist.gdata_sdist);
% hold off;
F = griddedInterpolant(xx', yy', h', 'cubic'); 
dfx = griddedInterpolant(xx',yy',px');                                    % H' wrt x
dfy = griddedInterpolant(xx',yy',py');                                    % H' wrt y

% contour(xx,yy, h, [0 0])
% axis equal;
% xlim([-1.6, 1.6]);
% ylim([-1, 1]);
% hold on

while(norm(x_state(1:2,:) - x_goal) >= 0.1)     % Solve until goal is reached
    
    [x_state, ~, ~, ~, ~] = bar.ROSinterface(sub,sub1);

%     Plt_data1 = [Plt_data1, [x_state(1,1); x_state(2,1)]];  % Plot the path of Robot
%     p1.XData = Plt_data1(1,:);
%     p1.YData = Plt_data1(2,:);
%     hold on

    
    u_des = -0.1*(x_goal- x_state(1:2,:));      % Nominal control
    [h, u] = nav_domain(x_state(1:2,:), dfx,dfy, F, u_des);
    i = i+1;
    iter = [iter, i];

    [h_bwdist, u_true] = barrier_bwdist(x_state(1:2,:), bwd, u_des);   % Obs is barrier measure. Barr is distance from goal measure. u is the control input.
    h_eval = [h_eval, h_bwdist];
    u = 2*u;
    bar.generateVelCmdsOmni(u, pub);

end
bar.shutdownROS(pub);
% save('squiggle_offline_ic1','Plt_data1')
save('BarrierEvaluation', 'h_eval')

%% Plot all the ellipsoidal/circular goals and obstacles
function PlotGoalsObstacles()

    spec(:,1) = [-0.46,-0.47,0.2741,0.3741]; % Specs of first obstacle
    spec(:,2) = [0.166,0.50,0.2741,0.2741];     % Specs of second obstacle
    spec(:,3) = [0.415,-0.64,0.15,0.15];     % Specs of second obstacle
    spec(:,4) = [ -1.025,0.49,0.2741,0.3741];     % Specs of second obstacle
    spec(:,5) = [1.135,-0.06,0.2741,0.3741];     % Specs of second obstacle

    P1 = [1/spec(3,1)^2 0; 0 1/(0.3741)^2];
    P2 = [1/spec(3,2)^2 0; 0 1/spec(3,2)^2];

    P3 = [1/spec(3,3)^2 0; 0 1/spec(3,3)^2];
    P4 = [1/spec(3,4)^2 0; 0 1/(0.3741)^2];
    P5 = [1/spec(3,5)^2 0; 0 1/(0.3741)^2];

    c1 = spec(1,1); c2 = spec(2,1); % Centers of obstacles
    c3 = spec(1,2); c4 = spec(2,2);

    c5 = spec(1,3); c6 = spec(2,3); 
    c7 = spec(1,4); c8 = spec(2,4);
    c9 = spec(1,5); c10 = spec(2,5);

    plot_ellipse(P1, c1, c2, 'k');
    hold on
    plot_ellipse(P2, c3, c4, 'k');
    hold on
    plot_ellipse(P3, c5, c6, 'k');
    hold on
    plot_ellipse(P4, c7, c8, 'k');
    hold on
    plot_ellipse(P5, c9, c10, 'k');
    hold on
    axis equal;
    xlim([-1.6, 1.6]);
    ylim([-1, 1]);
    set(gcf, 'color', 'w')
    
end

%% Plot ellipsoidal/circular level-sets
function plot_ellipse(P, a, b, c)

    theta = 0:0.00001:2*pi;
    x = (1/sqrt(P(1,1)))*cos(theta) + a;
    y = (1/sqrt(P(2,2)))*sin(theta) + b;
    plot(x,y,c,'LineWidth',3);
    hold on

end
