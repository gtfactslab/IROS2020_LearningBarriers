%============================= main ==================================
%
% @brief    Test file for an omnidirectional robot
%           Barrier computation is performed using an SVM kernel
%           Features: 1. Plots of true trajectory and SVM trajectory
%                     2. Comparison of point-wise control inputs from both
%                        controllers (RMSE error and a true signal values)
%                     3. Mismatch in the barrier function estimation
%============================= main ==================================
%
% @file     main.m
% @ingroup  Testing
% @author   Mohit Srinivasan, Amogh J. Dabholkar
%
% @date     2019/12/01
%
%!  Notes:
%!    set tabstop = 4, indent = 2. replace tab with spaces.
%
% @quitf
%============================= main ==================================
%==[0] Prep environment and working parms.
%
clear all;
close all;
clc;

bar = barrier.svm_lib;                          % Create SVM object

%==[1]  Implementation of barrier.
%

spec(:,1) = [-0.46,-0.47,0.2741,0.3741]; % Specs of first obstacle
spec(:,2) = [0.166,0.50,0.2741,0.2741];     % Specs of second obstacle
spec(:,3) = [0.415,-0.64,0.15,0.15];     % Specs of second obstacle
spec(:,4) = [ -1.025,0.49,0.2741,0.3741];     % Specs of second obstacle
spec(:,5) = [1.135,-0.06,0.2741,0.3741];     % Specs of second obstacle

dg = 0.01;                                  % 1cm grid size,
                                            % presuming units in meters

[bwd, grid] = barrier.interpolant.buildRobotarium(dg); % Build the interpolant object
bwd.specifyByRadii(spec);                   % Generate the interpolants.
% I = imread('Squiggle_2_trial.png');
% % h1 = axes;
% % image('XData', [-1.6, 1.6], 'YData', [-1, 1], 'CData', I);
% % set(h1, 'YDir', 'reverse');
% % axis([-1.6 1.6 -1 1])
% bW = imbinarize(I);
% bW = im2single(bW);
% bW = bW(:,:);
% bwd.squiggle(bW);

%==[2]  ROS Environment Setup
%
rosshutdown;
rosinit('http://localhost:11311')   % Initialize ROS master and global node.
sub = rossubscriber('/robot0/odom'); % Create a subscriber.
sub1 = rossubscriber('/robot0/laser_0'); % Subsribe to Lidar Data
pub = rospublisher('/robot0/cmd_vel','geometry_msgs/Twist');  % Create a publisher

[x_state, ~, ~, ~, ~] = bar.ROSinterface(sub,sub1);
msg = rosmessage(pub); % Empty message determined by the topic published by pub.
msg.Linear.X = 0.0;    % Initialize Robot with 0 velocities.
msg.Linear.Y = 0.0;
msg.Linear.Z = 0.0;
msg.Angular.X = 0.0;
msg.Angular.Y = 0.0;
msg.Angular.Z = 0.0;
send(pub,msg);
% ==[3]  Plot on MATLAB.

% PlotGoalsObstacles();
% 
% Plt_data1 = [];        % Plot the path of Robot
% Plt_data1 = [Plt_data1; x_state(1,1); x_state(2,1)];
% p1 = plot(Plt_data1(1), Plt_data1(2), 'k-.', 'LineWidth', 3);
% drawnow

X_Obsneg = [];
Y_Obsneg = [];
X_Obspos = [];
Y_Obspos = [];
X_neg = [];
Y_neg = [];
X_pos = [];
Y_pos = [];
gdatasdist = [];
x_goal = [2.7 - 1.6; 1.57 - 1];
% % x_goal = [2.7 - 1.6; 1.9 - 1];
% plot(x_goal(1), x_goal(2),'rX')             % Plot goal location
% plot(x_state(1), x_state(2), 'gX')          % PLot initial condition
% axis([-1.6 1.6 -1 1])
% hold on    

h_eval = [];
iter = [];
i = 0;
controller =  1;                            % Choose type of controller:
                                            % 1: SVM Kernel controller
                                            % (black)
                                            % 2: Bwdist controller (red)

accum = 2;                                  % Choose type of sampling:
                                            % 1: Keep track of past samples
                                            % 2: Keep local samples

cScale = 0.2;
 

while(norm(x_state(1:2,:) - x_goal) >= 0.1)     % Solve until goal is reached
    
    [x_state, x_rot, y_rot, x_rotplus, y_rotplus] = bar.ROSinterface(sub,sub1);

%     Plt_data1 = [Plt_data1, [x_state(1,1); x_state(2,1)]];  % Plot the path of Robot
%     p1.XData = Plt_data1(1,:);
%     p1.YData = Plt_data1(2,:);
%     hold on
   
%     X_neg = [X_neg, x_rot];               % All negative samples
%     Y_neg = [Y_neg, y_rot];               % from LiDAR data
% 
%     X_pos = [X_pos, x_rotplus];           % All positive samples
%     Y_pos = [Y_pos, y_rotplus];           % projected from negatives

    if accum == 1
        
        X_Obsneg = [X_Obsneg, x_rot];               % All negative samples
        Y_Obsneg = [Y_Obsneg, y_rot];               % from LiDAR data

        X_Obspos = [X_Obspos, x_rotplus];           % All positive samples
        Y_Obspos = [Y_Obspos, y_rotplus];           % projected from negatives
        
    else                                                
        X_Obsneg = x_rot;                           % All negative samples
        Y_Obsneg = y_rot;                           % from LiDAR data

        X_Obspos = x_rotplus;                       % All positive samples
        Y_Obspos = y_rotplus;                       % projected from negatives
    end
    
%     samp = [X_neg,X_pos;Y_neg,Y_pos];   % Training data
    samp = [X_Obsneg,X_Obspos;Y_Obsneg,Y_Obspos];   % Training data
    classif = [ones(1,length(X_Obsneg)),zeros(1,length(X_Obspos))];
    u_des = -0.1*(x_goal- x_state(1:2,:));      % Nominal control
    if isempty(samp)
        
        u = u_des;
           
    elseif controller == 1
        
        [gridFine,grids,bord_samp,shapeFine,shapeGrid] = bar.buildSTDR(x_state); 
        dist = bar.evalSVM(grids, gridFine, shapeFine, bord_samp, samp, classif);
        
%       Plot zero-level set of SVM-kernel estimated barrier function
%         contour(bar.xx, bar.yy, dist, [0 0])
%         axis equal;
%         xlim([-bar.grid_x, bar.grid_x]);
%         ylim([-bar.grid_y, bar.grid_y]);

        [h, A, B, u] = bar.barrier_constraint(x_state(1:2,:), u_des);
        display(h);
        
        if norm(u) <= 0.01                  % Check for deadlock
            accum = 2;                      
        else
            accum = 1;
        end
        
        i = i+1;
        iter = [iter, i];

        [h_bwdist, u_true] = barrier_bwdist(x_state(1:2,:), bwd, u_des);   % Obs is barrier measure. Barr is distance from goal measure. u is the control input.
        h_eval = [h_eval, h_bwdist];
    else    
        
        [h_bwdist, u] = barrier_bwdist(x_state(1:2,:), bwd, u_des);   % Obs is barrier measure. Barr is distance from goal measure. u is the control input.
        i = i+1;
        iter = [iter, i];

    end
    
    u = cScale*u;
    bar.generateVelCmdsOmni(u, pub);
    
end

%==[6]  ROS Shutdown
bar.shutdownROS(pub);
if controller == 2
%     save('bwdist_5obs_ic1.mat','Plt_data1');
    save('BarrierEvaluation', 'h_eval')
else
%     save('squiggle_online_ic1.mat','Plt_data1');
    save('BarrierEvaluation', 'h_eval')
end

% save('samp_2.mat', 'samp', 'X_pos', 'X_neg', 'Y_pos', 'Y_neg');
% save('samp_squiggle_2.mat', 'samp', 'X_pos', 'X_neg', 'Y_pos', 'Y_neg');

% contour(bar.xx, bar.yy, dist, [0 0])
% axis equal;
% xlim([-bar.grid_x, bar.grid_x]);
% ylim([-bar.grid_y, bar.grid_y]);
% hold on

%==[7]  Plotting Function Definitions
%

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

%% TODO: Modify this to reflect new workspace dimensions
function plot_Lp(nObs)

    obst.power = 2*4;            % p = 8 in L_p part. Always even!

    switch (nObs)
          case 1
              obst.pose(1:3, 1) = [ 0.8; 1.3; -pi/4 ];
              obst.dims(1:2, 1) = [ 1; 0.3 ]/2;

          case 2
              obst.pose(1:3, 1) = [ 0.8; 1.3; -pi/4 ]; % Default Number of Obstacles
              obst.dims(1:2, 1) = [ 1; 0.3 ]/2;

              obst.pose(1:3, 2) = [ 1.6; 0.3; pi/20 ];
              obst.dims(1:2, 2) = [ 0.5; 0.3 ]/2;

          case 3
              obst.pose(1:3, 1) = [ 0.8; 1.3; -pi/4 ];
              obst.dims(1:2, 1) = [ 1; 0.3 ]/2;

              obst.pose(1:3, 2) = [ 1.6; 0.3; pi/20 ];
              obst.dims(1:2, 2) = [ 0.5; 0.3 ]/2;

              obst.pose(1:3, 3) = [ 2.2; 1.4; pi/2 ];
              obst.dims(1:2, 3) = [ 0.4; 0.4 ]/2;


          case 4
              obst.pose(1:3, 1) = [ 0.8; 1.3; -pi/4 ];
              obst.dims(1:2, 1) = [ 1; 0.3 ]/2;

              obst.pose(1:3, 2) = [ 1.6; 0.3; pi/20 ];
              obst.dims(1:2, 2) = [ 0.5; 0.3 ]/2;

              obst.pose(1:3, 3) = [ 2.2; 1.4; pi/2 ];
              obst.dims(1:2, 3) = [ 0.4; 0.4 ]/2;

              obst.pose(1:3, 4) = [ 2.6; 0.6; pi/3 ];
              obst.dims(1:2, 4) = [ 0.5; 0.2 ]/2;


          case 5
              obst.pose(1:3, 1) = [ 0.8; 1.3; -pi/4 ];
              obst.dims(1:2, 1) = [ 1; 0.3 ]/2;

              obst.pose(1:3, 2) = [ 1.6; 0.3; pi/20 ];
              obst.dims(1:2, 2) = [ 0.5; 0.3 ]/2;

              obst.pose(1:3, 3) = [ 2.2; 1.4; pi/2 ];
              obst.dims(1:2, 3) = [ 0.4; 0.4 ]/2;

              obst.pose(1:3, 4) = [ 2.6; 0.6; pi/3 ];
              obst.dims(1:2, 4) = [ 0.5; 0.2 ]/2;

              obst.pose(1:3, 5) = [ 0.6; 0.5; pi/10 ];
              obst.dims(1:2, 5) = [ 0.4; 0.6 ]/2;
    end

    R = zeros(2, 2, nObs);          % Will hold rotation to apply.

    obst_Lp_box = {};

    for ii = 1:nObs

      % Compute corners for each obstacle: first in obstacle (body) frame,
      % then in world frame.
      obsCorners_b = [obst.dims(1, ii), obst.dims(1, ii), ...
          -obst.dims(1, ii), -obst.dims(1, ii), obst.dims(1, ii) ; ...
          obst.dims(2, ii), -obst.dims(2, ii), ...
          -obst.dims(2, ii), obst.dims(2, ii), obst.dims(2, ii)];

      R(:, :, ii) = [cos(obst.pose(3, ii)), -sin(obst.pose(3, ii)) ; ...
          sin(obst.pose(3, ii)),  cos(obst.pose(3, ii)) ];

      obsCorners_s = R(:, :, ii)*obsCorners_b;

      % Plot Lp square obstacle by iterating around the boundary according to
      % rays at specified angles (phi). Use phi and obstacle length to solve
      % for the boundary radius distance along that ray. Gives boundary in the
      % obstacle frame.
      phi = 0:pi/200:2*pi;
      rad_search = 0:0.1:100;
      Lp_box = zeros(2, length(phi));

      epNorm = (cos(phi)/obst.dims(1, ii)).^(obst.power) ...
          + (sin(phi)/obst.dims(2, ii)).^(obst.power);
      rho    = nthroot(1./epNorm,(obst.power));

      % Transform boundary points to world frame.
      Lp_box = R(:,:,ii)*[rho.*cos(phi) ; rho.*sin(phi)];
      obst_Lp_box{ii} = bsxfun(@plus, Lp_box, obst.pose([1:2],ii));

      plot(obst_Lp_box{ii}(1, :), obst_Lp_box{ii}(2, :), 'k-','LineWidth',3); % Plot Lp obstacle
      hold on;
    end
    xlim([-1.6, 1.6]); % Modified limits on X Axis
    ylim([-1 1]); % Modified limits on Y Axis

end
