%========================== barrier.svm ==========================
%
% @class    barrier.svm
%
% @brief    A class which uses certain default functions of MATLAB for
%           computation of barrier function and its gradient point-wise
%           in time.
%
%           Included in this class is also a member function which
%           generates the control input from a Quadratic Program (QP) from
%           the SVM-Kernel estimated barrier function.
%
%========================== barrier.svm ==========================
%
% @file     svm.m
%
% @author   Amogh J. Dabholkar, adabholkar6@gatech.edu
% @date     2019/12/01
%
% @note
%   set indent to 2 spaces.
%   set tab to 4 spaces (with conversion to spaces).
%
% @classf   barrier
%========================== barrier.interpolant ==========================
classdef svm_lib < barrier.base

properties (SetAccess = protected, GetAccess = public)
  F;         %< H(x)
  dfx;       %< H' wrt x
  dfy;       %< H' wrt y
  xx;        %< X Gridding
  yy;        %< Y Gridding
  grid_x;    %< X domain length
  grid_y;    %< Y domain length
  dg;        %< Grid Spacing
  sigma;     %< kernel sigma
  fineSigma  %< Fine kernel sigma
  kernSVM;   %< SVM Parameters
  kernArgs;  %< SVM Parameters
  nObs;      %< Number of obstacles
  extras;    %< Build Parameters
end

methods
    
%================================ svm ===============================
% 
% @brief  Constructor for the SVM-Barrier class.
%
function this = svm_lib()

  this.grid_x = 1.6; 
  this.grid_y = 1;
  this.dg = 0.01;
  this.sigma = (75*0.01)^2;                           % Bandwidth
  this.fineSigma = (25*0.01)^2;                       % Finer bandwidth

  %--[1.2] SVM params
  this.kernSVM  = 'rbf';      % SVM type (i.e. boundary)
%  this.kernArgs = {'PolynomialOrder',2};

  %--[1.3] Build example.
  this.nObs = 2;
  this.extras.grid = true;
  this.extras.gridType = 'inside';
  this.extras.gridFine = true;
  this.extras.boundary = true;
  this.extras.bRad     = sqrt(this.sigma);

end

%============================= generateFromArray =============================
%
% @brief  Define the interpolant.
%
function generateFromArray(this, gdatasdist)

if all(size(gdatasdist) == size(this.xx))
  gdatasdist = -0.1*(gdatasdist - 0.35);
  this.F = griddedInterpolant(this.xx',this.yy',gdatasdist','cubic');   % H(x)

  [px,py] = gradient(gdatasdist);

  this.dfx = griddedInterpolant(this.xx',this.yy', px', 'cubic');          % H' wrt x
  this.dfy = griddedInterpolant(this.xx',this.yy', py', 'cubic');          % H' wrt y

else
  error('Array dimensions do not align with the expected size')
end
end

%========================= specifyByLpConstraint ==========================
%
% @brief  Define the obstacles as Lp boxes.
%
function specifyByLpConstraint(this)
%TODO : Remove hard-coding.
if (~exist('this.nObs'))
  prompt = 'Enter number of obstacles: ';
  this.nObs = input(prompt);
  fprintf('Using nObs = %d.\n', this.nObs);
else
  fprintf('Using nObs = %d.\n', this.nObs);
end

obst.power = 2*4;            % p = 8 in L_p part. Always even!

switch (this.nObs)
  case 1
      obst.pose(1:3, 1) = [ -0.8; 0.3; -pi/4 ];
      obst.dims(1:2, 1) = [ 1; 0.3 ]/2;

  case 2
      obst.pose(1:3, 1) = [ -0.8; 0.3; -pi/4 ]; % Default Number of Obstacles
      obst.dims(1:2, 1) = [ 1; 0.3 ]/2;

      obst.pose(1:3, 2) = [ 0; -0.4; pi/20 ];
      obst.dims(1:2, 2) = [ 0.5; 0.3 ]/2;

  case 3
      obst.pose(1:3, 1) = [ -0.8; 0.3; -pi/4 ];
      obst.dims(1:2, 1) = [ 1; 0.3 ]/2;

      obst.pose(1:3, 2) = [ 0; -0.4; pi/20 ];
      obst.dims(1:2, 2) = [ 0.5; 0.3 ]/2;

      obst.pose(1:3, 3) = [ 0.8; 0.4; pi/2 ];
      obst.dims(1:2, 3) = [ 0.4; 0.4 ]/2;


  case 4
      obst.pose(1:3, 1) = [ -0.8; 0.3; -pi/4 ];
      obst.dims(1:2, 1) = [ 1; 0.3 ]/2;

      obst.pose(1:3, 2) = [ 0; -0.4; pi/20 ];
      obst.dims(1:2, 2) = [ 0.5; 0.3 ]/2;

      obst.pose(1:3, 3) = [ 0.8; 0.4; pi/2 ];
      obst.dims(1:2, 3) = [ 0.4; 0.4 ]/2;

      obst.pose(1:3, 4) = [ 1; -0.2; pi/3 ];
      obst.dims(1:2, 4) = [ 0.5; 0.2 ]/2;


  case 5
      obst.pose(1:3, 1) = [ -0.8; 0.3; -pi/4 ];
      obst.dims(1:2, 1) = [ 1; 0.3 ]/2;

      obst.pose(1:3, 2) = [ 0; -0.4; pi/20 ];
      obst.dims(1:2, 2) = [ 0.5; 0.3 ]/2;

      obst.pose(1:3, 3) = [ 0.8; 0.4; pi/2 ];
      obst.dims(1:2, 3) = [ 0.4; 0.4 ]/2;

      obst.pose(1:3, 4) = [ 1; -0.2; pi/3 ];
      obst.dims(1:2, 4) = [ 0.5; 0.2 ]/2;

      obst.pose(1:3, 5) = [ -1; -0.5; pi/10 ];
      obst.dims(1:2, 5) = [ 0.4; 0.6 ]/2;
end

R = zeros(2, 2, this.nObs);          % Will hold rotation to apply.

obst_Lp_box = {};

for ii = 1:this.nObs

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

%       plot(obst_Lp_box{ii}(1, :), obst_Lp_box{ii}(2, :), 'k-'); % Plot Lp obstacle
end

% Define function that plots the obstacles.
plotObs = @(ob, cspec, i)eval(['for i=1:length(ob),', ...
                        'plot(ob{i}(1,:),ob{i}(2,:),cspec), end']);

end
%============================ evalBoundary ============================
%
% @brief  Evaluates the boundary from the data points.
%
function [cError,evalGrid,theMap] = evalBoundary(this,gridFine,grid,bord_samp,samp,classif,shapeFine,shapeGrid)

evalGrid  = kernelFunction.grid2centers(gridFine);

theKernel = kernelGaussian(this.sigma);  % Instantiate kernel object.
if (isfield(this.extras,'boundary') && this.extras.boundary)
    [gridCent] = kernelFunction.grid2centers(grid);
    gridCent = [gridCent , bord_samp];
    sparseMap = kernMach.buildFromCenters(theKernel, gridCent);
    gCard = size(gridCent, 2);
else
    sparseMap = kernMach.buildFromGrid(theKernel, grid);
    gCard = prod(shapeGrid);
end
theMap = multiKernMach(sparseMap);
tic
theMap.estimate(double(classif), samp);
time1 = toc;
yOut = transpose(reshape(theMap.eval(evalGrid), shapeFine));
ySamp = theMap.eval(samp);
cSamp = (ySamp > 0.5);
cError = xor(cSamp, classif);
mapCenters = theMap.centers;

end
%============================ revalBoundary ============================
%
% @brief  Re-evaluates the boundary from the data points.
%
function revalBoundary(this,cError,samp,classif,shapeFine,evalGrid,theMap)

errCent = samp(:,cError);
errClass = classif(cError);

fineKern  = kernelGaussian(this.fineSigma);  % Instantiate kernel object.
errMap = kernMach(fineKern, errCent);
theMap.addMachine(errMap);


tic
theMap.estimateFrom(2, double(classif), samp);
time1 = toc;


yOut = transpose(reshape(theMap.eval(evalGrid), shapeFine));
ySamp = theMap.eval(samp);
cSamp = (ySamp > 0.5);
cError = xor(cSamp, classif);
mapCenters = theMap.centers;

end
%============================ evalSVM ============================
%
% @brief  Calculate the actual signed distance field using SVM.
%
function gdata_sdist = evalSVM(this,grid,gridFine,shapeFine,bord_samp,samp,classif)

theKernel = kernelGaussian(this.sigma);  % Instantiate kernel object.
if (isfield(this.extras,'boundary') && this.extras.boundary)
    [gridCent] = kernelFunction.grid2centers(grid);
    gridCent = [gridCent , bord_samp];
    sparseMap = kernMach.buildFromCenters(theKernel, gridCent);
    gCard = size(gridCent, 2);
else
    sparseMap = kernMach.buildFromGrid(theKernel, grid);
    gCard = prod(shapeGrid);
end
theMap = multiKernMach(sparseMap);

if (exist('bord_samp'))
    qBord = theMap.maptoH(bord_samp);
end

tic
qSamp = theMap.maptoH(samp);
time2 = toc;


ySVM = 2*(double(classif) - 0.5);
yTau = 0;
mapCenters = theMap.centers;
tic
% svmach = fitcsvm(qSamp',ySVM', 'Standardize',false, ...
%     'KernelFunction',this.kernSVM, 'BoxConstraint', 1000);
b = [0 1; 2 0];
svmach = fitcsvm(qSamp',ySVM', 'Standardize',false, ...
    'KernelFunction',this.kernSVM, 'Cost', b);
timeSVM = toc;
display(timeSVM);
[this.xx, this.yy] = meshgrid(gridFine{:});
gdata_kpca_coeff = theMap.maptoH( transpose([this.xx(:), this.yy(:)]) );
tic
[gdata_label, gdata_score] = predict(svmach, gdata_kpca_coeff');
timePredict = toc;
% display(timePredict);
gdata_posneg = reshape(gdata_label > yTau, flip(shapeFine));
cSamp = predict(svmach, qSamp') > yTau;
cError = xor(cSamp', classif);
fprintf('Missclassifications from SVM = %d.\n\n', nnz(cError));

gdata_dist  = reshape(max(gdata_score, [], 2), flip(shapeFine));
gdata_sdist = 2*(0.5-gdata_posneg).*gdata_dist;
this.F = imgaussfilt(gdata_sdist, 2);
this.F = this.F*this.dg;
% gdatasdist = [gdatasdist, gdata_dist];
this.generateFromArray(this.F);
end
%============================ barrier_constraint ============================
%
% @brief  Barrier Certificates & control generation
%
function [h, A, B, u] = barrier_constraint(this, x, u_des)
    
gamma = 0.3;
H = eye(2);

h = this.evaluate(x);   % Barrier Value

grad_h = this.gradient(x); % Gradient of Barrier

A = -grad_h';

B = gamma*h^3; %+ grad_ht;

opts = optimoptions(@quadprog, 'Display', 'off');
u = quadprog(H, -u_des, A, B, [], [], [], [], [], opts);

end

%============================ ROSinterface ============================
%
% @brief  Computations and initializations relevant to ROS and STDR.
%
function [x_state, x_rot, y_rot, x_rotplus, y_rotplus] = ROSinterface(this, sub, sub1)

    odom = receive(sub);                                         % 'receive' to get data from the subscriber.
    x_state(1,1)= odom.Pose.Pose.Position.X- this.grid_x;        % X Co-ordinate
    x_state(2,1)= odom.Pose.Pose.Position.Y- this.grid_y;        % Y Co-ordinate
    quat = odom.Pose.Pose.Orientation;
    angles = quat2eul([quat.W quat.X quat.Y quat.Z]);
    x_state(3,1) = wrapToPi(angles(1));                                    % Orientation

    laser = receive(sub1);         % 'receive' to get data from the subscriber.
    cart = laser.readCartesian;    % Convert Lidar data to cartesian coordinates
    
    if isempty(cart)
       x_rot = [];
       y_rot = [];
       x_rotplus = [];
       y_rotplus = [];
    else
        laser.Ranges = laser.Ranges - 0.04;
        cartplus = laser.readCartesian;
        [x_rot,y_rot] = this.DataTransform(cart,x_state); % Transformation of frames.
        [x_rotplus,y_rotplus] = this.DataTransform(cartplus,x_state);

%         plot(x_rot, y_rot, 'ro');  % Plot Lidar Data
%         axis([-1.6 1.6 -1 1])
%         % axis([0 this.grid_x 0 this.grid_y])
%         hold on;
%         plot(x_rotplus, y_rotplus, 'go');
%         axis([-1.6 1.6 -1 1])

    end
    
end
%========================== generateVelCmdsOmni ===========================
%
% @brief  Generates veolcity commands.
%
function generateVelCmdsOmni(this, u, pub)
    
%TODO : Sort out the input & output arguments.
msg = rosmessage(pub); % Empty message determined by the topic published by pub.
msg.Linear.X = u(1);                % Set Linear speed
msg.Linear.Y = u(2);               % Set Angular speed
send(pub,msg);                      % Send the message via the publisher
end
%=========================== generateVelCmdsUni ===========================
%
% @brief  Generates veolcity commands.
%
function generateVelCmdsUni(this, u, pub)
    
%TODO : Sort out the input & output arguments.
msg = rosmessage(pub); % Empty message determined by the topic published by pub.
msg.Linear.X = u(1);                % Set Linear speed
msg.Angular.Z = u(2);               % Set Angular speed
send(pub,msg);                      % Send the message via the publisher
end
%============================ shutdownROS ============================
%
% @brief  ROS Shutdown.
%
function shutdownROS(this, pub)
%TODO : Sort out the input & output arguments.
msg = rosmessage(pub); % Empty message determined by the topic published by pub.
msg.Linear.X = 0.0; % Set Velocities to 0 before shutting down.
msg.Linear.Y = 0.0;
msg.Linear.Z = 0.0;
msg.Angular.X = 0.0;
msg.Angular.Y = 0.0;
msg.Angular.Z = 0.0;
send(pub,msg);
rosshutdown;   % Shutting down global node
end
%============================ DataTransform ============================
%
% @brief  Transform data between 2 frames of references.
%
function [x_rot,y_rot] = DataTransform(this, cart, x_state)

cartf = cart + [x_state(1,1) x_state(2,1)];
cart_rnd = round(cartf,2);
x_cart = cart_rnd(:,1)';
y_cart = cart_rnd(:,2)';

x_center = x_state(1,1);
y_center = x_state(2,1);

% theta = pi + x_state(3,1);                % For Unicycle robot
theta = pi;                                 % For Omnidirectional robot

R = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];

a = [1 0 x_center;0 1 y_center; 0 0 1];
c = [1 0 -x_center;0 1 -y_center; 0 0 1];
M = a*R*c;
if size(cart_rnd,1) > 0
    for i=1:size(cart_rnd,1)
        rot(:,i) = M*[x_cart(i) y_cart(i) 1]';
    end
end
x_rot = rot(1, :);
y_rot = rot(2, :);

end
%============================ ROStf ============================
%
% @brief  Transform data between 2 frames of references.
%
function [x_rot,y_rot] = ROStf(this, cart,pt,tftree)

cart_rnd = round(cart,2);
x_cart = cart_rnd(:,1)';
y_cart = cart_rnd(:,2)';

if size(cart_rnd,1) > 0
    for i=1:size(cart_rnd,1)
        pt.Point.X = x_cart(i);
        pt.Point.Y = y_cart(i);
        pt.Point.Z = 0;
        tfpt = transform(tftree, 'map_static', pt);
        rot(1,i)= tfpt.Point.X;
        rot(2,i)= tfpt.Point.Y;
    end
end

x_rot = rot(1,:);
y_rot = rot(2,:);
x_rot = x_rot - this.grid_x;
y_rot = y_rot - this.grid_y;
end
%============================ evaluate ============================
%
% @brief  Barrier computation
%
% @param[in]  x   Location in the 2D world to evaluate barrier function at.
%
% The implementation accepts a single point or a column-wise matrix of
% points to evaluate at.
%
function h = evaluate(this, x_state)

h = double(this.F(x_state(1,:), x_state(2,:)));  

end
%============================= gradient ============================
%
% @brief  Compute the point-wise gradient of barrier function.
%
% @param[in]  x       Location(s) to evaluate \f$\nabla h\f$ at (column-wise).
% @param[out] gradH   Gradient evaluation (column-wise).
%
% The implementation accepts a single point or a column-wise matrix of
% points to evaluate at.
%
function [gradH] = gradient(this, x_state)

gradH = zeros(size(x_state));     % Because the gradient is also 2x1.

gradH(1,:) = double(this.dfx(x_state(1,:), x_state(2,:)));    
gradH(2,:) = double(this.dfy(x_state(1,:), x_state(2,:)));    

end

%============================ differential ===========================
%
% @brief  Compute the point-wise differential of the barrier function.
%
% @param[in]  x   State(s) to evaluate (column-wise).
%
% @param[out] dh  Differentials evaluated at indicated state points (row-wise).
%
function [dh] = differential(this,x)

dh = zeros(fliplr(size(x)));     % Because the gradient is also 2x1.

dh(:,1) = double(this.dfx(x(1,:), x(2,:)));    
dh(:,2) = double(this.dfy(x(1,:), x(2,:)));    

end
%============================ buildSTDR ===========================
%
% @brief  Builds the STDR example based on the changing grid.
%
function[gridFine,grid,bord_samp,shapeFine,shapeGrid] = buildSTDR(this,x_state)
    
  if (isfield(this.extras,'boundary') && this.extras.boundary)
    bRad = [];
    if (isfield(this.extras,'bRad'))
      bRad = this.extras.bRad;
    elseif (exist('sigma'))
      bRad = sqrt(this.sigma)/2;
    end
    if (~isempty('bRad'))
      % First, the horizontal boundaries.
      dgB = bRad;
      qbord = -this.grid_x:3.2/dgB:this.grid_x; % Modified for Robotarium
      bord_samp = [ [qbord; zeros(size(qbord))] , ...
                    [qbord; this.grid_x*ones(size(qbord))] ];
    
      % Second, the vertical boundaries.
      qbord = -this.grid_y:2/dgB:this.grid_y; % Modified for Robotarium
      qbord([1,end]) = [];
      bord_samp = [ bord_samp , ...
                    [ [zeros(size(qbord)) ; qbord] , ...
                      [this.grid_y*ones(size(qbord)) ; qbord] ] ];
    else
      disp('sigma not defined, so not providing set of boundary points.');
    end
  end

  if (isfield(this.extras,'grid') && this.extras.grid && exist('sigma')) % Specs chosen in testResidualRobotarium.m
    if (~isfield(this.extras, 'gridType'))
      this.extras.gridType = 'boundary';
      disp(['Default grid type is: ' this.extras.gridType]);
    end
    dgC = sqrt(this.sigma)/2;
    switch (this.extras.gridType)
      case {'inside','interior'}
        dOff = dgC/2;
        cOff = 0;
      case 'outside'
        dOff = -dgC/2;
        cOff = 2*0.01;
      case 'boundary'
        dOff = 0;
        cOff = 2*0.01;
    end
    grid{1} = linspace(-this.grid_x, this.grid_x, ...
                              cOff + round(2*this.grid_x/dgC) );
    grid{2} = linspace(-this.grid_y, this.grid_y, ...
                              cOff + round(2*this.grid_y/dgC) ); % Sizing the Grid
    shapeGrid = [size(grid{1},2), size(grid{2},2)];
  end

  if (isfield(this.extras,'gridFine') && this.extras.gridFine) 
    gridFine{1} = linspace( -this.grid_x, this.grid_x, round(2*this.grid_x/this.dg)+1 );
    gridFine{2} = linspace( -this.grid_y, this.grid_y, round(2*this.grid_y/this.dg)+1 );
    shapeFine = [size(gridFine{1},2), size(gridFine{2},2)];
  end

end
end
end
