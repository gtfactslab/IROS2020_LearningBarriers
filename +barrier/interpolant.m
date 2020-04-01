%========================== barrier.interpolant ==========================
%
% @class    barrier.interpolant
%
% @brief    A class which uses certain default functions of MATLAB for
%           computation of barrier function and its point-wise gradient.
%
% In the controls context here, barrier synthesis means to compute the
% barrier function and its gradient.
% 
%
%========================== barrier.interpolant ==========================

%
% @file     interpolant.m
%
% @author   Amogh J. Dabholkar,   adabholkar6@gatech.edu
% @date     2019/09/16
%
% @note
%   set indent to 2 spaces.
%   set tab to 4 spaces (with conversion to spaces).
%
% @classf   barrier
%========================== barrier.interpolant ==========================
classdef interpolant < barrier.base

properties (SetAccess = protected, GetAccess = public)
  F;         %< H(x)
  dfx;       %< dH/dx
  dfy;       %< dH/dy
  grid_x;    %< X grid domain
  grid_y;    %< Y Grid domain
  dg;        %< Grid Spacing

end

methods
   
  %================================ interpolant ===============================
  % 
  % @brief  Constructor for the interpolant barrier class.
  %
  % @param[in]  grid    Grid specification structure (x, y, dg).
  %
  function this = interpolant(grid)

  this.grid_x = grid.x; 
  this.grid_y = grid.y;
  this.dg = grid.dg;

  end

  %============================= generateFromArray =============================
  %
  % @brief  Define the interpolant.
  %
  % @param[in]  D   2D obstacle map defined over world domain.
  %
  function generateFromArray(this,D)

  if all(size(D) == size(this.grid_x))
    this.F = griddedInterpolant(this.grid_x',this.grid_y',D','cubic'); % h(x)

    [px,py] = gradient(D);

    this.dfx = griddedInterpolant(this.grid_x',this.grid_y',px'); % dh/dx 
    this.dfy = griddedInterpolant(this.grid_x',this.grid_y',py'); % dh/dy
  else
    error('Array dimensions do not align with the expected size')
  end
  
  end

  %============================= specifyByRadii =============================
  %
  % @brief  Define the distance field.
  %
  % @param[in]  spec    Specifications of obstacle centers and radii (4xNObs).
  %
  function specifyByRadii(this, spec)

  bwP = ((this.grid_x-spec(1,1))/spec(3,1)).^2 ...
                               + ((this.grid_y-spec(2,1))/spec(4,1)).^2 >= 1;

  if size(spec,2) > 1
    for i = 2:size(spec,2)                 % For n obstacles
      bwp = ((this.grid_x-spec(1,i))/spec(3,i)).^2 ...
                               + ((this.grid_y-spec(2,i))/spec(4,i)).^2 >= 1;
      bwP = bwP & bwp;
    end
  end

  Dp = bwdist(~bwP); 
  Dn = bwdist( bwP);
  D = (Dp-0.5*bwP) - (Dn-0.5*(~bwP));  % Distance field +outside, -inside.

  D = imgaussfilt(D, 2);        % Smooth to recover roundedness of boundary.
  D = D*this.dg;                % Factor in grid size.

  this.generateFromArray(D);
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
  function h = evaluate(this,x)
      
  h = double(this.F(x(1,:), x(2,:)));  
      
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
  function [gradH] = gradient(this,x_state)

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
  % @param[out] dh  Differentials evaluated at indicated points (row-wise).
  %
  function [dh] = differential(this,x)

  dh = zeros(fliplr(size(x)));     % Because the differential is 1x2.

  dh(:,1) = double(this.dfx(x(1,:), x(2,:)));    
  dh(:,2) = double(this.dfy(x(1,:), x(2,:)));    

  end
end


methods(Static)

  %============================ buildRobotarium ============================
  %
  % @brief  Builds the Robotarium world based on the grid spacing and 
  %         returns an instantiated barrier class.
  %
  % @param[in]  dg   Gridding resolution factor.
  %
  % @param[out] bar  Instantiated interpolant class.
  % @param[out] grid A structure consisting of the dimensions and grid spacing
  %                  of the Robotarium fixed world.
  %
  function [bar, grid] = buildRobotarium(dg)

    lx = -1.6:dg:1.6;               % Robotarium Domain in x and y 
    ly = -1.0:dg:1.0;
    [grid.x,grid.y] = meshgrid(lx,ly);

    grid.dg = dg;
    grid.lx = lx;
    grid.ly = ly;
    bar = barrier.interpolant(grid);

  end 
  
end

end

%
%========================== barrier.interpolant ==========================
