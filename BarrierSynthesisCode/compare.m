close all;
clear all;
clc;

%%  Load all the trajectories
offlineSVM = load('offlineSVM_5obs_ic1.mat','Plt_data1');
onlineSVM = load('onlineSVM_5obs_ic1.mat','Plt_data1');
bwdist = load('bwdist_5obs_ic1.mat','Plt_data1');

TrajOffline = offlineSVM.Plt_data1;
TrajOnline = onlineSVM.Plt_data1;
Trajbwdist = bwdist.Plt_data1;

%% Compute the resampled Traj object
least = min([size(Trajbwdist,2) size(TrajOffline,2) size(TrajOnline,2)]);

Traj.bwdist(1,:) = interp1(Trajbwdist(1,:),linspace(1,size(Trajbwdist,2),least));
Traj.bwdist(2,:) = interp1(Trajbwdist(2,:),linspace(1,size(Trajbwdist,2),least));
Traj.Offline(1,:) = interp1(TrajOffline(1,:),linspace(1,size(TrajOffline,2),least));
Traj.Offline(2,:) = interp1(TrajOffline(2,:),linspace(1,size(TrajOffline,2),least));
Traj.Online(1,:) = interp1(TrajOnline(1,:),linspace(1,size(TrajOnline,2),least));
Traj.Online(2,:) = interp1(TrajOnline(2,:),linspace(1,size(TrajOnline,2),least));

%% Testing the resampling by plotting
p1 = plot(Traj.bwdist(1,:), Traj.bwdist(2,:), 'b-.', 'LineWidth', 3);
axis([-1.6 1.6 -1 1]);
hold on

p1 = plot(Traj.Offline(1,:), Traj.Offline(2,:), 'r-.', 'LineWidth', 3);
axis([-1.6 1.6 -1 1]);
hold on

p1 = plot(Traj.Online(1,:), Traj.Online(2,:), 'g-.', 'LineWidth', 3);
axis([-1.6 1.6 -1 1]);
hold on
title('Sampled Trajectories From Three Different Controllers', 'interpreter', 'latex');
set(gcf, 'color', 'w')
legend('BWdist Trajectory', 'Offline SVM Trajectory', 'Online SVM TRajectory', 'interpreter', 'latex');

%% Frechet Distance Calculation
cm2_off_bw = frechet(Trajbwdist(1,:)', Trajbwdist(2,:)', TrajOffline(1,:)', TrajOffline(2,:)');
cm2_on_bw = frechet(Trajbwdist(1,:)', Trajbwdist(2,:)', TrajOnline(1,:)', TrajOnline(2,:)');
cm2_off_on = frechet(TrajOnline(1,:)', TrajOnline(2,:)', TrajOffline(1,:)', TrajOffline(2,:)');

%% Correlation Coefficient Calculation
r_online = corr2(Traj.bwdist, Traj.Online);
r_offline = corr2(Traj.bwdist, Traj.Offline);
r_online_offline = corr2(Traj.Online, Traj.Offline);

grad_bwdist = gradient(Traj.bwdist);
grad_offline = gradient(Traj.Offline);
grad_online = gradient(Traj.Online);

d_online = max(norm(Traj.Online - Traj.bwdist)) + max(norm(grad_online - grad_bwdist));
display(d_online);

d_offline = max(norm(Traj.Offline - Traj.bwdist)) + max(norm(grad_offline - grad_bwdist));
display(d_offline);

d_offline_online = max(norm(Traj.Offline - Traj.Online)) + max(norm(grad_offline - grad_online));
display(d_offline_online);

%% Function definition - HausdorffDist
function [hd, D] = HausdorffDist(P,Q,lmf,dv)
% Calculates the Hausdorff Distance between P and Q
%
% hd = HausdorffDist(P,Q)
% [hd D] = HausdorffDist(P,Q)
% [hd D] = HausdorffDist(...,lmf)
% [hd D] = HausdorffDist(...,[],'visualize')
%
% Calculates the Hausdorff Distance, hd, between two sets of points, P and
% Q (which could be two trajectories). Sets P and Q must be matrices with
% an equal number of columns (dimensions), though not necessarily an equal
% number of rows (observations).
%
% The Directional Hausdorff Distance (dhd) is defined as:
% dhd(P,Q) = max p c P [ min q c Q [ ||p-q|| ] ].
% Intuitively dhd finds the point p from the set P that is farthest from
% any point in Q and measures the distance from p to its nearest neighbor
% in Q.
% 
% The Hausdorff Distance is defined as max{dhd(P,Q),dhd(Q,P)}
%
% D is the matrix of distances where D(n,m) is the distance of the nth
% point in P from the mth point in Q.
%
% lmf: If the size of P and Q are very large, the matrix of distances
% between them, D, will be too large to store in memory. Therefore, the
% function will check your available memory and not build the D matrix if
% it will exceed your available memory and instead use a faster version of
% the code. If this occurs, D will be returned as the empty matrix. You may
% force the code to forgo the D matrix even for small P and Q by calling the
% function with the optional 3rd lmf variable set to 1. You may also force
% the function to return the D matrix by setting lmf to 0. lmf set to []
% allows the code to automatically choose which mode is appropriate.
%
% Including the 'vis' or 'visualize' option plots the input data of
% dimension 1, 2 or 3 if the small dataset algorithm is used.
%
% Performance Note: Including the lmf input increases the speed of the
% algorithm by avoiding the overhead associated with checking memory
% availability. For the lmf=0 case, this may represent a sizeable
% percentage of the entire run-time.
%
%

% %%% ZCD Oct 2009 %%%
% Edits ZCD: Added the matrix of distances as an output. Fixed bug that
%   would cause an error if one of the sets was a single point. Removed
%   excess calls to "size" and "length". - May 2010
% Edits ZCD: Allowed for comparisons of N-dimensions. - June 2010
% Edits ZCD: Added large matrix mode to avoid insufficient memory errors
%   and a user input to control this mode. - April 2012
% Edits ZCD: Using bsxfun rather than repmat in large matrix mode to
%   increase performance speeds. [update recommended by Roel H on MFX] -
%   October 2012
% Edits ZCD: Added a plotting function for visualization - October 2012
%

sP = size(P); sQ = size(Q);

if ~(sP(2)==sQ(2))
    error('Inputs P and Q must have the same number of columns')
end

if nargin > 2 && ~isempty(lmf)
    % the user has specified the large matrix flag one way or the other
    largeMat = lmf;     
    if ~(largeMat==1 || largeMat==0)
        error('3rd ''lmf'' input must be 0 or 1')
    end
else
    largeMat = 0;   % assume this is a small matrix until we check
    % If the result is too large, we will not be able to build the matrix of
    % differences, we must loop.
    if sP(1)*sQ(1) > 2e6
        % ok, the resulting matrix or P-to-Q distances will be really big, lets
        % check if our memory can handle the space we'll need
        memSpecs = memory;          % load in memory specifications
        varSpecs = whos('P','Q');   % load in variable memory specs
        sf = 10;                    % build in a saftey factor of 10 so we don't run out of memory for sure
        if prod([varSpecs.bytes]./[sP(2) sQ(2)]) > memSpecs.MaxPossibleArrayBytes/sf
            largeMat = 1;   % we have now concluded this is a large matrix situation
        end
    end
end

if largeMat
% we cannot save all distances, so loop through every point saving only
% those that are the best value so far

maxP = 0;           % initialize our max value
% loop through all points in P looking for maxes
for p = 1:sP(1)
    % calculate the minimum distance from points in P to Q
    minP = min(sum( bsxfun(@minus,P(p,:),Q).^2, 2));
    if minP>maxP
        % we've discovered a new largest minimum for P
        maxP = minP;
    end
end
% repeat for points in Q
maxQ = 0;
for q = 1:sQ(1)
    minQ = min(sum( bsxfun(@minus,Q(q,:),P).^2, 2));
    if minQ>maxQ
        maxQ = minQ;
    end
end
hd = sqrt(max([maxP maxQ]));
D = [];
    
else
% we have enough memory to build the distance matrix, so use this code
    
% obtain all possible point comparisons
iP = repmat(1:sP(1),[1,sQ(1)])';
iQ = repmat(1:sQ(1),[sP(1),1]);
combos = [iP,iQ(:)];

% get distances for each point combination
cP=P(combos(:,1),:); cQ=Q(combos(:,2),:);
dists = sqrt(sum((cP - cQ).^2,2));

% Now create a matrix of distances where D(n,m) is the distance of the nth
% point in P from the mth point in Q. The maximum distance from any point
% in Q from P will be max(D,[],1) and the maximum distance from any point
% in P from Q will be max(D,[],2);
D = reshape(dists,sP(1),[]);

% Obtain the value of the point, p, in P with the largest minimum distance
% to any point in Q.
vp = max(min(D,[],2));
% Obtain the value of the point, q, in Q with the largets minimum distance
% to any point in P.
vq = max(min(D,[],1));

hd = max(vp,vq);

end





% --- visualization ---
if nargin==4 && any(strcmpi({'v','vis','visual','visualize','visualization'},dv))
    if largeMat == 1 || sP(2)>3
        warning('MATLAB:actionNotTaken',...
            'Visualization failed because data sets were too large or data dimensionality > 3.')
        return;
    end
    % visualize the data
    figure
    subplot(1,2,1)
    hold on
    axis equal
    
    % need some data for plotting
    [mp ixp] = min(D,[],2);
    [mq ixq] = min(D,[],1);
    [mp ixpp] = max(mp);
    [mq ixqq] = max(mq);
    [m ix] = max([mq mp]);
    if ix==2
        ixhd = [ixp(ixpp) ixpp];
        xyf = [Q(ixhd(1),:); P(ixhd(2),:)];
    else
        ixhd = [ixqq ixq(ixqq)];
        xyf = [Q(ixhd(1),:); P(ixhd(2),:)];
    end
    
    % -- plot data --
    if sP(2) == 2
        h(1) = plot(P(:,1),P(:,2),'bx','markersize',10,'linewidth',3);
        h(2) = plot(Q(:,1),Q(:,2),'ro','markersize',8,'linewidth',2.5);
        % draw all minimum distances from set P
        Xp = [P(1:sP(1),1) Q(ixp,1)];
        Yp = [P(1:sP(1),2) Q(ixp,2)];
        plot(Xp',Yp','-b');
        % draw all minimum distances from set Q
        Xq = [P(ixq,1) Q(1:sQ(1),1)];
        Yq = [P(ixq,2) Q(1:sQ(1),2)];
        plot(Xq',Yq','-r');
        
        % denote the hausdorff distance
        h(3) = plot(xyf(:,1),xyf(:,2),'-ks','markersize',12,'linewidth',2);
        uistack(fliplr(h),'top')
        xlabel('Dim 1'); ylabel('Dim 2');
        title(['Hausdorff Distance = ' num2str(m)])
        legend(h,{'P','Q','Hausdorff Dist'},'location','best')
        
    elseif sP(2) == 1   
        ofst = hd/2;    % plotting offset
        h(1) = plot(P(:,1),ones(1,sP(1)),'bx','markersize',10,'linewidth',3);
        h(2) = plot(Q(:,1),ones(1,sQ(1))+ofst,'ro','markersize',8,'linewidth',2.5);
        % draw all minimum distances from set P
        Xp = [P(1:sP(1)) Q(ixp)];
        Yp = [ones(sP(1),1) ones(sP(1),1)+ofst];
        plot(Xp',Yp','-b');
        % draw all minimum distances from set Q
        Xq = [P(ixq) Q(1:sQ(1))];
        Yq = [ones(sQ(1),1) ones(sQ(1),1)+ofst];
        plot(Xq',Yq','-r');
        
        % denote the hausdorff distance
        h(3) = plot(xyf(:,1),[1+ofst;1],'-ks','markersize',12,'linewidth',2);
        uistack(fliplr(h),'top')
        xlabel('Dim 1'); ylabel('visualization offset');
        set(gca,'ytick',[])
        title(['Hausdorff Distance = ' num2str(m)])
        legend(h,{'P','Q','Hausdorff Dist'},'location','best')
        
    elseif sP(2) == 3
        h(1) = plot3(P(:,1),P(:,2),P(:,3),'bx','markersize',10,'linewidth',3);
        h(2) = plot3(Q(:,1),Q(:,2),Q(:,3),'ro','markersize',8,'linewidth',2.5);
        % draw all minimum distances from set P
        Xp = [P(1:sP(1),1) Q(ixp,1)];
        Yp = [P(1:sP(1),2) Q(ixp,2)];
        Zp = [P(1:sP(1),3) Q(ixp,3)];
        plot3(Xp',Yp',Zp','-b');
        % draw all minimum distances from set Q
        Xq = [P(ixq,1) Q(1:sQ(1),1)];
        Yq = [P(ixq,2) Q(1:sQ(1),2)];
        Zq = [P(ixq,3) Q(1:sQ(1),3)];
        plot3(Xq',Yq',Zq','-r');
        
        % denote the hausdorff distance
        h(3) = plot3(xyf(:,1),xyf(:,2),xyf(:,3),'-ks','markersize',12,'linewidth',2);
        uistack(fliplr(h),'top')
        xlabel('Dim 1'); ylabel('Dim 2'); zlabel('Dim 3');
        title(['Hausdorff Distance = ' num2str(m)])
        legend(h,{'P','Q','Hausdorff Dist'},'location','best')
        

    end
    
    subplot(1,2,2)
    % add offset because pcolor ignores final rows and columns
    [X Y] = meshgrid(1:(sQ(1)+1),1:(sP(1)+1));
    hpc = pcolor(X-0.5,Y-0.5,[[D; D(end,:)] [D(:,end); 0]]);
    set(hpc,'edgealpha',0.25)
    xlabel('ordered points in Q (o)')
    ylabel('ordered points in P (x)')
    title({'Distance (color) between points in P and Q';...
        'Hausdorff distance outlined in white'})
    colorbar('location','SouthOutside')
    
    hold on
    % bug: does not draw when hd is the very last point
    rectangle('position',[ixhd(1)-0.5 ixhd(2)-0.5 1 1],...
        'edgecolor','w','linewidth',2);
    
end
end

%% Function definition- Frechet distance
function f = frechet(X1,Y1,X2,Y2,varargin)

%get path point length
L1=length(X1);
L2=length(X2);

%check vector lengths
if or(L1~=length(Y1),L2~=length(Y2))
    error('Paired input vectors (Xi,Yi) must be the same length.')
end

%check for column inputs
if or(or(size(X1,1)<=1,size(Y1,1)<=1),or(size(X2,1)<=1,size(Y2,1)<=1))
    error('Input vectors must be column vectors.')
end

%create maxtrix forms
X1_mat=ones(L2,1)*X1';
Y1_mat=ones(L2,1)*Y1';
X2_mat=X2*ones(1,L1);
Y2_mat=Y2*ones(1,L1);

%calculate frechet distance matrix
frechet1=sqrt((X1_mat-X2_mat).^2+(Y1_mat-Y2_mat).^2);
fmin=min(frechet1(:));
fmax=max(frechet1(:));

%handle resolution
if ~isempty(varargin)
    res=varargin{1};
    if res<=0
        error('The resolution parameter must be greater than zero.')
    elseif ((fmax-fmin)/res)>10000
        warning('Given these two curves, and that resolution, this might take a while.')
    elseif res>=(fmax-fmin)
        warning('The resolution is too low given these curves to compute anything meaningful.')
        f=fmax;
        return
    end
else
    res=(fmax-fmin)/1000;
end

%compute frechet distance
for q3=fmin:res:fmax
    im1=bwlabel(frechet1<=q3);
    
    %get region number of beginning and end points
    if and(im1(1,1)~=0,im1(1,1)==im1(end,end))
        f=q3;
        break
    end
end

end
%% Function definition - FrechetDist
function [cm, cSq] = DiscreteFrechetDist(P,Q,dfcn)
% Calculates the discrete Frechet distance between curves P and Q
%
% [cm, cSq] = DiscreteFrechetDist(P,Q)
% [cm, cSq] = DiscreteFrechetDist(...,dfcn)
%
% P and Q are two sets of points that define polygonal curves with rows of
% vertices (data points) and columns of dimensionality. The points along
% the curves are taken to be in the order as they appear in P and Q.
%
% Returned in cm is the discrete Frechet distance, aka the coupling
% measure, which is zero when P equals Q and grows positively as the curves
% become more dissimilar.
%
% The optional dfcn argument allows the user to specify a function with
% which to calculate distance between points in P and Q. If not provided,
% the L2 norm is used.
%
% The secondary output, cSq, is the coupling sequence, that is, the
% sequence of steps along each curve that must be followed to achieve the
% minimum coupling distance, cm. The output is returned in the form of a
% matrix with column 1 being the index of each point in P and column 2
% being the index of each point in Q. (NOTE: the coupling sequence is not
% unique in general)
%
% Explanation:
% The Frechet distance is a measure of similarity between to curves, P and
% Q. It is defined as the minimum cord-length sufficient to join a point
% traveling forward along P and one traveling forward along Q, although the
% rate of travel for either point may not necessarily be uniform.
%
% The Frechet distance, FD, is not in general computable for any given
% continuous P and Q. However, the discrete Frechet Distance, also called
% the coupling measure, cm, is a metric that acts on the endpoints of
% curves represented as polygonal chains. The magnitude of the coupling
% measure is bounded by FD plus the length of the longest segment in either
% P or Q,  and approaches FD in the limit of sampling P and Q.
%
% This function implements the algorithm to calculate discrete Frechet
% distance outlined in:
% T. Eiter and H. Mannila. Computing discrete Frechet distance. Technical
% Report 94/64, Christian Doppler Laboratory, Vienna University of
% Technology, 1994.
% 
%
%
% EXAMPLE:
% % create data
% t = 0:pi/8:2*pi;
% y = linspace(1,3,6);
% P = [(2:7)' y']+0.3.*randn(6,2);
% Q = [t' sin(t')]+2+0.3.*randn(length(t),2);
% [cm, cSq] = DiscreteFrechetDist(P,Q);
% % plot result
% figure
% plot(Q(:,1),Q(:,2),'o-r','linewidth',3,'markerfacecolor','r')
% hold on
% plot(P(:,1),P(:,2),'o-b','linewidth',3,'markerfacecolor','b')
% title(['Discrete Frechet Distance of curves P and Q: ' num2str(cm)])
% legend('Q','P','location','best')
% line([2 cm+2],[0.5 0.5],'color','m','linewidth',2)
% text(2,0.4,'dFD length')
% for i=1:length(cSq)
%   line([P(cSq(i,1),1) Q(cSq(i,2),1)],...
%        [P(cSq(i,1),2) Q(cSq(i,2),2)],...
%        'color',[0 0 0]+(i/length(cSq)/1.35));
% end
% axis equal
% % display the coupling sequence along with each distance between points
% disp([cSq sqrt(sum((P(cSq(:,1),:) - Q(cSq(:,2),:)).^2,2))])
%
%
% 
% %%% ZCD June 2011 %%%
% %%% edits ZCD May 2013: 1) remove excess arguments to internal functions
% and persistence for speed, 2) added example, 3) allowed for user defined
% distance function, 4) added aditional output option for coupling sequence
%


% size of the data curves
sP = size(P);
sQ = size(Q);

% check validity of inputs
if sP(2)~=sQ(2)
    error('Curves P and Q must be of the same dimension')
elseif sP(1)==0
    cm = 0;
    return;
end

% initialize CA to a matrix of -1s
CA = ones(sP(1),sQ(1)).*-1;

% distance function
if nargin==2
    dfcn = @(u,v) sqrt(sum( (u-v).^2 ));
end

% final coupling measure value
cm = c(sP(1),sQ(1));

% obtain coupling measure via backtracking procedure
if nargout==2
    cSq = zeros(sQ(1)+sP(1)+1,2);    % coupling sequence
    CApad = [ones(1,sQ(1)+1)*inf; [ones(sP(1),1)*inf CA]];  % pad CA
    Pi=sP(1)+1; Qi=sQ(1)+1; count=1;  % counting variables
    while Pi~=2 || Qi~=2
        % step down CA gradient
        [v,ix] = min([CApad(Pi-1,Qi) CApad(Pi-1,Qi-1) CApad(Pi,Qi-1)]);
        if ix==1
            cSq(count,:) = [Pi-1 Qi];
            Pi=Pi-1;
        elseif ix==2
            cSq(count,:) = [Pi-1 Qi-1];
            Pi=Pi-1; Qi=Qi-1;
        elseif ix==3
            cSq(count,:) = [Pi Qi-1];
            Qi=Qi-1;
        end
        count=count+1;
    end
    % format output: remove extra zeroes, reverse order, subtract off
    % padding value, and add in the last point
    cSq = [flipud(cSq(1:find(cSq(:,1)==0,1,'first')-1,:))-1; sP(1) sQ(1)];
end


% debug
% assignin('base','CAw',CA)

function CAij = c(i,j)
    % coupling search function
    if CA(i,j)>-1
        % don't update CA in this case
        CAij = CA(i,j);
    elseif i==1 && j==1
        CA(i,j) = dfcn(P(1,:),Q(1,:));     % update the CA permanent
        CAij = CA(i,j);                    % set the current relevant value
    elseif i>1 && j==1
        CA(i,j) = max( c(i-1,1), dfcn(P(i,:),Q(1,:)) );
        CAij = CA(i,j);
    elseif i==1 && j>1
        CA(i,j) = max( c(1,j-1), dfcn(P(1,:),Q(j,:)) );
        CAij = CA(i,j);
    elseif i>1 && j>1
        CA(i,j) = max( min([c(i-1,j), c(i-1,j-1), c(i,j-1)]),...
            dfcn(P(i,:),Q(j,:)) );
        CAij = CA(i,j);
    else
        CA(i,j) = inf;
    end
end     % end function, c

end     % end main function
