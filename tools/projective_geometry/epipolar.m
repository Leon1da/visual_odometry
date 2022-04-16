source "tools/utilities/geometry_helpers_3d.m"

%projects a world point onto an image plane
function [p_img, is_outside] = projectPoint(K,p) 
  p_img=K*p;
  p_img*=1./p_img(3);
  r=2*K(1,3);
  c=2*K(2,3);
  is_outside=(p(3)<0 || p_img(1)<0 ||p_img(1)>c || p_img(2)<0
             || p_img(2)>r);
endfunction;


%projects a set of world points onto the image plane
function [P_img, is_outside] = projectPoints(K,P)
  P=K*P;
  P_img=P(1:2,:)./P(3,:);
  r=2*K(1,3);
  c=2*K(2,3);
  is_outside= ((P(3,:)<0)
  | (P_img(1,:)<0)
  | (P_img(1,:)>c)
  | (P_img(2,:)<0)
  | (P_img(2,:)>r));
endfunction;

function [P1, P2] = pruneProjectedPointPairs(P1, is_outside1,  P2, is_outside2)
  if (size(P1,2)!=size(P2,2))
    disp("points do not match");
  endif;
  n=size(P1,2);
  k=1;
  for (i=1:n)
    if (is_outside1(i) || is_outside2(i))
      continue;
    endif;
    P1(:,k)=P1(:,i);
    P2(:,k)=P2(:,i);
    ++k;
  endfor;
  P1=P1(:,1:k);
  P2=P2(:,1:k);
endfunction;

function [P1, P2, indices] = pruneProjectedPointPairsIndices(P1, is_outside1,  P2, is_outside2)
  if (size(P1,2)!=size(P2,2))
    disp("points do not match");
  endif;
  n=size(P1,2);
  k=1;
  for (i=1:n)
    if (is_outside1(i) || is_outside2(i))
      continue;
    endif;
    indices(:, k) = i;
    P1(:,k)=P1(:,i);
    P2(:,k)=P2(:,i);
    ++k;
  endfor;
  P1=P1(:,1:k);
  P2=P2(:,1:k);
endfunction;

function  E = transform2essential(X)
  E=X(1:3,1:3)'*skew(X(1:3,4));
endfunction;


function [X1, X2] = essential2transform(E)
  W=[0, -1,  0;
     1,  0,  0;
     0,  0,  1];

  [U,S,V]=svd(E);
  R1=V*W*U';
  if (det(R1)<0) #right handed condition
    [U,S,V]=svd(-E);
    R1=V*W*U';
  endif;
  #1st solution for the rotation
  X1=eye(4);
  X1(1:3,1:3)=R1;
  t_cross=R1*E;
  X1(1:3,4)= [t_cross(3,2)-t_cross(2,3);
              t_cross(1,3)-t_cross(3,1);
              t_cross(2,1)-t_cross(1,2)];

  #2nd solution for the rotation
  R2=V*W'*U';
  X2=eye(4);
  X2(1:3,1:3)=R2;
  t_cross=R2*E;
  X2(1:3,4)= [t_cross(3,2)-t_cross(2,3);
              t_cross(1,3)-t_cross(3,1);
              t_cross(2,1)-t_cross(1,2)];
endfunction;

function [X1,X2] = fundamental2transform(K,F)
  E=K'*F*K;
  [X1,X2] = essential2transform(E);
endfunction;

function F = transform2fundamental(K,X)
  iK = inv(K);
  F=iK'*transform2essential(X)*iK;
endfunction


#estimate fundamental matrix
function F = estimateFundamentalSimple(P1_img, P2_img)
  H=zeros(9,9);
  n_points=size(P1_img,2);
  for (i=1:n_points)
    p1_img=[P1_img(:,i); 1];
    p2_img=[P2_img(:,i); 1];
    A=reshape(p1_img*p2_img',1,9);
    H+=A'*A;
  endfor;
  [V,lambda]=eig(H);
  F=reshape(V(:,1),3,3);
endfunction

#computes a preconditioning matrix, that scales the points around the
#origin
#A(1:2,1:2): inverse sigma of the points
#A(1:2,3)  : -inverse sigma*mean
function A = computePreconditioningMatrix(P_img)
  n_points=size(P_img,2);
  P_img=[P_img; ones(1,n_points)];
  s=sum(P_img,2);
  s2=P_img*P_img';
  mu=s/n_points;
  sigma=s2/n_points-mu*mu';
  A=eye(3);
  A(1:2,1:2)=inv(chol(sigma(1:2,1:2)));
  A(1:2,3)=-A(1:2,1:2)*mu(1:2);
endfunction

#estimate fundamental matrix
function F = estimateFundamental(P1_img, P2_img, use_preconditioning=false)
  n_points=size(P1_img,2);
  A1=A2=eye(3);
  if (use_preconditioning)
    A1=computePreconditioningMatrix(P1_img);
    A2=computePreconditioningMatrix(P2_img);
  endif;
  AP1=A1(1:2,1:2)*P1_img+repmat(A1(1:2,3),1,n_points);
  AP2=A2(1:2,1:2)*P2_img+repmat(A2(1:2,3),1,n_points);
  Fa=estimateFundamentalSimple(AP1,AP2);
  F=A1'*Fa*A2;
endfunction


% triangulates a point, passing through two lines
% one passing through the origin, and having
% direction vector d1
% one passing through a point p2, and having
% direction d2
function [success, p, e] = triangulatePoint(p2, d1, d2)
  p=zeros(3,1);
  success=false;
  e=-1;
                      
  D=[-d1, d2];         #assemble system matrix to find ascissa 
  s=-(D'*D)\(D'*p2);          #s: ascissa of closest point on p1 and p2
  if (s(1)<0 || s(2)<0)
    return;
  endif;
  success=true;
  p1_triangulated=d1*s(1);   # point on 1st line
  p2_triangulated=d2*s(2)+p2; # point on 2nd line
  e=norm(p1_triangulated-p2_triangulated); #difference between the points
  p=0.5*(p1_triangulated+p2_triangulated);               #midpoint
endfunction;

# triangulates a batch of points in image coords,
# X: is the pose of the world w.r.t the 2nd camera
function [n_success, P, errors] = triangulatePoints(K,X,P1_img,P2_img)
  #initialize vars
  n_points=size(P1_img,2);
  P=zeros(3,n_points);
  errors=zeros(1,n_points);
  n_success= 0;
  
  #inverse transform
  iX=inv(X);
  #inverse camera matrix
  iK=inv(K);
  #inverse rotation * inverse camera matix
  iRiK=iX(1:3,1:3)*iK;

  #express the points in camera coordinates
  #rotate the direction vector of P2 in world frame
  D1_cam=iK*[P1_img; ones(1,n_points)];
  D2_cam=iRiK*[P2_img; ones(1,n_points)];
  p2=iX(1:3,4);
  for (i=1:n_points)
    p1_cam=D1_cam(:,i);
    p2_cam=D2_cam(:,i);
    [success, p, e]=triangulatePoint(p2, p1_cam, p2_cam);
    if (success==true)
      ++n_success;
      P(:,n_success)=p;
      errors(n_success)=e;
    endif;
  endfor;
endfunction;

#estimates a transform from a set of image projections
#with camera matrix K
function [X, P] = estimateTransform(K, P1_img, P2_img, use_preconditioning=false)  
 
  # compute fundamental
  F=estimateFundamental(P1_img,P2_img, use_preconditioning);

  # extract essential
  E=K'*F*K;
  
  #extract transforms from essential
  [X1,X2]=essential2transform(E);
  
  X=X1;
  n_in_front=0;
  #for each transform pick the best
  X_test=X1;
  [n_test, P_test]=triangulatePoints(K, X_test, P1_img, P2_img);
  if (n_test>n_in_front)
    X=X_test;
    n_in_front=n_test;
    P = P_test;
  endif;
  X_test(1:3,4)=-X_test(1:3,4);
  [n_test, P_test]=triangulatePoints(K, X_test,P1_img, P2_img);
  if (n_test>n_in_front)
    X=X_test;
    n_in_front=n_test;
    P = P_test;
  endif;
  X_test=X2;
  [n_test, P_test]=triangulatePoints(K, X_test,P1_img, P2_img);
  if (n_test>n_in_front)
    X=X_test;
    n_in_front=n_test;
    P = P_test;
  endif;
  X_test(1:3,4)=-X_test(1:3,4);
  [n_test, P_test]=triangulatePoints(K, X_test,P1_img, P2_img);
  if (n_test>n_in_front)
    X=X_test;
    n_in_front=n_test;
    P = P_test;
  endif;
##  n_in_front
endfunction;
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
function [X, P] = estimateTransform2(K, P1_img, P2_img)
  
  [E] = eightPoint(P1_img, P2_img, K, K);
  
  [X1,X2]=essential2transform(E);
    
  X=X1;
  n_in_front=0;
  #for each transform pick the best
  X_test=X1;
  [n_test, P_test]=triangulatePoints(K, X_test, P1_img, P2_img);
  if (n_test>n_in_front)
    X=X_test;
    n_in_front=n_test;
    P = P_test;
  endif;
  X_test(1:3,4)=-X_test(1:3,4);
  [n_test, P_test]=triangulatePoints(K, X_test,P1_img, P2_img);
  if (n_test>n_in_front)
    X=X_test;
    n_in_front=n_test;
    P = P_test;
  endif;
  X_test=X2;
  [n_test, P_test]=triangulatePoints(K, X_test,P1_img, P2_img);
  if (n_test>n_in_front)
    X=X_test;
    n_in_front=n_test;
    P = P_test;
  endif;
  X_test(1:3,4)=-X_test(1:3,4);
  [n_test, P_test]=triangulatePoints(K, X_test,P1_img, P2_img);
  if (n_test>n_in_front)
    X=X_test;
    n_in_front=n_test;
    P = P_test;
  endif;
##  n_in_front
endfunction;


function Nmatrix = getNormMat2d(x)
% Function: Compute the normalization matrix of 2d points in 
%           homogeneous coordinate
% Normalization criteria:
%   1. Move the centroidto the origin
%   2. Average distance of  points to the centroid is sqrt(2) 
% 
% Usage:
% 
%       Nmatrix = getNormMat(x)
%   where:
%       Nmatrix - the normalization matrix
%       x - input data, dim: 3xN
% 
% Institute: Australian National University
% Author: Zhen Zhang
% Last modified: 11 Apr. 2018

% Get the centroid
centroid = mean(x, 2);
% Compute the distance to the centroid
dist = sqrt(sum((x - repmat(centroid, 1, size(x, 2))) .^ 2, 1));
% Get the mean distance
mean_dist = mean(dist);
% Craft normalization matrix
Nmatrix = [sqrt(2) / mean_dist, 0, -sqrt(2) / mean_dist * centroid(1);...
           0, sqrt(2) / mean_dist, -sqrt(2) / mean_dist * centroid(2);...
           0, 0, 1];

endfunction; 
  
function [eMatrix] = eightPoint(P1_img, P2_img, K1, K2)
% Function Introdution:
% Given a set of correspondences between two images and the intrisic matrix
% of the calibrated camera for both views, compute the essential matrix
% associated with the epipolar geometry using eight points
%
% Inputs:
% matchedPoints1 - the coordinates of matched features in the first image,
%   expressed in inhomogeneous coordinate. The size is of Nx2, where N is the
%   number of matched points
% matchedPoints2 - same as above
% K1 - the intrisic matrix of the calibrated camera from the first view
% K2 - the intrisic matrix of the calibrated camera from the second view
%
% Outputs:
% eMatrix - the computed essential matrix
%
% Author: Frederic Zhang
% Last modified: 5 Jun. 2018
% Version: 3.0

% 8-point algorithm

% Error checking
matchedPoints1 = P1_img';
matchedPoints2 = P2_img';
[n1, c1] = size(matchedPoints1);
[n2, c2] = size(matchedPoints2);
if((c1 ~= 2) || (c2 ~= 2))
    error('Points are not formated with correct number of coordinates.');
endif;
if((n1 < 8) || (n2 < 8))
    error('There are not enough points to carry out the operation.');
endif;;

% Arrange data
p1 = transpose([matchedPoints1(1: 8, :), ones(8, 1)]);
p2 = transpose([matchedPoints2(1: 8, :), ones(8, 1)]);
norm1 = getNormMat2d(p1);
norm2 = getNormMat2d(p2);

% Normalisation
p1 = norm1 * p1;
p2 = norm2 * p2;

p1 = transpose(p1 ./ repmat(p1(3, :), [3, 1]));
p2 = transpose(p2 ./ repmat(p2(3, :), [3, 1]));

x1 = p1(:, 1);
y1 = p1(:, 2);
x2 = p2(:, 1);
y2 = p2(:, 2);

% Craft matrix A
A = [x2 .* x1, x2 .* y1, x2, y2 .* x1, y2 .* y1, y2, x1, y1, ones(8, 1)];
% Perform SVD
[~, ~, V] = svd(A);
fMatrix = [V(1, 9), V(2, 9), V(3, 9); V(4, 9), V(5, 9), V(6, 9); V(7, 9), V(8, 9), V(9, 9)];
% Obtain fundamental matrix
[U, S, V] = svd(fMatrix);
fMatrix = U(:, 1) * S(1,1) * transpose(V(:, 1)) + U(:, 2) * S(2,2) * transpose(V(:, 2));
fMatrix = norm2' * fMatrix * norm1;

% Return essential matrix
eMatrix = K2' * fMatrix * K1;

endfunction;


function [n_success, P, indices, errors] = triangulatePointsMultipleViews(K, Xs, Ps, Ps_)
  num_landmarks = size(Ps, 1);
  num_poses = size(Xs, 3);
  
  n_success = 0;
  P = [];
  errors = [];
  indices = [];
  
  iK = inv(K);
  for i = 1:num_landmarks
    points = [];
    directions = [];
    for j = 1:num_poses
      if Ps_(i, j) == 1 # check measurement
        X = Xs(:, :, j); # camera pose
        c = X(1:3, 4); # camera center
        p_img = Ps(i, j, :); # image point
        p_hom = [p_img(1); p_img(2); 1];
        p_cam = iK * p_hom;
        p_dir = X(1:3, 1:3) * p_cam;
        p_dir = normalize(p_dir);
        
        points(:, end+1) = c;
        directions(:, end+1) = p_dir;
      endif;
    endfor;
    n_points = size(points, 2);
    n_dirs = size(directions, 2);
    if n_points > 1 && n_points == n_dirs
      [success, p, e] = triangulateMultipleViews(points, directions);        
      if (success==true)
        ++n_success;
        P(:,n_success)=p;
        errors(:,n_success)=e;
        indices(n_success) = i;
      endif;
    endif;
  endfor;
endfunction;



# triangulate a batches of points appearing into two frames P1_img and P2_img
# each i-th element of P1_img and P2_img refer to the same point seen into two different image planes
# K is the camera matrix (intrinsic parameters)
# X1 and X2 refer to the two camera poses (in the world frame) respectively 
function [n_success, P, errors] = triangulatePoints3(K, X1, X2, P1_img, P2_img)
  if (size(P1_img,2) != size(P2_img,2))
    disp("points do not match");
  endif;
  
  n_points=size(P1_img,2);
 
  P=zeros(3,n_points);
  errors=zeros(3,n_points);
  n_success= 0;
  
  iK = inv(K);
  R1 =  X1(1:3,1:3);
  R2 =  X2(1:3,1:3);
  P1_cam = iK * [P1_img; ones(1, n_points)];
  P2_cam = iK * [P2_img; ones(1, n_points)];
  P1_directions = R1 * P1_cam;
  P2_directions = R2 * P2_cam;
  
  p1 = X1(1:3, 4);
  p2 = X2(1:3, 4);

  for (i=1:n_points)
    d1 = P1_directions(:,i);
    d2 = P2_directions(:,i);
    d1 = normalize(d1);
    d2 = normalize(d2);
    
    [success, p, e]=triangulatePointOrthogonal(p1, p2, d1, d2); # type 2
    points = [p1, p2];
    directions = [d1, d2];    
    [success, p, e] = triangulateMultipleViews(points, directions);# type 3 (multipe views)
##    [p_img, success] = projectPoint(K,p); # check reprojection

##    plot_rays(points, directions);
    
    if (success==true)
      ++n_success;
      P(:,n_success)=p;
      errors(:,n_success)=e;
    endif;

  endfor; 
endfunction;

function [success, p, e]  = triangulateMultipleViews(points, directions)
  A = zeros(3,3);
  B = zeros(3,1);
  P = zeros(3,1);

  for i = 1:size(points,2)
    a = directions(1,i);
    b = directions(2,i);
    c = directions(3,i);

    x = points(1,i);
    y = points(2,i);
    z = points(3,i);

    A(1,1) += 1 - a*a;
    A(1,2) += -a*b;
    A(1,3) += -a*c;
    A(2,2) += 1 - b*b;
    A(2,3) += -b*c;
    A(3,3) += 1 - c*c;

    B(1,1) += (1-a*a)*x - a*b*y - a*c*z;
    B(2,1) += -a*b*x + (1-b*b)*y - b*c*z;
    B(3,1) += -a*c*x - b*c*y + (1-c*c)*z;
  endfor
  
  % A is symmetric
  A(2,1) = A(1,2);
  A(3,1) = A(1,3);
  A(3,2) = A(2,3);

  % Solve linear system
  p = A\B;
  success = 1;
  e = 0;
  
endfunction

  

function [success, p, e] = triangulatePointOrthogonal(p1, p2, d1, d2)
  A = [-dot(d1,d1), dot(d2,d1); -dot(d1,d2), dot(d2,d2)];
  b = [dot(p1,d1)-dot(p2,d1); -dot(p2,d2)+dot(p1,d2)];
  s = inv(A)*b;
  success = 1;
  p1_triangulated=p1 + d1*s(1);   # point on 1st line
  p2_triangulated=p2 + d2*s(2); # point on 2nd line
  e=norm(p1_triangulated-p2_triangulated); #difference between the points
  p=0.5*(p1_triangulated+p2_triangulated);               #midpoint
endfunction;

function scale = estimateRelativeScale(K, X1, X2, P)
  
  num_points = size(P,2);
  idx = (1:num_points);
  idx_perm = nchoosek (idx, 2)';
  num_perm = size(idx_perm, 2);
  
  rss = [];
  for i = 1:num_perm
    p1 = [P(:, idx_perm(1, i)); 1];
    p2 = [P(:, idx_perm(2, i)); 1];
    
    p11 = X1*p1; # frame cam 1
    p21 = X1*p2; # frame cam 1
     
    p12 = X2*p1; # frame cam 2
    p22 = X2*p2; # frame cam 2
  
    d1 = norm(p11 - p21);
    d2 = norm(p12 - p22);
  
    rs = d1/d2;
    rss(i) = rs;
    
  endfor;
  
  scale = median(rss);
   
endfunction;


# paper implementation
function scale = estimateRelativeScale2(P1, P2, indices)
  
  lambda_x = P2(1, :) * P1(1, :)' / (norm(P2(1, :)) * norm(P2(1, :)))
  lambda_y = P2(2, :) * P1(2, :)' / (norm(P2(2, :)) * norm(P2(2, :)))
  lambda_z = P2(3, :) * P1(3, :)' / (norm(P2(3, :)) * norm(P2(3, :)))
  
  scale = lambda_z;
  
endfunction;


function scale = estimateRelativeScale3(P1, P2, indices)
  
  num_points = size(indices,2);
  idx = (1:num_points);
  idx_perm = nchoosek (idx, 2)';
  num_perm = size(idx_perm, 2);
  
  rss = [];
  for i = 1:num_perm
    k = idx_perm(1, i);
    j = idx_perm(2, i);
    
    scale = norm(P1(:, k)-P1(:,j))/norm(P2(:, k)-P2(:, j));
  
    rss(i) = scale;
  endfor;
  
  
  rss;
  scale = median(rss);
  
##  i = 10
##  scale = norm(P1(:, i)-P1(:, i + 1))/norm(P2(:, i)-P2(:, i + 1))
  

endfunction;

function scale = estimateRelativeScale4(P1, P2, indices, X)
  R = X(1:3, 1:3)
  T = X(1:3, 4)
  num_points = size(indices,2);
  M = zeros(num_points, num_points);
  
  for i = 1:num_points
##    for j = 1:num_points
      M(i, i) = P2(:, i)'*R*P1(:, i);
      M(i, num_points) = P2(:, i)'*T; 
##    endfor;
  endfor;
##  M
  [V,D] = eigs(M'*M, 1, 'SM')
  median(V(1:num_points-1))
##  [V, lambda] = eig(M'*M)
  scale = M;
##  i = 10
##  scale = norm(P1(:, i)-P1(:, i + 1))/norm(P2(:, i)-P2(:, i + 1))
  

endfunction;

# p = ground truth points position
# p_hat = estimated points position
function [s, R, t] = closed_form_scale_estimation(p, p_hat)

  N = size(p, 2);

  mu = sum(p, 2)/N; 
  mu_hat = sum(p_hat, 2)/N;


  sigma = 0;
  sigma_hat = 0;
  for i = 1:N
    sigma = sigma + norm(p(:, i) - mu)*norm(p(:, i) - mu);
    sigma_hat = sigma_hat + norm(p_hat(:, i) - mu_hat)*norm(p_hat(:, i) - mu_hat);
  endfor;
  sigma = sigma/N;
  sigma_hat = sigma_hat/N;

  SIGMA = 0;
  for i = 1:N
    SIGMA = SIGMA + (p(:, i) - mu)*(p_hat(:, i) - mu_hat)';
  endfor;
  SIGMA = SIGMA/N;

  [U, D, V] = svd(SIGMA);

  if det(U)*det(V') < 0
    W = diag([1, 1,-1]);
  else 
    W = eye(3);
  endif;

  R = U*W*V';
  s = trace(D*W)/sigma_hat;
  t = mu - s*R*mu_hat;

endfunction;
