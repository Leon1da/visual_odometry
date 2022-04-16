close all
clear
clc

source "26_total_least_squares/projective_least_squares.m";
source "tools/projective_geometry/epipolar.m";
source "tools/utilities/data_utils.m";

# landmarks positions
[~, XL_true, ~] = read_landmark_positions();

# cameras pose
[XCam_true, XCam_guess] = camera_poses();

num_total_landmarks = size(XL_true, 2);
num_total_poses = size(XCam_true, 3);
num_total_measurements = 0;
measurement_dim = 2;

K = camera_matrix;

# init camera table
camera_table = zeros(4, 4, num_total_poses);
camera_table(:, :, 1) = eye(4)*camera_transformation; # first camera pose

# init point triangulated table
point_table = zeros(3, num_total_landmarks);
point_table_placeholder = zeros(1, num_total_landmarks);

# init measurement table
measurement_table = zeros(num_total_landmarks, num_total_poses, measurement_dim);
measurement_table_placeholder = zeros(num_total_landmarks, num_total_poses);

# read first measurement and fill the table
[meas_seq, ~, ~, meas_point_id, meas_point_coord, ~] = read_pose_measurement_information(1);
num_measurements = size(meas_point_id, 2);
num_total_measurements = num_total_measurements + num_measurements;
for i = 1:num_measurements
  landmark_id = meas_point_id(:, i);
  pose_id = meas_seq;
  measurement_coord = meas_point_coord(:, i);
  measurement_table(landmark_id, pose_id, :) = measurement_coord; 
  measurement_table_placeholder(landmark_id, pose_id) = 1;
endfor;

# read j-th measurement and fill the table
for j = 2:num_total_poses
  [meas_seq, ~, ~, meas_point_id, meas_point_coord, ~] = read_pose_measurement_information(j);
  P0 = [];
  P1 = [];
  idxs = [];
  
  num_curr_measurements = size(meas_point_id, 2);
  num_total_measurements = num_total_measurements + num_curr_measurements;
  for i = 1:num_curr_measurements
    landmark_id = meas_point_id(:, i);
    pose_id = meas_seq;
    measurement_coord = meas_point_coord(:, i);
    measurement_table(landmark_id, pose_id, :) = measurement_coord; 
    measurement_table_placeholder(landmark_id, pose_id) = 1;
    
    # get landmarks correspondences to estimate transform (eight points algorithm)
    if measurement_table_placeholder(landmark_id, pose_id - 1) == 1 # corrispondence found
      P0(:, end+1) = measurement_table(landmark_id, pose_id - 1, :);
      P1(:, end+1) = measurement_table(landmark_id, pose_id, :);
      idxs(end+1) = landmark_id;
    endif;
  endfor;
  # estimate transform using current and previous points clouds
  X_est = estimateTransform(K, P1, P0, true);
      
  absolute_scale = 18.0; # empirically determinaited
  T_est = X_est(1:3, 4)./absolute_scale;
  X_est(1:3, 4) = T_est;
    
  camera_table(:, :, j) = camera_table(:, :, j-1) * X_est;
 
  
  # triangulate points usign all measurements until now
  [n_success, P, indices, errors] = triangulatePointsMultipleViews(K, camera_table(:, :, 1:j), measurement_table, measurement_table_placeholder);
  P_ = get_landmark_by_indices(indices);
  rel_scale = P./P_;
  
  j
  
  for i = 1:n_success
    id = indices(i);
    point_table(:, id) = P(:, i);
    point_table_placeholder(id) = 1;  
  endfor;
  
 
  # mapping between real index and state index
  book_keeping_association = zeros(1, num_total_landmarks);
  index = 0;
  for n = 1:num_total_landmarks
    if point_table_placeholder(n) == 1 
      book_keeping_association(n) = ++index;    
    endif;
  endfor;
  
   
  Zp_ = Zp(:, num_total_measurements-num_measurements-num_curr_measurements + 1:num_total_measurements);
  projection_associations_ = projection_associations(:, num_total_measurements-num_measurements-num_curr_measurements + 1:num_total_measurements);
  num_measurements = num_curr_measurements
  
  Zp__ = [];
  projection_associations__ = [];
  for i = 1:size(projection_associations_, 2)
    pose_index = projection_associations_(1, i);
    landmark_index = projection_associations_(2, i);
    if book_keeping_association(landmark_index) != 0
      projection_associations__(:, end+1) = [pose_index-j+2, book_keeping_association(landmark_index)]';
      Zp__(:, end+1) = Zp_(:, i);
    endif;
  endfor;
  
  Zp_ = Zp__;
  projection_associations_ = projection_associations__;
  
  xl_true = P_;
  
  xr_guess = camera_table(:, :, j-1:j); 
  xl_guess = P;
  
  # optimization
  num_poses = size(xr_guess, 3)
  num_landmarks = size(xl_guess, 2)
  

  # one round PICP
  damping=0;
  kernel_threshold=1e3;
  num_iterations=1;

  [XR, XL, chi_stats_p, num_inliers_p, H, b] = doProjLS(xr_guess, xl_guess,
         Zp_, projection_associations_,
         num_poses,
         num_landmarks,
         num_iterations,
         damping,
         kernel_threshold);
  
  camera_table(:, :, j-1:j) = XR;
  
  
  # mapping between real index and state index
  index = 0;
  for n = 1:num_total_landmarks
    if book_keeping_association(n) != 0
      point_table(:, n) = XL(:, ++index);
    endif;
  endfor;
  
  # if last pose do bundle adjustment to refine the entire sequence (loop closure)
  if j == num_total_poses
    
    Zp_ = Zp(:, 1:num_total_measurements);
    projection_associations_ = projection_associations(:, 1:num_total_measurements);
    
    Zp__ = [];
    projection_associations__ = [];
    for i = 1:size(projection_associations_, 2)
      pose_index = projection_associations_(1, i);
      landmark_index = projection_associations_(2, i);
      if book_keeping_association(landmark_index) != 0
        projection_associations__(:, end+1) = [pose_index, book_keeping_association(landmark_index)]';
        Zp__(:, end+1) = Zp_(:, i);
      endif;
    endfor;
    
    Zp_ = Zp__;
    projection_associations_ = projection_associations__;
    
    xr_guess = camera_table; 
    xl_guess = XL;
        
    num_landmarks = size(xl_guess, 2)
    num_poses = size(xr_guess, 3)
    
    # do bundle adjustment
    damping=0;
    kernel_threshold=1e3;
    num_iterations=10;

    [XR, XL, chi_stats_p, num_inliers_p, H, b] = doProjLS(xr_guess, xl_guess,
           Zp_, projection_associations_,
           num_poses,
           num_landmarks,
           num_iterations,
           damping,
           kernel_threshold);
  endif;
  
endfor;



       
# Plot Landmarks
figure(1);
hold on;
grid;

subplot(1,2,1);
title("Landmarks");
plot3(xl_true(1,:), xl_true(2,:), xl_true(3,:),'b*',"linewidth",2);
hold on;
plot3(xl_guess(1,:), xl_guess(2,:), xl_guess(3,:),'ro',"linewidth",2);

subplot(1,2,2);
title("Landmarks after optimization");
plot3(xl_true(1,:), xl_true(2,:), xl_true(3,:),'b*',"linewidth",2);
hold on;
plot3(XL(1,:), XL(2,:), XL(3,:),'ro',"linewidth",2);

# Plot Poses
figure(3);
hold on;
grid;

subplot(1,2,1);
title("Poses");
plot3(XR_true(1,:), XR_true(2,:), XR_true(3,:),'b*',"linewidth",2);
hold on;
plot3(xr_guess(1,:), xr_guess(2,:), xr_guess(3,:),'ro',"linewidth",2);

subplot(1,2,2);
title("Poses after optimization");
plot3(XR_true(1,:), XR_true(2,:), XR_true(3,:),'b*',"linewidth",2);
hold on;
plot3(XR(1,:), XR(2,:), XR(3,:),'ro',"linewidth",2);


figure(5);
hold on;
grid;
title("chi evolution");

subplot(1,2,1);
plot(chi_stats_p, 'r-', "linewidth", 2);
legend("Chi Proj"); grid; xlabel("iterations");
subplot(1,2,2);
plot(num_inliers_p, 'b-', "linewidth", 2);
legend("#inliers");grid; xlabel("iterations");

figure(6);
title("H matrix");
H_ =  H./H;                      # NaN and 1 element
H_(isnan(H_))=0;                 # Nan to Zero
H_ = abs(ones(size(H_)) - H_);   # switch zero and one
H_ = flipud(H_);                 # switch rows
colormap(gray(64));
hold on;
image([0.5, size(H_,2)-0.5], [0.5, size(H_,1)-0.5], H_*64);
hold off;

