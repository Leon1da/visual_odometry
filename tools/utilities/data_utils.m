source "tools/utilities/data_loader.m"
source "tools/utilities/geometry_helpers_3d.m"

function [XR_true, XR_guess] = pose_vector2transform();
  
  
  global pose_id;
  global odometry_pose;
  global ground_truth_pose;  
   
  num_poses = size(pose_id, 2);
  
  XR_true = zeros(4, 4, num_poses);
  XR_guess = zeros(4, 4, num_poses);
  
  for(i=1:num_poses)
    # convert x y theta vector (3d) to x y z alpha beta gamma vector (6d)
    ground_truth_pose_vector3d = ground_truth_pose(:, i);
    XR_true(:,:,i) = vector2transform(ground_truth_pose_vector3d); # convert gt_odom vector to transform
    
    # convert x y theta vector (3d) to x y z alpha beta gamma vector (6d)
    odom_pose_vector3d = odometry_pose(:, i);
    XR_guess(:,:,i) = vector2transform(odom_pose_vector3d); # convert odom vector to transform
    
  endfor;

endfunction;


function [XCam_true, XCam_guess] = camera_poses()  
  global camera_transformation;
  [XR_true, XR_guess] = pose_vector2transform();
  num_poses = size(XR_true, 3);
  XCam_true = zeros(4, 4, num_poses);
  XCam_guess = zeros(4, 4, num_poses);
  for i = 1:num_poses
    XCam_true(:,:, i) = XR_true(:,:, i) * camera_transformation;
    XCam_guess(:,:, i) = XR_guess(:,:, i) * camera_transformation;  
  endfor;
endfunction;
# transform a 3d vector (x, y, theta) to a homogeneous transform
function X = vector2transform(vector3);

# num_projection_measurements = total number of landmark measurements (all landmark image meas for all poses)
# Zp = list of all the (u,v) image coordinates of the landmarks of all the measurements
# projection_associations = list of all the (seq_id = pose_id, landmark_id)  of all the measurements
  # vector6 = [x, y, z, alpha, beta, gamma]
  # vector3 = [x, y, theta]
  vector6 = [vector3(1), vector3(2), 0, 0, 0, vector3(3)]'; # two months of madonnas
  X = v2t(vector6);
endfunction;

function [num_projection_measurements, Zp, projection_associations] = projection_measurements();
  global measurements_filenames;
  num_projection_measurements = 0; 
  Zp = [,];
  projection_associations = [,];
  for(i=1:size(measurements_filenames, 1))
    filename = measurements_filenames(i,:);
    [meas_seq, meas_ground_truth_odom_pose, meas_odom_pose, meas_point_id, meas_point_coord, meas_point_appearance] = read_measurement_information(filename); 
    for(j=1:size(meas_point_id, 2))
      num_projection_measurements++;
      Zp(num_projection_measurements, :) = meas_point_coord(:,j);
      projection_associations(num_projection_measurements, :) = [meas_seq, meas_point_id(j)];
    endfor;
  endfor;
  
  Zp = Zp';
  projection_associations = projection_associations';
endfunction;

function [num_projection_measurements, Zp, projection_associations] = per_pose_projection_measurements(sequences, indices=[]);
  global measurements_filenames;
  num_projection_measurements = 0; 
  Zp = [,];
  projection_associations = [,];
  
  for i = 1:size(sequences, 2)
    current_num_projection_measurements = 0; 
    current_Zp = [,];
    current_projection_associations = [,];
    seq = sequences(i);
    filename = measurements_filenames(seq,:); # get sequence filename
    [meas_seq, meas_ground_truth_odom_pose, meas_odom_pose, meas_point_id, meas_point_coord, meas_point_appearance] = read_measurement_information(filename); 
    for(j=1:size(meas_point_id, 2))
      if size(indices, 1) > 0 && find(indices == meas_point_id(j)) 
        current_num_projection_measurements++;
        current_Zp(current_num_projection_measurements, :) = meas_point_coord(:,j);
        current_projection_associations(current_num_projection_measurements, :) = [meas_seq, meas_point_id(j)];
      endif;
    endfor;
    num_projection_measurements = num_projection_measurements + current_num_projection_measurements; 
    Zp = [Zp; current_Zp];
    projection_associations = [projection_associations; current_projection_associations];

  endfor;
  
  Zp = Zp';
  projection_associations = projection_associations';
endfunction;

function [num_poses, pose_associations] = compute_pose_association();
  global pose_id;
  num_poses = size(pose_id, 2);
  pose_associations=zeros(2,num_poses - 1);
  for (i=1:num_poses-1)
      pose_associations(:,i)=[pose_id(:, i), pose_id(i)+1]';
  endfor;
endfunction;


# indices = set of landmark index
# landmark = set of landmarks whose index is in indices
function [landmarks] = get_landmark_by_indices(indices)
  global landmark_position; 
  num_indices = size(indices, 2);
  landmarks = zeros(3, num_indices);
  for(i = 1:num_indices)
    landmarks(:, i) = landmark_position(:, indices(i));
  endfor;
endfunction;
function [P1_img, P2_img, indices] = image_point_associations(idx1, idx2, coord1, coord2);
  for (i = 1:size(idx1, 2))
    for (j = 1:size(idx2, 2))  
      if( idx1(:,i) == idx2(:,j)) # same index
        P1_img(:, end+1) = coord1(:, i);
        P2_img(:, end+1) = coord2(:, j);
        indices(end+1) = idx1(:,i); # = idx2(:,j)
        continue;
      endif;
    endfor;
  endfor;
endfunction;
  
  
function [P1_img, P2_img, indices] = sequence_point_associations(start_seq_id, end_seq_id);
  global measurements_filenames;
  fn1 = measurements_filenames(start_seq_id,:);
  [~, ~, ~, idx1, coord1, ~] = read_measurement_information(fn1); 
  fn2 = measurements_filenames(end_seq_id,:);
  [~, ~, ~, idx2, coord2, ~] = read_measurement_information(fn2); 
  [P1_img, P2_img, indices] = image_point_associations(idx1, idx2, coord1, coord2);

endfunction;


function [meas_seq, meas_ground_truth_odom_pose, meas_odom_pose, meas_point_id, meas_point_coord, meas_point_appearance] = read_pose_measurement_information(pose_id)
  global measurements_filenames;
  measurement_filename = measurements_filenames(pose_id, :);
  [meas_seq, meas_ground_truth_odom_pose, meas_odom_pose, meas_point_id, meas_point_coord, meas_point_appearance] = read_measurement_information(measurement_filename);   
endfunction;

##function [P_img, indices] = images_points_associations(idx, coords, dims)
##  P = 0;
##  indices = 0;
##  
##  start_ = 1
##  end_ = dims(1);
##  idx1 = idx(start_:end_);
##  coord1 = coords(start_:end_));
##  
##  start_ = end_;
##  end_ = start_ + dims(2); 
##  idx2 = idx(start_:end_);
##  coord2 = coords(start_:end_));
##  
##  [P1_img, P2_img, indices] = image_point_associations(idx1, idx2, coord1, coord2);
##  
##  start = 1;
##  for i = 1:size(dims,2)
##    dim = dims(i);
##    idx1 = idx(1:dim);
##    [P1_img, P2_img, indices] = image_point_associations(idx1, idx2, coord1, coord2);
##    start = start + dim;
##  endfor;
##  
##endfunction;
##
##function [P_img, indices] = sequences_points_associations(sequences)
##  
  
##  global measurements_filenames;
##  dims = [];
##  idx = [];
##  coords = [];
##  i = 1;
##  for i = 1:size(sequences, 2)
##    seq = sequences(i);
##    fn = measurements_filenames(seq,:);
##    [~, ~, ~, m_indices, m_coords, ~] = read_measurement_information(fn); 
##    m_dim = size(m_indices, 2);
##    dims = [dims, m_dim];
##    idx = [idx, m_indices];
##    coords = [coords, m_coords];
##  endfor;
##  idx
##  coords
##  dims
##  [P_img, indices] = images_points_associations(idx, coords, dims)
##endfunction;  
  

  
  
  

################################## TEST FUNCTIONS #########################################

# test SE(2) to SE(3)
[XR_true, XR_guess] = pose_vector2transform();

# test pose associations
[num_poses, pose_associations] = compute_pose_association();

# test projection measurements
[num_projection_measurements, Zp, projection_associations] = projection_measurements();

### test landmark matching between two frames
##[P1_img, P2_img, indices] = sequence_point_associations(1, 2);
##
### given a set of landamark index return their position
##indices = [1, 4, 9];
##[landmarks] = get_landmark_by_indices(indices);





