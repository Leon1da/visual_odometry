global trajectory_filename = '06-VisualOdometry/data/trajectoy.dat';
global world_filename = '06-VisualOdometry/data/world.dat';
global camera_filename = '06-VisualOdometry/data/camera.dat';
global dirname = '06-VisualOdometry/data/';
global measurements_filenames;

# robot poses
global pose_id;
global odometry_pose;
global ground_truth_pose;

# landmark position
global landmark_id;
global landmark_position;
global landmark_appearance;

# camera parameters
global camera_matrix;
global camera_transformation;
global camera_z_near;
global camera_z_far; 
global camera_image_width;
global camera_image_height;


# pose_id, landmark_id, meas_seq, meas_point_id start from zero (0)
# in order to align data indices to octave data table I add one (1) to all the indices

# read trajectory information
function [pose_id, odometry_pose, ground_truth_pose]=read_robot_poses();
  global trajectory_filename;
  trajectory_data = dlmread(trajectory_filename, ' ', 0,1);
  pose_id = trajectory_data(:, 1)';
  pose_id = pose_id.+1;
  odometry_pose = trajectory_data(:,2:4)';
  ground_truth_pose = trajectory_data(:, 5:7)';
   
endfunction;

# read world information
function [landmark_id, landmark_position, landmark_appearance]=read_landmark_positions(); 
  global world_filename;
  world_data = dlmread(world_filename);
  landmark_id = world_data(:, 1)';
  landmark_id = landmark_id.+ 1;
  landmark_position = world_data(:, 2:4)';
  landmark_appearance = world_data(:, 5:14)';
endfunction;

# read camera information
function [camera_matrix, camera_transformation, camera_z_near, camera_z_far, camera_image_width, camera_image_height]=read_camera_parameters(); 
  global camera_filename;
  camera_data = dlmread(camera_filename);
  camera_matrix = camera_data(2:4, 1:3);
  camera_transformation = camera_data(6:9, 1:4);
  camera_z_near = camera_data(10, 2);
  camera_z_far = camera_data(11, 2);
  camera_image_width = camera_data(12, 2);
  camera_image_height = camera_data(13, 2); 
endfunction;



# read measurment information
function [meas_seq, meas_ground_truth_odom_pose, meas_odom_pose, meas_point_id, meas_point_coord, meas_point_appearance] = read_measurement_information(measurement_filename); 
  meas_seq = dlmread(measurement_filename, ' ', [0, 1, 0, 1]);
  meas_seq = meas_seq.+1;
  measurement_data_odom = dlmread(measurement_filename, ' ', [1, 1, 2, 3]);
  meas_ground_truth_odom_pose = measurement_data_odom(1, :)';
  meas_odom_pose = measurement_data_odom(2,:)';
  
  measurement_data_points = dlmread(measurement_filename, '', 3, 1);
  meas_point_id = measurement_data_points(:, 2)';
  meas_point_id = meas_point_id.+ 1;
  meas_point_coord = measurement_data_points(:, 3:4)'; # u,v col,row
  meas_point_appearance = measurement_data_points(:, 5:14)';

endfunction;


# return a list of measurement file names in dirname  
function [measurements_filenames] = get_measurement_filenames();
   
# given odometry pose and gt pose vectors (3-dim: x y theta) 
# return the corrispondent vector of homogeneous tranformation (4x4 matrix)
  global dirname;
  [filenames] = readdir(dirname);
  measurements_filenames = [""];
  for(i=1:size(filenames))
    name = filenames(i){1,1};
    find = strfind(name,'meas');
    if (find == 1)
      measurements_filenames(end+1,:) = strcat(dirname, name);
    endif;
  endfor;
endfunction;


###################################### TEST FUNCTIONS #########################################
# test read robot data
[pose_id, odometry_pose, ground_truth_pose] = read_robot_poses();

# test read landmark data
[landmark_id, landmark_position, landmark_appearance] = read_landmark_positions();

# test read camera file
[camera_matrix, camera_transformation, camera_z_near, camera_z_far, camera_image_width, camera_image_height] = read_camera_parameters(); 

# test get all measurement files
[measurements_filenames] = get_measurement_filenames();

# test read file measurements
num_measurements = size(measurements_filenames, 1);
for i = 1:num_measurements
  measurement_filename = measurements_filenames(i, :);
  [meas_seq, meas_ground_truth_odom_pose, meas_odom_pose, meas_point_id, meas_point_coord, meas_point_appearance] = read_measurement_information(measurement_filename); 
endfor
