source "./26_total_least_squares/total_least_squares.m"

function [XR, XL, chi_stats_p, num_inliers_p, H, b] = doProjLS(XR, XL,
	     Zp, projection_associations,
	     num_poses,
	     num_landmarks,
	     num_iterations,
	     damping,
	     kernel_threshold)
       
       global landmark_dim;
       
      ######################################## LANDMARK MEASUREMENTS ######################################## 
      # generate an ideal number of landmark measurements
      # each pose observes each landmark
      num_landmark_measurements=num_poses*num_landmarks;
      Zl=zeros(landmark_dim,num_landmark_measurements);
      landmark_associations=zeros(2,num_landmark_measurements);
      
      ######################################## POSE MEASUREMENTS ######################################## 

      # generate an odometry trajectory for the robot
      num_pose_measurements=num_poses-1;
      Zr=zeros(4,4,num_pose_measurements);
      pose_associations=zeros(2,num_pose_measurements);
      
      ## uncomment the following to suppress pose-landmark measurements
      Zl=zeros(3,0);

      # uncomment the following to suppress pose-pose measurements
      Zr=zeros(4,4,0);

       [XR, XL, chi_stats_l, num_inliers_l, chi_stats_p, num_inliers_p,chi_stats_r, num_inliers_r, H, b] = doTotalLS(XR, XL,
	     Zl, landmark_associations,
	     Zp, projection_associations,
	     Zr, pose_associations,
	     num_poses,
	     num_landmarks,
	     num_iterations,
	     damping,
	     kernel_threshold);
       
    
endfunction;

 
