close all
clear
clc
source "tools/utilities/data_utils.m";
source "tools/projective_geometry/epipolar.m";
source "tools/visualization/visualization.m";


XCam_true = zeros(4, 4, num_poses);
XCam_guess = zeros(4, 4, num_poses);
for i = 1:num_poses
  XCam_true(:,:, i) = XR_true(:,:, i) * camera_transformation;
  XCam_guess(:,:, i) = XR_guess(:,:, i) * camera_transformation;  
endfor;


start_seq_id = 1;
end_seq_id = start_seq_id + 1;

noise = 0;
if noise
  Xs = XCam_guess(:,:, start_seq_id);
  Xe = XCam_guess(:,:, end_seq_id);
else
  Xs = XCam_true(:,:, start_seq_id);
  Xe = XCam_true(:,:, end_seq_id);
  # plot_reference_frame(Xs);
endif;

K = camera_matrix;
[P1_img, P2_img, indices] = sequence_point_associations(start_seq_id, end_seq_id);


[X_est, P_est] = estimateTransform(K, P2_img, P1_img, true);

P_true = get_landmark_by_indices(indices);

Xs
Xe


absolute_scale = 50.0
X1 = Xs
X_est(1:3, 4) = X_est(1:3, 4)./absolute_scale
X2 = X1*X_est

[n_success, P, errors] = triangulatePoints3(K, X1, X2, P1_img, P2_img);

plot_points(P_true, color="b*"); # true
plot_points(P, color="ro"); # triangulated
##plot_points(P_est, color="go"); # estimated

scale = P./P_true;




