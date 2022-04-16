clear;
source "tools/utilities/data_utils.m";
source "tools/projective_geometry/epipolar.m";
source "tools/visualization/visualization.m";


XCam_true = zeros(4, 4, num_poses);
XCam_guess = zeros(4, 4, num_poses);
for i = 1:num_poses
  XCam_true(:,:, i) = XR_true(:,:, i) * camera_transformation;
  XCam_guess(:,:, i) = XR_guess(:,:, i) * camera_transformation;  
endfor;


offset = 1;
start_seq_id = 1;
end_seq_id = start_seq_id + offset;

K = camera_matrix;
[P1_img, P2_img, indices] = sequence_point_associations(start_seq_id, end_seq_id);

[X_est] = estimateTransform(K, P1_img, P2_img, true);

P_ = get_landmark_by_indices(indices);

noise = 1;
if noise
  Xs = XCam_guess(:,:, start_seq_id)
  Xe = XCam_guess(:,:, end_seq_id)
else
  Xs = XCam_true(:,:, start_seq_id)
  Xe = XCam_true(:,:, end_seq_id)
endif;

Xe_est = Xs*X_est

X1 = Xs;
X2 = Xe_est;
##X1 = Xs;
##X2 = Xe;
[n_success, P, errors] = triangulatePoints3(K, X2, X1, P1_img, P2_img);

plot_points(P_, color="b*"); # true
plot_points(P, color="r*"); # triangulated

scale = P./P_;
