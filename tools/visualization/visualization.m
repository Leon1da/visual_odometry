# plots a set of rays passing from point (points(:, i)) and having direction (directions(:, i))
function plot_rays(points, directions, color='black', scale=10)
  num_points = size(points, 2);
  num_directions = size(directions, 2);
  
  if num_points != num_directions
    disp("Error plot_rays: dimension mismatch.")
  endif;
  
  for i = 1:num_points
    p = points(:, i);
    d = directions(:, i);
    q = quiver3(p(1), p(2), p(3), d(1), d(2), d(3));
    hold on;
    set(q, 'Color', color);
    set(q, 'ShowArrowHead', 'off');
    set(q, 'AutoScale', 'on');
    set(q, 'AutoScaleFactor', scale); 
    
  endfor;

endfunction;

# plots a set of points (points)
function plot_points(points, color='m*')
  num_points = size(points, 2);
  for i = 1:num_points  
    p = points(:, i);
    plot_point(p, color);
  endfor;
endfunction;

function plot_point(point, color="m*") 
    plot3(point(1), point(2), point(3), color,"linewidth", 1);
    hold on;
endfunction;

function plot_reference_frames(rfs)
  num_rfs = size(rfs, 3);
  for i = 1:num_rfs
    rf = rfs(:, :, i);
    plot_reference_frame(rf);
  endfor;
endfunction;
  
function plot_reference_frame(rf)
    
  t = rf(1:3, 4);
  R = rf(1:3,1:3);
  x = R(1:3, 1);
  y = R(1:3, 2);
  z = R(1:3, 3);

  plot_rays([t], [x], 'red', 1);
  plot_rays([t], [y], 'green', 1);
  plot_rays([t], [z], 'blue', 1);
  
endfunction;
