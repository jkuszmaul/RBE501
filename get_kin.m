% params: A Nx4 matrix of the D-H parameters. May be symbolic. (d, theta, a, alpha).
% t: symbolic variable which params are a function of.
% poses: A 4x4xN matrix of homogeneous transformation matrices for each link.
function [poses, vels, jacobians] = get_kin(params, t)
  params = params(t);
  N = size(params, 1);
  N = N(1);
  if N < 1
    return
  end
  poses = cell(N, 1);

  poses{1} = dh_mat(params(1, 1), params(1, 2), params(1, 3), params(1, 4));
  for i = 2:N
    poses{i} = poses{i-1} * dh_mat(params(i, 1), params(i, 2), params(i, 3), params(i, 4));
  end
  vars = [t];
  isrot = zeros([6 1]);
  for i = 1:N
    if size(symvar(params(i, 1)))
      % d is variable
      vars(i, 1) = params(i, 1);
      isrot(i) = 0;
    else
      % theta is variable
      vars(i, 1) = params(i, 2);
      isrot(i) = 1;
    end
  end
  vels = cell(N, 1);
  for i = 1:N
    vels{i} = simplify(diff(poses{i}, t));
  end
  zs = sym('zs', [3 N]);
  for i = 1:N
    Ri = poses{i};
    Ri = Ri(1:3, 1:3);
    zs(:, i) = Ri * [0;0;1];%* diff(vars(i), t);
  end
  jacobians = cell(N, 1);
  for i = 1:N
    jacobian = sym('jacobian', [6 i]);
    if isrot(1)
      jacobian(1:3, 1) = cross([0;0;1], poses{i}(1:3, 4));
      jacobian(4:6, 1) = [0;0;1];
    else
      jacobian(1:3, 1) = zs(:, 1);
      jacobian(4:6, 1) = zeros(3, 1);
    end
    for j = 2:i
      if isrot(j)
        jacobian(1:3, j) = cross(zs(:, j-1), poses{i}(1:3, 4) - poses{j-1}(1:3, 4));
        jacobian(4:6, j) = zs(:, j-1);
      else
        jacobian(1:3, j) = transpose(zs(:, j-1));
        jacobian(4:6, j) = zeros(3, 1);
      end
    end
    jacobians{i} = simplify(jacobian);
  end
end
