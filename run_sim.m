% model: As returned from gen_test_model.
% us: matrix where each column is the control inputs for a given timestep.
% Note: Initial states are assumed to be 0.
function [poses, vels, absvels] = run_sim(model, us, dt)
  poses = [];
  vels = [];
  absvels = [];
  N = size(us, 1);
  poses(:, 1) = zeros(N, 1);
  vels(:, 1) = zeros(N, 1);
  absvels(:, 1) = zeros(6, 1);
  for u = us
    [np, nv, na, vabs] = model(poses(:, end), vels(:, end), u, dt);
    poses(:, end+1) = np;
    vels(:, end+1) = nv;
    absvels(:, end+1) = vabs;
  end
end
