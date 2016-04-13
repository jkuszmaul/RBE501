[step, kin] = gen_test_model();
u = zeros(3, 3500);
dt = 0.03;
[poses, vels, absvels] = run_sim(step, u, dt);

absp = kin(poses(:, 1));
Plot = plot(absp(1, :), absp(2, :));
for pos = poses
  absp = kin(pos);
  set(Plot, 'XData', absp(1, :));
  set(Plot, 'YData', absp(2, :));
  xlim([-3 3]);
  ylim([-3 2]);
  drawnow
end
