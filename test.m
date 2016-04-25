function test()
  %[step, kin, jac, inv_vel, inv_dyn] = gen_test_model();
  [step, kin, jac, inv_vel, inv_dyn] = test_model_6dof();
  dt = 0.03;
  goal_vel = [5.1;-10;-15];
  goal_vel = [4;3; 4; 2; 4;-5];
  maxtau = 100;

  [goalp, goalv] = get_goal(goal_vel);

  q0 = zeros(3, 1);
  dq0 = zeros(3, 1);

  time_iter = 50;
  traj_fun = get_poly_traj(q0, dq0, goalp, goalv, time_iter * dt);

  Kp = eye(3);
  Kv = eye(3);
  num_iter = 200;
  %opts = optimset('MaxFunEvals', 5000 * 150);
  opts = optimoptions(@fmincon, 'MaxFunEvals', num_iter * 3 * time_iter,...
                      'TolX', 0.001, 'TolFun', 0.001,...
                      'Algorithm', 'interior-point');
  A = [eye(3); -eye(3)];
  b = maxtau * [ones(3); -ones(3)];
  u0 = gen_taus(q0, dq0, time_iter, dt);
  u0
  u = fmincon(@obj_fun, u0, [], [], [], [], -maxtau, maxtau, [], opts);
  [poses, vels, absvels] = run_sim(step, u, dt);
  vels
  absvels
  u

  absp = kin(poses(:, 1));
  Plot = plot3(absp(1, :), absp(2, :), absp(3, :), '-o');
  while true
    for pos = poses
      absp = kin(pos);
      set(Plot, 'XData', absp(1, :));
      set(Plot, 'YData', absp(2, :));
      set(Plot, 'ZData', absp(3, :));
      xlim([-3 3]);
      ylim([-3 3]);
      zlim([-2 2]);
      drawnow
    end
    pause(1)
  end

  function minval = obj_fun(u)
    [poses, vels, absvels] = run_sim(step, u, dt);
    minval = absvels(:, end) - goal_vel;
    minval;
    if isnan(minval)
      minval = 0
      return
    end
    minval = norm(minval)
    return
    for vel = absvels
      diff = vel - goal_vel;
      if norm(minval) > norm(diff)
        minval = norm(diff);
      end
    end
    minval
  end

  function w = wrapToPi(t)
    w = eval(t);
    while (max(w) > pi)
      w(w>pi) = w(w>pi) - 2*pi;
    end
    while (min(w) < -pi)
      w(w<-pi) = w(w<-pi) + 2*pi;
    end
  end

  function [q, dq] = get_goal(v)
  % Returns a valid goal pose for a goal velocity.
    [qv, dqv] = inv_vel(v)
    if size(symvar(dqv))
      q = [2;2;2];
      dq = subs(dqv, qv, q);
    else
      q = wrapToPi(qv)
      dq = wrapToPi(dqv)
    end
  end

  function [tau] = get_control(t, q, dq)
    [q_d, dq_d, ddq_d] = traj_fun(t);
    eq = q - q_d;
    edq = dq - dq_d;
    tau = inv_dyn(q, dq, ddq_d - Kp * eq - Kv * edq);
    tau = eval(tau);
  end

  function [taus] = gen_taus(q0, dq0, num_iter, dt)
    q = q0;
    dq = dq0;
    taus = [];
    for i = 1:num_iter
      tau = get_control(i * dt, q, dq);
      tau(tau > maxtau) = maxtau;
      tau(tau < -maxtau) = -maxtau;
      taus(:, i) = tau;
      [q, dq, ddq, vabs] = step(q, dq, taus(:, i), dt);
    end
    jac(q) * dq
  end
end

