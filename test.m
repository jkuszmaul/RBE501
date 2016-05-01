function [u, vels, absvels] = test()
  [step, kin, jac, inv_vel, inv_dyn] = gen_test_model();
  goal_vel = .8 * [5;-10;0;0;0;-15];
  maxtau = [0;120;60];

  %[step, kin, jac, inv_vel, inv_dyn] = test_model_6dof();
  %goal_vel = 1.5*[4;3; 4; 2; 4;-5];
  %maxtau = [150;120;100];

  dt = 0.03;

  [goalp, goalv] = get_goal(goal_vel);

  num_joints = 3;

  q0 = zeros(num_joints, 1);
  dq0 = zeros(num_joints, 1);

  time_iter = 50;
  traj_fun = get_poly_traj(q0, dq0, goalp, goalv, time_iter * dt);

  lb = repmat(-maxtau, 1, time_iter);
  ub = repmat(maxtau, 1, time_iter);

  Kp = eye(3);
  Kv = eye(3);
  num_iter = 200;
  %opts = optimset('MaxFunEvals', 5000 * 150);
  opts = optimoptions(@fmincon, 'MaxFunEvals', num_iter * 3 * time_iter,...
                      'TolX', 0.004, 'TolFun', 0.004,...
                      'Algorithm', 'interior-point');
  u0 = gen_taus(q0, dq0, time_iter, dt);
  u0
  [u, fval] = fmincon(@obj_fun, u0, [], [], [], [], lb, ub, [], opts);
  [poses, vels, absvels] = run_sim(step, u, dt);
  vels
  absvels
  u
  if fval > 0.0001
    disp 'Failed to find solution'
    exit
  end
  plot(transpose(u))
  legend('Joint 1 Torque', 'Joint 2 Torque', 'Joint 3 Torque')
  xlabel('Timesteps')
  ylabel('Control inputs')
  pause

  absp = kin(poses(:, 1));
  plot3d = false
  Plot = 0
  if plot3d
    Plot = plot3(absp(1, :), absp(2, :), absp(3, :), '-o');
  else
    Plot = plot(absp(1, :), absp(2, :), '-o');
  end
  for i = 1:200
    for pos = poses
      absp = kin(pos);
      set(Plot, 'XData', absp(1, :));
      set(Plot, 'YData', absp(2, :));
      if plot3d
        set(Plot, 'ZData', absp(3, :));
        zlim([-2 2]);
      end
      xlim([-3 3]);
      ylim([-3 3]);
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
      q = 2 * ones(size(qv));
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
      for j = 1:size(tau, 1)
        if tau(j, :) > maxtau(j, :)
          tau(j, :) = maxtau(j, :);
        end
        if tau(j, :) < -maxtau(j, :)
          tau(j, :) = -maxtau(j, :);
        end
      end
      taus(:, i) = tau;
      [q, dq, ddq, vabs] = step(q, dq, taus(:, i), dt);
    end
    jac(q) * dq
  end
end

