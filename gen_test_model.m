function [step_fun, kin, jacfun, inv_vel_kin, inv_dyn] = gen_test_model()
  t = sym('t');
  t = sym(t, 'real');
  % Create params.
  % For simplicity, we use a fully-actuated planar arm.
  % This makes the problem easier to think about (although it is easy to
  % generalize it to 3 dimensions and 6 DOF).
  % The more complex issue is under-actuated arms.
  % Use a RRR planar arm.
  t1 = sym('t1(t)');
  t2 = sym('t2(t)');
  t3 = sym('t3(t)');
  params(t) = [0 t1 1 0;
               0 t2 1 0;
               0 t3 1 0];
  masses = [1;1;1];
  q = [t1; t2; t3];
  inertias = {eye(3),eye(3),eye(3)};
  B = 0.1;

  [poses, vels, jacobians] = get_kin(params, t);
  [inv_dyn, accel] = get_dyn(params, t, masses, inertias);
  step_fun = @sim_step;
  kin = @get_poses;
  qvar = sym('qvar', [3 1]);
  jac3 = subs(jacobians{3}, q, qvar);
  jacfun = matlabFunction(jac3, 'Vars', {qvar});
  inv_vel_kin = @inv_vel;

  function [qv, dqv] = inv_vel(v)
    v = [v(1:2, :); v(6)];
    dq = sym('dq', [3 1]);
    velsubs = simplify(jac3 * dq);
    % Because we are planar, ignore Z-vel and X/Y Rot.
    velsubs = [velsubs(1:2); velsubs(6)]
    %[qv1, qv2, qv3]
    [dqv1, dqv2, dqv3] = solve(velsubs == v, dq(1), dq(2), dq(3))
    dqv = [dqv1; dqv2; dqv3];
    qv = qvar;
  end

  function absp = get_poses(p)
    % Returns a matrix where each column is the x, y coordinates of
    % a given link.
    absp = [0;0;0];
    for i = 1:size(p, 1)
      pmat = subs(poses{i}, q, p);
      absp(:, end+1) = pmat(1:3, 4);
    end
  end

  function [np, nv, na, vabs] = sim_step(p, v, u, dt)
    tau = u - B*v;
    na = accel(p, v, tau);
    nv = v + na * dt;
    np = p + v * dt + 0.5 * na * dt^2;
    %vabs = subs(jacobians{3}, q, np) * nv;
    vabs = jacfun(np) * nv;
    % For 2D, we only care about X/Y translation and Z rotation.
    %vabs = [vabs(1:2);vabs(6)];
  end
end
