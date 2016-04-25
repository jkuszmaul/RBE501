function [step_fun, kin, jacfun, inv_vel_kin, inv_dyn] = test_model_6dof()
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
  params(t) = [0 t1 0 pi/2;
               0 t2 1 0;
               0 t3 1 0];
  masses = [1;1;1];
  q = [t1; t2; t3];
  inertias = {eye(3),eye(3),eye(3)};


  [poses, vels, jacobians] = get_kin(params, t);
  [inv_dyn, accel] = get_dyn(params, t, masses, inertias);
  step_fun = @sim_step;
  kin = @get_poses;
  qvar = sym('qvar', size(q));
  jac_end = subs(jacobians{size(q, 1)}, q, qvar);
  jacfun = matlabFunction(jac_end, 'Vars', {qvar});
  inv_vel_kin = @inv_vel;

  function [qv, dqv] = inv_vel(v)
    dq = sym('dq', size(qvar));
    velsubs = simplify(jac_end * dq);
    ans = solve(velsubs == v, 'PrincipalValue', true);
    qv = [ans.qvar1(1); ans.qvar2(1); ans.qvar3(1)];
    dqv = [ans.dq1(1); ans.dq2(1); ans.dq3(1)];
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
    na = accel(p, v, u);
    nv = v + na * dt;
    np = p + v * dt + 0.5 * na * dt^2;
    %vabs = subs(jacobians{3}, q, np) * nv;
    vabs = jacfun(np) * nv;
  end
end
