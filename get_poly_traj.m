% t: Total time it should take.
% traf_fun: Function taking time (0-t) and returning [q_d, dq_d, ddq_d].
function [traj_fun] = get_poly_traj(q0, dq0, qf, dqf, t)

  N = size(q0, 1);
  funs = cell(N, 1);
  for i = 1:size(q0, 1)
    funs{i} = traj_single(q0(i), dq0(i), qf(i), dqf(i), t);
  end

  traj_fun = @eval_traj;

  function [q_d, dq_d, ddq_d] = eval_traj(t)
    q_d = [];
    dq_d = [];
    ddq_d = [];

    for i = 1:N
      fun = funs{i};
      state = fun(t);
      q_d(end+1, 1) = state(1);
      dq_d(end+1, 1) = state(2);
      ddq_d(end+1, 1) = state(3);
    end
  end

  function fun = traj_single(p0, v0, pf, vf, tf)
    b=[p0;v0;pf;vf];

    A = [1 0 0 0;
         0 1 0 0;
         1 tf tf^2 tf^3;
         0 1 2*tf 3*tf^2];
    a = A\b;

    fun = @(t) [a(1) + a(2) * t + a(3) * t^2 + a(4) * t^3;
                a(2) + 2 * a(3) * t + 3 * a(4) * t^2;
                2 * a(3) + 6 * a(4) * t].';
  end
end
