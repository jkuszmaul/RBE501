function [taufun, accelfun] = get_dyn(params, t, masses, inertias)
  % We will do the following:
  % 1) Call get_kin to get the kinematics.
  % 2) Use the kinematics to get the velocities at each joint.
  % 3) Use the masses/moments of inertias to get the total energy.
  % 4) Assume that gravity is in the y-direction.
  % 5) Use the Euler-Lagrange equations.
  grav_axis = 2;
  [poses, vels, jacobians] = get_kin(params, t);
  g = sym('g');
  params = params(t);
  N = size(params, 1);
  if N < 1
    return
  end
  vars = [t];
  for i = 1:N
    if size(symvar(params(i, 1)))
      % d is variable
      vars(i, 1) = params(i, 1);
    else
      % theta is variable
      vars(i, 1) = params(i, 2);
    end
  end

  K = 0;
  P = 0;
  dvars = diff(vars, t);
  for i = 1:N
    vel = jacobians{i} * dvars(1:i);
    linvel = vel(1:3);
    angvel = vel(4:6);
    Klin = 0.5 * masses(i) * dot(linvel, linvel);
    Kang = 0.5 * angvel.' * inertias{i} * angvel;
    K = K + Klin + Kang;
    P = P + g * masses(i) * poses{i}(grav_axis, 4);
  end
  L = K - P;
  q = sym('q', [N 1]);
  q = sym(q, 'real');
  dq = sym('dq', [N 1]);
  dq = sym(dq, 'real');
  ddq = sym('ddq', [N 1]);
  ddq = sym(ddq, 'real');
  ddvars = diff(dvars, t);
  Lsubs = subs(subs(L, dvars, dq), vars, q);
  dLdq = subs(subs(jacobian(Lsubs, q), q, vars), dq, dvars);
  dLddq = subs(subs(jacobian(Lsubs, dq), q, vars), dq, dvars);
  tau = simplify(diff(dLddq, t) - dLdq);
  tau = reshape(tau, [], 1);
  M = simplify(subs(jacobian(subs(tau, ddvars, ddq), ddq), vars, q));
  G = subs(simplify(diff(tau, g) * g), vars, q);
  C = subs(subs(simplify(tau - M * ddvars - G), dvars, dq), vars, q);
  G = simplify(subs(G, g, 9.8));
  Minv = simplify(M^-1);
  tau = subs(tau, g, 9.8);
  tausym = sym('tausym', [N 1]);
  taufun = @(pos, vel, accel) subs(subs(subs(tau, ddvars, accel), dvars, vel), vars, pos);
  accelvar = simplify(Minv * (tausym - C - G));
  accelfun = matlabFunction(accelvar, 'Vars', {q, dq, tausym});
end
