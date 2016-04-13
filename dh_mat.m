function mat = dh_mat(d, theta, a, alpha)
  mat = to_htrans([0 0 d]') * to_htrans(rotz(theta)) * to_htrans([a 0 0]') * to_htrans(rotx(alpha));
end
