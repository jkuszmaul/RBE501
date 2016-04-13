function H = to_htrans(t)
  % t: If 3x3, assumed to be a rotation matrix
  % and so [R 0
  %         0 1] is returned.
  % If 3x1, assumed to be translation vector
  % and so [I t
  %         0 1] is returned.
  H = t;
  if size(t) == [3 1]
    H = [eye(3) t; 0 0 0 1];
  else if size(t) == [3 3]
    H = [t zeros(3, 1); 0 0 0 1];
  end
end
