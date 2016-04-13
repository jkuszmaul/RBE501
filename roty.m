function rot = rotx(theta)
  rot = [cos(theta) 0 sin(theta);
         0 1 0;
         -sin(theta) 0 cos(theta)];
end
