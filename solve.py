#!/usr/bin/python
import numpy as np
import scipy.optimize

# 2-link arm parameters.
l1 = 1
l2 = 1
m1 = 1
m2 = 1
grav = 0
a = m1 * l1 ** 2 + m2 * l1 ** 2
b = m2 * l2 ** 2
c = m2 * l1 * l2
d = grav * (m1 * l1 + m2 * l1)
e = grav * m2 * l2

# Calculate dynamics for a given state.

def M(x):
  return np.matrix([
    [a + b + 2 * c * np.cos(x[1, 0]), b + c * np.cos(x[1, 0])],
    [b + c * np.cos(x[1, 0])        , b]])

def C(x):
  drag = 0 * np.matrix([[1., 0.],[0., 1.]])
  ret = drag + np.matrix([
    [-c * np.sin(x[1, 0]) * x[3, 0], -c * np.sin(x[1, 0]) * (x[2, 0] + x[3, 0])],
    [c * np.sin(x[1, 0]) * x[2, 0] , 0]])
  #print ret
  return ret

def G(x):
  return np.matrix([
    [-d * np.sin(x[0, 0]) - e * np.sin(x[0, 0] + x[1, 0])],
    [e * np.sin(x[0, 0] + x[1, 0])]])

# Cost function.
def func(inputs):
  dt = 0.01
  x = np.matrix([[0.],[0.],[0.],[0.]])
  goal_vel = np.matrix([[1.],[1.]])
  dist = 0.
  cost = dist
  print "New"
  print inputs
  for i in xrange(len(inputs) / 2):
    tau = np.matrix(inputs[2*i:2 * (i + 1)])
    tau = np.transpose(tau)
    #print tau
    M_inv = np.linalg.inv(M(x))
    #print M_inv
    qddot = M_inv * (tau - C(x) * x[2:4, :] - G(x))
    #print qddot.tolist()
    x[0:2, :] += x[2:4, :] * dt
    x[2:4, :] += qddot * dt
    dist = np.linalg.norm(goal_vel - x[2:4, :])
    #print x.tolist()
    if (dist < 1e-4):
      print x
      print cost
      return cost
    cost += dist
  print x
  #return dist * np.ones(inputs.shape)
  print cost
  return cost

#print scipy.optimize.root(func, np.zeros((200, 2)))
num_steps = 100
print scipy.optimize.fmin_slsqp(
    func, np.zeros(num_steps * 2),
    bounds = [(-10, 10), (-10, 10)] * num_steps,
    iter=1000)
