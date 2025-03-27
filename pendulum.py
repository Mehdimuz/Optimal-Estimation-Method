# Mehdi Muzaffari - MAE 674 Project 1

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

dt = 0.01   
num_steps = 10000
time = np.arange(0, num_steps) * dt

theta1 = np.pi/2
theta1_dot = 1
theta2 = np.pi/6
theta2_dot = 1

g = 9.81
l1 = 1.0
l2 = 1.0
m1 = 1.0
m2 = 1.0

def dynamics(theta1, theta1_dot, theta2, theta2_dot):
    
    theta1_ddot = (
    -g * (2 * m1 + m2) * np.sin(theta1)
    - m2 * g * np.sin(theta1 - 2 * theta2)
    - 2 * np.sin(theta1 - theta2) * m2 * (
    theta2_dot**2 * l2 + theta1_dot**2 * l1 * np.cos(theta1 - theta2))) / (l1 * (2 * m1 + m2 - m2 * np.cos(2 * theta1 - 2 * theta2)))

    theta2_ddot = (
    2 * np.sin(theta1 - theta2) * (
    theta1_dot**2 * l1 * (m1 + m2)
    + g * (m1 + m2) * np.cos(theta1)
    + theta2_dot**2 * l2 * m2 * np.cos(theta1 - theta2))) / (l2 * (2 * m1 + m2 - m2 * np.cos(2 * theta1 - 2 * theta2)))

    return theta1_ddot, theta2_ddot


theta1_history = np.zeros(num_steps)
theta2_history = np.zeros(num_steps)


theta1_history[0] = theta1
theta2_history[0] = theta2


for k in range(1, num_steps):

    theta1_ddot, theta2_ddot = dynamics(theta1, theta1_dot, theta2, theta2_dot)

    theta1_dot = theta1_dot + theta1_ddot * dt
    theta1 = theta1 + theta1_dot * dt

    theta2_dot = theta2_dot + theta2_ddot * dt
    theta2 = theta2 + theta2_dot * dt

    theta1_history[k] = theta1
    theta2_history[k] = theta2


x1 = l1 * np.sin(theta1_history)
y1 = -l1 * np.cos(theta1_history)
x2 = x1 + l2 * np.sin(theta2_history)
y2 = y1 - l2 * np.cos(theta2_history)


fig, ax = plt.subplots()
ax.set_xlim(-l1 - l2 - 0.5, l1 + l2 + 0.5)
ax.set_ylim(-l1 - l2 - 0.5, l1 + l2 + 0.5)
ax.set_aspect('equal')


line, = ax.plot([], [], '-o', lw=2, markersize=8, label="Pendulums")
trail1, = ax.plot([], [], 'blue', lw=0.5, label="Pendulum 1 Trajectory")
trail2, = ax.plot([], [], 'red', lw=0.5, label="Pendulum 2 Trajectory")


def init():
    line.set_data([], [])
    trail1.set_data([], [])
    trail2.set_data([], [])
    return line, trail1, trail2


def update(frame):
    current_x = [0, x1[frame], x2[frame]]
    current_y = [0, y1[frame], y2[frame]]

    line.set_data(current_x, current_y)

    trail1.set_data(x1[:frame], y1[:frame])
    trail2.set_data(x2[:frame], y2[:frame])

    return line, trail1, trail2

ani = animation.FuncAnimation(fig, update, frames=num_steps, init_func=init, interval=dt * 1000, blit=True)


plt.legend()
plt.show()


