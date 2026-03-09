from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "build" / "python"))

import sai_model_py as sm

robot_file = Path(__file__).resolve().parent / "rprbot.urdf"
print(f"Loading robot file: {robot_file}")

robot = sm.SaiModel(str(robot_file))
ee_link = "link2"
ee_pos_in_link = np.array([0.0, 0.0, 1.0])

position = robot.position(ee_link, ee_pos_in_link)
velocity = robot.linear_velocity(ee_link, ee_pos_in_link)
rotation = robot.rotation(ee_link)
J = robot.J(ee_link, ee_pos_in_link)
gravity = robot.joint_gravity_vector()

print("\nInitial configuration")
print("q:", robot.q)
print("position:", position)
print("velocity:", velocity)
print("rotation:\n", rotation)
print("J:\n", J)
print("M:\n", robot.M)
print("gravity:", gravity)

new_q = np.array([np.pi / 2.0, 1.0, np.pi / 2.0])
new_dq = np.array([0.0, 1.0, np.pi / 12.0])
robot.q = new_q
robot.dq = new_dq

print("\nBefore update_kinematics()")
print("q:", robot.q)
print("position:", robot.position(ee_link, ee_pos_in_link))
print("M:\n", robot.M)

robot.update_kinematics()
print("\nAfter update_kinematics()")
print("position:", robot.position(ee_link, ee_pos_in_link))
print("velocity:", robot.linear_velocity(ee_link, ee_pos_in_link))
print("M:\n", robot.M)

robot.update_model()
print("\nAfter update_model()")
print("M:\n", robot.M)

J_task = J[2:3, :]
op = robot.operational_space_matrices(J_task)
print("\nOperational space matrices")
print("Lambda:\n", op.Lambda)
print("Jbar:\n", op.Jbar)
print("N:\n", op.N)
