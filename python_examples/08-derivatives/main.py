from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "build" / "python"))

import sai_model_py as sm

robot_file = Path(__file__).resolve().parent / "resources" / "human.urdf"
print(f"Loading robot file: {robot_file}")

robot = sm.SaiModel(str(robot_file))
robot.update_model()

link_name = "right_hand"
pos_in_link = np.array([0.1, 0.2, 0.3])

start_time = time.perf_counter()
dJdq = robot.get_jacobian_derivative(link_name, pos_in_link)
elapsed_ms = (time.perf_counter() - start_time) * 1000.0
print(f"Jacobian derivative computation time: {elapsed_ms:.6f} ms")

dMdq = robot.get_mass_matrix_derivative()
dbdq = robot.get_dynamic_bias_derivative_wrt_q(True)
dbddq = robot.get_dynamic_bias_derivative_wrt_dq(True)
dgdq = robot.get_gravity_derivative()

print("len(dJdq):", len(dJdq), "shape dJdq[0]:", dJdq[0].shape if dJdq else None)
print("len(dMdq):", len(dMdq), "shape dMdq[0]:", dMdq[0].shape if dMdq else None)
print("dbdq shape:", dbdq.shape)
print("dbddq shape:", dbddq.shape)
print("dgdq shape:", dgdq.shape)
