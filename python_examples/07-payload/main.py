from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "build" / "python"))

import sai_model_py as sm

robot_file = Path(__file__).resolve().parent / "resources" / "rprbot.urdf"
print(f"Loading robot file: {robot_file}")

robot = sm.SaiModel(str(robot_file))
robot.q = np.zeros(robot.dof)
robot.update_model()

print("gravity torques without load:", robot.joint_gravity_vector())

robot.add_load(
    "link2",
    1.0,
    np.array([0.0, 1.0, 0.0]),
    np.diag([0.2, 0.3, 0.4]),
    np.eye(4),
    "load",
)
print("gravity torques with load:", robot.joint_gravity_vector())

robot.remove_load("load")
print("gravity torques without load:", robot.joint_gravity_vector())
