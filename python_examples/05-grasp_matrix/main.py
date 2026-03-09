from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "build" / "python"))

import sai_model_py as sm

robot_file = Path(__file__).resolve().parent / "linkage.urdf"
linkage = sm.SaiModel(str(robot_file))

q = np.array([54.7356 / 180.0 * np.pi, 54.7356 / 180.0 * np.pi, 54.7356 / 180.0 * np.pi])
linkage.q = q
linkage.update_model()

linkage.add_environmental_contact(
    "link0", np.array([1.0, 0.0, 0.0]), np.eye(3), sm.ContactType.SurfaceContact
)
linkage.add_environmental_contact(
    "link1", np.array([1.0, 0.0, 0.0]), np.eye(3), sm.ContactType.PointContact
)

data = linkage.environmental_grasp_matrix_at_geometric_center()
print("--------------------------------------------")
print("2 contacts")
print("--------------------------------------------")
print("center point:", data.resultant_point)
print("Grasp matrix:\n", data.G)
print("Grasp matrix inverse:\n", data.G_inv)
print("R:\n", data.R)

linkage.add_environmental_contact(
    "link2", np.array([1.0, 0.0, 0.0]), np.eye(3), sm.ContactType.PointContact
)
data = linkage.environmental_grasp_matrix_at_geometric_center()
print("--------------------------------------------")
print("3 contacts")
print("--------------------------------------------")
print("center point:", data.resultant_point)
print("Grasp matrix:\n", data.G)
print("Grasp matrix inverse:\n", data.G_inv)
print("R:\n", data.R)
