from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "build" / "python"))

import sai_model_py as sm

robot_file = Path(__file__).resolve().parent / "pbot.urdf"
print(f"Loading robot file: {robot_file}")

robot = sm.SaiModel(str(robot_file), True)
print("\nrobot degrees of freedom:", robot.dof)
print("robot coordinates:", robot.q)

mass_params = robot.get_link_mass_params("link0")
print("link 0 mass properties:")
print("mass:", mass_params.mass)
print("center of mass in link:", mass_params.com_pos)
print("inertia tensor:\n", mass_params.inertia)
