from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "build" / "python"))

import sai_model_py as sm

puma_file = ROOT / "urdf_models" / "puma" / "puma.urdf"
kuka_file = ROOT / "urdf_models" / "iiwa7" / "kuka_iiwa.urdf"

puma_robot = sm.SaiModel(str(puma_file))
kuka_robot = sm.SaiModel(str(kuka_file))

print("\npuma degrees of freedom:", puma_robot.dof)
print("puma q:", puma_robot.q)
print("puma jacobian at end-effector:\n", puma_robot.J("end-effector"))

print("\nkuka degrees of freedom:", kuka_robot.dof)
print("kuka q:", kuka_robot.q)
print("kuka jacobian at end-effector:\n", kuka_robot.J("end-effector"))
