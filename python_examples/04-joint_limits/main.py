from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "build" / "python"))

import sai_model_py as sm

robot_file = Path(__file__).resolve().parent / "rpspbot.urdf"
print(f"Loading robot file: {robot_file}")
robot = sm.SaiModel(str(robot_file))

print("\nJoint limits are parsed automatically in the URDF file.\n")
for limit in robot.joint_limits():
    print(
        f"Joint name: {limit.joint_name} - joint index: {limit.joint_index} "
        f"- joint name: {robot.joint_name(limit.joint_index)} "
        f"- lower limit: {limit.position_lower} - upper limit: {limit.position_upper} "
        f"- velocity: {limit.velocity} effort: {limit.effort}"
    )
