from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "build" / "python"))

import sai_model_py as sm

robot_file = Path(__file__).resolve().parent / "rpspsbot.urdf"
print(f"Loading robot file: {robot_file}")
robot = sm.SaiModel(str(robot_file))

print("\nLinks:")
robot.display_links()
print("\nJoints:")
robot.display_joints()

print("\nSpherical joints:")
for joint in robot.spherical_joints():
    print(f"name: {joint.joint_name} - index: {joint.index} - w_index: {joint.w_index}")

print("\nq_size:", robot.q_size)
print("dof:", robot.dof)
print("q:", robot.q)
print("dq:", robot.dq)

print("\nJoint names:", robot.joint_names())
print("\nJoint name by q index:")
for i in range(robot.q_size):
    print(f"joint: {i} - name: {robot.joint_name(i)}")

eef_link_name = "link5"
print("\nEnd-effector position:", robot.position(eef_link_name))

sph_joint_name = robot.spherical_joints()[0].joint_name
axis = np.array([1.0, 0.0, 1.0])
axis = axis / np.linalg.norm(axis)
angle = np.pi / 6.0
quat_xyzw = np.array(
    [
        axis[0] * np.sin(angle / 2.0),
        axis[1] * np.sin(angle / 2.0),
        axis[2] * np.sin(angle / 2.0),
        np.cos(angle / 2.0),
    ]
)
robot.set_spherical_quat(sph_joint_name, quat_xyzw)
print("set spherical quat:", quat_xyzw)
print("model spherical quat:", robot.spherical_quat(sph_joint_name))
print("new q:", robot.q)

robot.update_kinematics()
print("new end-effector position:", robot.position(eef_link_name))
