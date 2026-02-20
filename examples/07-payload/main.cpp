// 07-payload:
// example of adding and removing a load

#include <SaiModel.h>

#include <iostream>
#include <random>

using namespace std;

const string robot_fname = "resources/rprbot.urdf";

int main(int argc, char** argv) {
	cout << "Loading robot file: " << robot_fname << endl;

    auto robot = new SaiModel::SaiModel(robot_fname);
    VectorXd q_init = VectorXd::Zero(robot->dof());
    robot->setQ(q_init);
    robot->updateModel();

    // gravity torques without load
    cout << "gravity torques without load: " << robot->jointGravityVector().transpose() << "\n";

	// add load to robot
    std::string load_name = "load";
    std::string ee_link = "link2";
	double mass = 1;
	Vector3d com = Vector3d(0, 1, 0);
	Matrix3d inertia = Vector3d(0.2, 0.3, 0.4).asDiagonal();
	robot->addLoad(ee_link, mass, com, inertia, Affine3d::Identity(), load_name);

    // gravity torques with load
    cout << "gravity torques with load: " << robot->jointGravityVector().transpose() << "\n";

    // remove load
    robot->removeLoad(load_name);

    // gravity torques without load
    cout << "gravity torques without load: " << robot->jointGravityVector().transpose() << "\n";

}