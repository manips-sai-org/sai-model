// 08-derivatives:
// example of computing rbdl derivatives

#include <SaiModel.h>

#include <iostream>
#include <random>

using namespace std;

const string robot_fname = "resources/human.urdf";

int main(int argc, char** argv) {
	cout << "Loading robot file: " << robot_fname << endl;

    auto robot = new SaiModel::SaiModel(robot_fname);
    robot->updateModel();

    // kinematic hessian
    string link_name = "right_hand";
    Vector3d pos_in_link = Vector3d(0.1, 0.2, 0.3);
    std::vector<MatrixXd> dJdq = robot->getJacobianDerivative(link_name, pos_in_link);

    // mass matrix derivative
    std::vector<MatrixXd> dMdq = robot->getMassMatrixDerivative();

    // dynamic bias derivative (coriolis + centrifugal + gravity (where gravity is an option))
    bool eval_with_gravity = true;  // include gravity 
    MatrixXd dbdq = robot->getDynamicBiasDerivativeWrtQ(eval_with_gravity);  // wrt joint position
    MatrixXd dbddq = robot->getDynamicBiasDerivativeWrtDq(eval_with_gravity);  // wrt joint velocity

    // gravity derivative
    MatrixXd dgdq = robot->getGravityDerivative();
}