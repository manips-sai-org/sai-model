// 09-muscle:
// example of computing muscle jacobian

#include <SaiModel.h>

#include <iostream>
#include <random>

#include <chrono>
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

using namespace std;

const string robot_fname = "resources/human.urdf";
const string muscle_fname = "resources/muscles_modified.xml";

int main(int argc, char** argv) {
	cout << "Loading robot file: " << robot_fname << endl;

    auto robot = new SaiModel::SaiModel(robot_fname);
    robot->updateModel();

    // parse muscles
    robot->addMuscleSystem(muscle_fname, "main");

    // compute muscle jacobian
    auto t1 = high_resolution_clock::now();
    auto L = robot->computeMuscleJacobian();
    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;
    std::cout <<  ms_double.count() << " ms\n";

    std::cout << "L matrix size: " << L.rows() << ", " << L.cols() << "\n";
    std::cout << "L norm: " << L.norm() << "\n";
}
