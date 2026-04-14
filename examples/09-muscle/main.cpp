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
const string muscle_fname = "resources/muscles_fixed.xml";

bool isNotFullRank(const Eigen::MatrixXd& A, double tol = 1e-9) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A);
    int rank = (svd.singularValues().array() > tol).count();
    return rank < std::min(A.rows(), A.cols());
}

double minSingularValue(const Eigen::MatrixXd& A) {
    if (A.size() == 0) {
        throw std::invalid_argument("Matrix is empty.");
    }

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const Eigen::VectorXd& s = svd.singularValues();

    return s(s.size() - 1);
}

struct ThinSVDComponents {
    MatrixXd U_s;
    VectorXd S_s;
    MatrixXd V_s;
};

ThinSVDComponents computeSingularComponents(
    const MatrixXd& A,
    double tol = 1e-6)
{
    Eigen::JacobiSVD<MatrixXd> svd(
        A, Eigen::ComputeThinU | Eigen::ComputeThinV);

    const VectorXd& S = svd.singularValues();

    // determine rank based on tolerance
    int r = 0;
    for (int i = 0; i < S.size(); ++i) {
        if (S(i) > tol) r++;
        else break;  // singular values are sorted
    }

    return {
        svd.matrixU().rightCols(r),
        S.tail(r),
        svd.matrixV().rightCols(r)
    };
}

int main(int argc, char** argv) {
	cout << "Loading robot file: " << robot_fname << endl;

    auto robot = new SaiModel::SaiModel(robot_fname);
    robot->updateModel();

    // parse muscles
    robot->addMuscleSystem(muscle_fname, "main");

    // get number of muscles
    int n_muscles = robot->getNumMuscles();

    // compute muscle jacobian
    {
        auto t1 = high_resolution_clock::now();
        auto L = robot->computeMuscleJacobian();
        std::cout << "L is not full rank status: " << isNotFullRank(L) << "\n";
        auto t2 = high_resolution_clock::now();
        duration<double, std::milli> ms_double = t2 - t1;
        std::cout <<  ms_double.count() << " ms\n";

        std::cout << "L matrix size: " << L.rows() << ", " << L.cols() << "\n";
        std::cout << "L norm: " << L.norm() << "\n";
        std::cout << "L min singular value: " << minSingularValue(L) << "\n";

        // svd
        auto svd_results = computeSingularComponents(L);
        std::cout << "task range: " << svd_results.S_s.transpose() << "\n";
    }

    // compute muscle jacobian derivative
    {
        auto t1 = high_resolution_clock::now();
        auto dLdq = robot->computeMuscleJacobianDerivative();
        auto t2 = high_resolution_clock::now();
        duration<double, std::milli> ms_double = t2 - t1;
        std::cout <<  ms_double.count() << " ms\n";
    }

    // compute muscle jacobian inverse 
    {
        auto t1 = high_resolution_clock::now();
        auto W = robot->computeMuscleCapacityMatrix();
        std::cout << "W is not full rank status: " << isNotFullRank(W) << "\n";
        auto L_inv = robot->computeMuscleJacobianInverse(W);
        auto t2 = high_resolution_clock::now();
        duration<double, std::milli> ms_double = t2 - t1;
        std::cout << ms_double.count() << " ms\n";
        std::cout << "L_inv norm: " << L_inv.norm() << "\n";
    }

    // compute muscle jacobian inverse derivative
    {
        auto t1 = high_resolution_clock::now();
        auto W = robot->computeMuscleCapacityMatrix();
        auto dLinvdq = robot->computeMuscleJacobianInverseDerivative(W);
        auto t2 = high_resolution_clock::now();
        duration<double, std::milli> ms_double = t2 - t1;
        std::cout <<  ms_double.count() << " ms\n";
    }
}
