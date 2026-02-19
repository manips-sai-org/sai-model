/**
 * @file main.cpp
 * @author William Chong (wmchong@stanford.edu)
 * @brief 
 * @version 0.1
 * @date 2026-02-17
 * 
 * @copyright Copyright (c) 2026
 * 
 */

#include "Sai2Model.h"
#include <adrbdl/adrbdl.h>

#include <random>

#include <chrono>
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

const std::string robot_fname = "./resources/full_body_rev2.urdf";

template<typename Derived>
Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1>
randomVectorBetween(const Eigen::MatrixBase<Derived>& minVec,
                    const Eigen::MatrixBase<Derived>& maxVec)
{
    using Scalar = typename Derived::Scalar;

    if (minVec.size() != maxVec.size())
        throw std::invalid_argument("Vectors must have the same size.");

    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> result(minVec.size());

    static std::random_device rd;
    static std::mt19937 gen(rd());

    for (int i = 0; i < minVec.size(); ++i)
    {
        std::uniform_real_distribution<Scalar> dist(minVec[i], maxVec[i]);
        result[i] = dist(gen);
    }

    return result;
}

int main() {

    auto robot = std::make_shared<Sai2Model::Sai2Model>(robot_fname);
    auto ad_robot = std::make_shared<AutoDiffRigidBodyDynamics::Model>(robot_fname);

    const std::string link_name = "right_hand";
    const Vector3d pos_in_link = Vector3d(0.1, 0.2, 0.3);

    auto joint_limits = robot->jointLimits();
    VectorXd q_lower(robot->dof()), q_upper(robot->dof());
    for (int i = 0; i < robot->dof(); ++i) {
        q_lower(i) = joint_limits[i].position_lower;
        q_upper(i) = joint_limits[i].position_upper;
    }

    if (true) {
        // test jacobian derivative
        int n_samples = 1e4;
        double ed_avg = 0;
        double ad_avg = 0;
        for (int i = 0; i < n_samples; ++i) {
            VectorXd q_rand = randomVectorBetween(q_lower, q_upper);
            robot->setQ(q_rand);
            robot->updateKinematics();
            auto t1 = high_resolution_clock::now();
            std::vector<MatrixXd> dJdq_test = robot->getJacobianDerivative(link_name, pos_in_link);
            auto t2 = high_resolution_clock::now();
            duration<double, std::milli> ms_double = t2 - t1;
            ed_avg += (1./n_samples) * ms_double.count();
            
            VectorXdual q_rand_dual = q_rand.cast<dual>();
            t1 = high_resolution_clock::now();
            std::vector<MatrixXd> dJdq_truth = AutoDiffRigidBodyDynamics::jacobianDerivative(ad_robot, q_rand_dual, link_name, pos_in_link, AutoDiffRigidBodyDynamics::BOTH);
            t2 = high_resolution_clock::now();
            ms_double = t2 - t1;
            ad_avg += (1./n_samples) * ms_double.count();

            double error = 0;
            for (int j = 0; j < robot->dof(); ++j) {
                error += (dJdq_test[j] - dJdq_truth[j]).norm();
            }
            if (error > 1e-8) {
                throw runtime_error("dJdq doesn't match");
            }
        }

        std::cout << "ed avg: " << ed_avg << "\n";
        std::cout << "ad avg: " << ad_avg << "\n";
        std::cout << "dJdq test passed\n";
    }

    if (true)
    {
        // test mass matrix derivative
        int n_samples = 1e4;
        double ed_avg = 0;
        double ad_avg = 0;
        for (int i = 0; i < n_samples; ++i) {
            VectorXd q_rand = randomVectorBetween(q_lower, q_upper);
            robot->setQ(q_rand);
            robot->updateModel();
            
            auto t1 = high_resolution_clock::now();
            std::vector<MatrixXd> dMdq_test = robot->getMassMatrixDerivative();
            auto t2 = high_resolution_clock::now();
            duration<double, std::milli> ms_double = t2 - t1;
            ed_avg += (1./n_samples) * ms_double.count();
            
            t1 = high_resolution_clock::now();
            VectorXdual q_rand_dual = q_rand.cast<dual>();
            std::vector<MatrixXd> dMdq_truth = AutoDiffRigidBodyDynamics::massMatrixDerivative(ad_robot, q_rand_dual);
            t2 = high_resolution_clock::now();
            ms_double = t2 - t1;
            ad_avg += (1./n_samples) * ms_double.count();

            double error = 0;
            for (int j = 0; j < robot->dof(); ++j) {
                error += (dMdq_test[j] - dMdq_truth[j]).norm();
            }
            if (error > 1e-8) {
                std::cout << i << "\n";
                for (int j = 0; j < robot->dof(); ++j) {
                    if ((dMdq_test[j] - dMdq_truth[j]).norm() > 1e-8) {
                        std::cout << "j: " << j << "\n";
                        // std::cout << dMdq_test[j] << "\n---\n";
                        // std::cout << dMdq_truth[j] << "\n---\n";
                        break;
                    }
                }
                throw runtime_error("dMdq doesn't match");
            }
        }

        std::cout << "ed avg: " << ed_avg << "\n";
        std::cout << "ad avg: " << ad_avg << "\n";
        std::cout << "dMdq test passed\n";
    }
    
    {
        // test dynamic bias derivatives
        // test mass matrix derivative
        int n_samples = 1e4;
        double ed_avg = 0;
        double ad_avg = 0;
        for (int i = 0; i < n_samples; ++i) {
            VectorXd q_rand = randomVectorBetween(q_lower, q_upper);
            VectorXd dq_rand = VectorXd::Random(robot->dof());
            robot->setQ(q_rand);
            robot->setDq(dq_rand);
            robot->updateModel();
            
            auto t1 = high_resolution_clock::now();
            auto [dbdq_test, dbddq_test] = robot->getDynamicBiasDerivative();
            auto dgdq_test = robot->getGravityDerivative();
            auto t2 = high_resolution_clock::now();
            duration<double, std::milli> ms_double = t2 - t1;
            ed_avg += (1./n_samples) * ms_double.count();
            
            t1 = high_resolution_clock::now();
            VectorXdual q_rand_dual = q_rand.cast<dual>();
            VectorXdual dq_rand_dual = dq_rand.cast<dual>();
            // std::vector<MatrixXd> dMdq_truth = AutoDiffRigidBodyDynamics::massMatrixDerivative(ad_robot, q_rand_dual);
            auto dbdq_truth = AutoDiffRigidBodyDynamics::jacDynamicBiasWrtQ(ad_robot, q_rand_dual, dq_rand_dual);
            auto dbddq_truth = AutoDiffRigidBodyDynamics::jacDynamicBiasWrtDq(ad_robot, q_rand_dual, dq_rand_dual);
            auto dgdq_truth = AutoDiffRigidBodyDynamics::jacGravityBias(ad_robot, q_rand_dual);
            t2 = high_resolution_clock::now();
            ms_double = t2 - t1;
            ad_avg += (1./n_samples) * ms_double.count();

            double error = (dbdq_test - dbdq_truth).norm() + (dbddq_test - dbddq_truth).norm() + (dgdq_truth - dgdq_truth).norm();
            if (error > 1e-8) {
                std::cout << i << "\n";
                throw runtime_error("dynamic bias derivative doesn't match");
            }
        }

        std::cout << "ed avg: " << ed_avg << "\n";
        std::cout << "ad avg: " << ad_avg << "\n";
        std::cout << "dMdq test passed\n";
    }

}