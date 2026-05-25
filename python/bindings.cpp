#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "SaiModel.h"

namespace py = pybind11;

namespace {

Eigen::Matrix4d affineToMatrix(const Eigen::Affine3d& T) { return T.matrix(); }

Eigen::Affine3d matrixToAffine(const Eigen::Matrix4d& T) {
	if (T.rows() != 4 || T.cols() != 4) {
		throw std::invalid_argument("Expected a 4x4 homogeneous transform matrix");
	}
	return Eigen::Affine3d(T);
}

Eigen::Vector4d quatToXyzw(const Eigen::Quaterniond& q) {
	return Eigen::Vector4d(q.x(), q.y(), q.z(), q.w());
}

Eigen::Quaterniond quatFromXyzw(const Eigen::Vector4d& q) {
	return Eigen::Quaterniond(q(3), q(0), q(1), q(2));
}

}  // namespace

PYBIND11_MODULE(sai_model_py, m) {
	m.doc() = "Python bindings for the SaiModel C++ library";

	py::enum_<SaiModel::ContactType>(m, "ContactType")
		.value("PointContact", SaiModel::ContactType::PointContact)
		.value("SurfaceContact", SaiModel::ContactType::SurfaceContact)
		.export_values();

	py::class_<SaiModel::JointLimit>(m, "JointLimit")
		.def_readonly("joint_name", &SaiModel::JointLimit::joint_name)
		.def_readonly("joint_index", &SaiModel::JointLimit::joint_index)
		.def_readonly("position_lower", &SaiModel::JointLimit::position_lower)
		.def_readonly("position_upper", &SaiModel::JointLimit::position_upper)
		.def_readonly("velocity", &SaiModel::JointLimit::velocity)
		.def_readonly("effort", &SaiModel::JointLimit::effort);

	py::class_<SaiModel::SphericalJointDescription>(m, "SphericalJointDescription")
		.def_readonly("joint_name", &SaiModel::SphericalJointDescription::joint_name)
		.def_readonly("parent_link_name",
					  &SaiModel::SphericalJointDescription::parent_link_name)
		.def_readonly("child_link_name",
					  &SaiModel::SphericalJointDescription::child_link_name)
		.def_readonly("index", &SaiModel::SphericalJointDescription::index)
		.def_readonly("w_index", &SaiModel::SphericalJointDescription::w_index);

	py::class_<SaiModel::LinkMassParams>(m, "LinkMassParams")
		.def_readonly("mass", &SaiModel::LinkMassParams::mass)
		.def_readonly("com_pos", &SaiModel::LinkMassParams::com_pos)
		.def_readonly("inertia", &SaiModel::LinkMassParams::inertia)
		.def_readonly("link_name", &SaiModel::LinkMassParams::link_name);

	py::class_<SaiModel::OpSpaceMatrices>(m, "OpSpaceMatrices")
		.def_readonly("J", &SaiModel::OpSpaceMatrices::J)
		.def_readonly("Lambda", &SaiModel::OpSpaceMatrices::Lambda)
		.def_readonly("Jbar", &SaiModel::OpSpaceMatrices::Jbar)
		.def_readonly("N", &SaiModel::OpSpaceMatrices::N);

	py::class_<SaiModel::GraspMatrixData>(m, "GraspMatrixData")
		.def_readonly("G", &SaiModel::GraspMatrixData::G)
		.def_readonly("G_inv", &SaiModel::GraspMatrixData::G_inv)
		.def_readonly("R", &SaiModel::GraspMatrixData::R)
		.def_readonly("resultant_point", &SaiModel::GraspMatrixData::resultant_point);

	py::class_<SaiModel::ContactModel>(m, "ContactModel")
		.def(py::init<const std::string&, const Eigen::Vector3d&,
					  const Eigen::Matrix3d&, const SaiModel::ContactType&>(),
			 py::arg("link_name"), py::arg("pos"), py::arg("orientation"),
			 py::arg("contact_type"))
		.def_readwrite("contact_link_name", &SaiModel::ContactModel::contact_link_name)
		.def_readwrite("contact_position", &SaiModel::ContactModel::contact_position)
		.def_readwrite("contact_orientation",
					   &SaiModel::ContactModel::contact_orientation)
		.def_readwrite("contact_type", &SaiModel::ContactModel::contact_type);

	py::class_<SaiModel::SaiModel>(m, "SaiModel")
		.def(py::init<const std::string, bool>(), py::arg("urdf_path"),
			 py::arg("verbose") = false)
		.def_property("q", &SaiModel::SaiModel::q, &SaiModel::SaiModel::setQ)
		.def_property("dq", &SaiModel::SaiModel::dq, &SaiModel::SaiModel::setDq)
		.def_property("ddq", &SaiModel::SaiModel::ddq, &SaiModel::SaiModel::setDdq)
		.def_property_readonly("M", &SaiModel::SaiModel::M)
		.def_property_readonly("M_inv", &SaiModel::SaiModel::MInv)
		.def_property_readonly("dof", &SaiModel::SaiModel::dof)
		.def_property_readonly("q_size", &SaiModel::SaiModel::qSize)
		.def_property(
			"T_robot_base",
			[](const SaiModel::SaiModel& model) {
				return affineToMatrix(model.TRobotBase());
			},
			[](SaiModel::SaiModel& model, const Eigen::Matrix4d& T) {
				model.setTRobotBase(matrixToAffine(T));
			})
		.def_property("world_gravity", &SaiModel::SaiModel::worldGravity,
					  &SaiModel::SaiModel::setWorldGravity)
		.def("update_kinematics", &SaiModel::SaiModel::updateKinematics)
		.def("update_model", py::overload_cast<>(&SaiModel::SaiModel::updateModel))
		.def("update_model_with_mass_matrix",
			 py::overload_cast<const Eigen::MatrixXd&>(
				 &SaiModel::SaiModel::updateModel),
			 py::arg("M"))
		.def("joint_names", &SaiModel::SaiModel::jointNames)
		.def("is_link_in_robot", &SaiModel::SaiModel::isLinkInRobot,
			 py::arg("link_name"))
		.def("joint_limits", &SaiModel::SaiModel::jointLimits)
		.def("joint_limits_position_lower",
			 &SaiModel::SaiModel::jointLimitsPositionLower)
		.def("joint_limits_position_upper",
			 &SaiModel::SaiModel::jointLimitsPositionUpper)
		.def("spherical_joints", &SaiModel::SaiModel::sphericalJoints)
		.def("spherical_quat",
			 [](const SaiModel::SaiModel& model, const std::string& joint_name) {
				 return quatToXyzw(model.sphericalQuat(joint_name));
			 },
			 py::arg("joint_name"))
		.def("set_spherical_quat",
			 [](SaiModel::SaiModel& model, const std::string& joint_name,
				const Eigen::Vector4d& q_xyzw) {
				 model.setSphericalQuat(joint_name, quatFromXyzw(q_xyzw));
			 },
			 py::arg("joint_name"), py::arg("q_xyzw"))
		.def("joint_index", &SaiModel::SaiModel::jointIndex, py::arg("joint_name"))
		.def("spherical_joint_index_w", &SaiModel::SaiModel::sphericalJointIndexW,
			 py::arg("joint_name"))
		.def("joint_name", &SaiModel::SaiModel::jointName, py::arg("joint_id"))
		.def("child_link_name", &SaiModel::SaiModel::childLinkName,
			 py::arg("joint_name"))
		.def("parent_link_name", &SaiModel::SaiModel::parentLinkName,
			 py::arg("joint_name"))
		.def("joint_gravity_vector", &SaiModel::SaiModel::jointGravityVector)
		.def("coriolis_force", &SaiModel::SaiModel::coriolisForce)
		.def("coriolis_plus_gravity", &SaiModel::SaiModel::coriolisPlusGravity)
		.def("factorized_christoffel_matrix",
			 &SaiModel::SaiModel::factorizedChristoffelMatrix)
		.def("J", &SaiModel::SaiModel::J, py::arg("link_name"),
			 py::arg("pos_in_link") = Eigen::Vector3d::Zero())
		.def("J_world_frame", &SaiModel::SaiModel::JWorldFrame, py::arg("link_name"),
			 py::arg("pos_in_link") = Eigen::Vector3d::Zero())
		.def("J_local_frame", &SaiModel::SaiModel::JLocalFrame, py::arg("link_name"),
			 py::arg("pos_in_link") = Eigen::Vector3d::Zero(),
			 py::arg("rot_in_link") = Eigen::Matrix3d::Identity())
		.def("Jv", &SaiModel::SaiModel::Jv, py::arg("link_name"),
			 py::arg("pos_in_link") = Eigen::Vector3d::Zero())
		.def("Jv_world_frame", &SaiModel::SaiModel::JvWorldFrame,
			 py::arg("link_name"),
			 py::arg("pos_in_link") = Eigen::Vector3d::Zero())
		.def("Jv_local_frame", &SaiModel::SaiModel::JvLocalFrame,
			 py::arg("link_name"),
			 py::arg("pos_in_link") = Eigen::Vector3d::Zero(),
			 py::arg("rot_in_link") = Eigen::Matrix3d::Identity())
		.def("Jw", &SaiModel::SaiModel::Jw, py::arg("link_name"))
		.def("Jw_world_frame", &SaiModel::SaiModel::JwWorldFrame,
			 py::arg("link_name"))
		.def("Jw_local_frame", &SaiModel::SaiModel::JwLocalFrame,
			 py::arg("link_name"),
			 py::arg("rot_in_link") = Eigen::Matrix3d::Identity())
		.def("compute_inverse_kinematics",
			 py::overload_cast<const std::vector<std::string>&,
							   const std::vector<Eigen::Vector3d>&,
							   const std::vector<Eigen::Vector3d>&>(
				 &SaiModel::SaiModel::computeInverseKinematics),
			 py::arg("link_names"), py::arg("pos_in_links"),
			 py::arg("desired_pos_world_frame"))
		.def("compute_inverse_kinematics_frames",
			 [](SaiModel::SaiModel& model,
				const std::vector<std::string>& link_names,
				const std::vector<Eigen::Matrix4d>& frames_in_links,
				const std::vector<Eigen::Matrix4d>& desired_frames_world) {
				 std::vector<Eigen::Affine3d> frames_in_links_affine;
				 std::vector<Eigen::Affine3d> desired_frames_world_affine;
				 frames_in_links_affine.reserve(frames_in_links.size());
				 desired_frames_world_affine.reserve(desired_frames_world.size());
				 for (const auto& T : frames_in_links) {
					 frames_in_links_affine.push_back(matrixToAffine(T));
				 }
				 for (const auto& T : desired_frames_world) {
					 desired_frames_world_affine.push_back(matrixToAffine(T));
				 }
				 return model.computeInverseKinematics(
					 link_names, frames_in_links_affine, desired_frames_world_affine);
			 },
			 py::arg("link_names"), py::arg("frames_in_links"),
			 py::arg("desired_frames_world"))
		.def("transform",
			 [](const SaiModel::SaiModel& model, const std::string& link_name,
				const Eigen::Vector3d& pos_in_link,
				const Eigen::Matrix3d& rot_in_link) {
				 return affineToMatrix(
					 model.transform(link_name, pos_in_link, rot_in_link));
			 },
			 py::arg("link_name"), py::arg("pos_in_link") = Eigen::Vector3d::Zero(),
			 py::arg("rot_in_link") = Eigen::Matrix3d::Identity())
		.def("transform_in_world",
			 [](const SaiModel::SaiModel& model, const std::string& link_name,
				const Eigen::Vector3d& pos_in_link,
				const Eigen::Matrix3d& rot_in_link) {
				 return affineToMatrix(
					 model.transformInWorld(link_name, pos_in_link, rot_in_link));
			 },
			 py::arg("link_name"), py::arg("pos_in_link") = Eigen::Vector3d::Zero(),
			 py::arg("rot_in_link") = Eigen::Matrix3d::Identity())
		.def("velocity6d", &SaiModel::SaiModel::velocity6d, py::arg("link_name"),
			 py::arg("pos_in_link") = Eigen::Vector3d::Zero())
		.def("velocity6d_in_world", &SaiModel::SaiModel::velocity6dInWorld,
			 py::arg("link_name"), py::arg("pos_in_link") = Eigen::Vector3d::Zero())
		.def("acceleration6d", &SaiModel::SaiModel::acceleration6d,
			 py::arg("link_name"), py::arg("pos_in_link") = Eigen::Vector3d::Zero())
		.def("acceleration6d_in_world", &SaiModel::SaiModel::acceleration6dInWorld,
			 py::arg("link_name"), py::arg("pos_in_link") = Eigen::Vector3d::Zero())
		.def("position", &SaiModel::SaiModel::position, py::arg("link_name"),
			 py::arg("pos_in_link") = Eigen::Vector3d::Zero())
		.def("position_in_world", &SaiModel::SaiModel::positionInWorld,
			 py::arg("link_name"), py::arg("pos_in_link") = Eigen::Vector3d::Zero())
		.def("linear_velocity", &SaiModel::SaiModel::linearVelocity,
			 py::arg("link_name"), py::arg("pos_in_link") = Eigen::Vector3d::Zero())
		.def("linear_velocity_in_world", &SaiModel::SaiModel::linearVelocityInWorld,
			 py::arg("link_name"), py::arg("pos_in_link") = Eigen::Vector3d::Zero())
		.def("linear_acceleration", &SaiModel::SaiModel::linearAcceleration,
			 py::arg("link_name"), py::arg("pos_in_link") = Eigen::Vector3d::Zero())
		.def("linear_acceleration_in_world",
			 &SaiModel::SaiModel::linearAccelerationInWorld, py::arg("link_name"),
			 py::arg("pos_in_link") = Eigen::Vector3d::Zero())
		.def("rotation", &SaiModel::SaiModel::rotation, py::arg("link_name"),
			 py::arg("rot_in_link") = Eigen::Matrix3d::Identity())
		.def("rotation_in_world", &SaiModel::SaiModel::rotationInWorld,
			 py::arg("link_name"),
			 py::arg("rot_in_link") = Eigen::Matrix3d::Identity())
		.def("angular_velocity", &SaiModel::SaiModel::angularVelocity,
			 py::arg("link_name"))
		.def("angular_velocity_in_world", &SaiModel::SaiModel::angularVelocityInWorld,
			 py::arg("link_name"))
		.def("angular_acceleration", &SaiModel::SaiModel::angularAcceleration,
			 py::arg("link_name"))
		.def("angular_acceleration_in_world",
			 &SaiModel::SaiModel::angularAccelerationInWorld, py::arg("link_name"))
		.def("get_link_mass_params", &SaiModel::SaiModel::getLinkMassParams,
			 py::arg("link_name"))
		.def("com_position", &SaiModel::SaiModel::comPosition)
		.def("com_jacobian", &SaiModel::SaiModel::comJacobian)
		.def("task_inertia_matrix", &SaiModel::SaiModel::taskInertiaMatrix,
			 py::arg("task_jacobian"))
		.def("task_inertia_matrix_with_pseudo_inv",
			 &SaiModel::SaiModel::taskInertiaMatrixWithPseudoInv,
			 py::arg("task_jacobian"))
		.def("dyn_consistent_inverse_jacobian",
			 &SaiModel::SaiModel::dynConsistentInverseJacobian,
			 py::arg("task_jacobian"))
		.def("nullspace_matrix", &SaiModel::SaiModel::nullspaceMatrix,
			 py::arg("task_jacobian"))
		.def("operational_space_matrices",
			 &SaiModel::SaiModel::operationalSpaceMatrices,
			 py::arg("task_jacobian"))
		.def("add_environmental_contact", &SaiModel::SaiModel::addEnvironmentalContact,
			 py::arg("link"), py::arg("pos_in_link") = Eigen::Vector3d::Zero(),
			 py::arg("orientation") = Eigen::Matrix3d::Identity(),
			 py::arg("contact_type") = SaiModel::ContactType::SurfaceContact)
		.def("add_manipulation_contact", &SaiModel::SaiModel::addManipulationContact,
			 py::arg("link"), py::arg("pos_in_link") = Eigen::Vector3d::Zero(),
			 py::arg("orientation") = Eigen::Matrix3d::Identity(),
			 py::arg("contact_type") = SaiModel::ContactType::SurfaceContact)
		.def("update_environmental_contact",
			 &SaiModel::SaiModel::updateEnvironmentalContact, py::arg("link"),
			 py::arg("pos_in_link"), py::arg("orientation"),
			 py::arg("contact_type"))
		.def("update_manipulation_contact",
			 &SaiModel::SaiModel::updateManipulationContact, py::arg("link"),
			 py::arg("pos_in_link"), py::arg("orientation"),
			 py::arg("contact_type"))
		.def("delete_environmental_contact",
			 &SaiModel::SaiModel::deleteEnvironmentalContact, py::arg("link_name"))
		.def("delete_manipulation_contact",
			 &SaiModel::SaiModel::deleteManipulationContact, py::arg("link_name"))
		.def("manipulation_grasp_matrix", &SaiModel::SaiModel::manipulationGraspMatrix,
			 py::arg("center_point"), py::arg("resultant_in_world_frame") = false,
			 py::arg("contact_forces_in_local_frames") = false)
		.def("manipulation_grasp_matrix_at_geometric_center",
			 &SaiModel::SaiModel::manipulationGraspMatrixAtGeometricCenter,
			 py::arg("resultant_in_world_frame") = false,
			 py::arg("contact_forces_in_local_frames") = false)
		.def("environmental_grasp_matrix",
			 &SaiModel::SaiModel::environmentalGraspMatrix,
			 py::arg("center_point"), py::arg("resultant_in_world_frame") = false,
			 py::arg("contact_forces_in_local_frames") = false)
		.def("environmental_grasp_matrix_at_geometric_center",
			 &SaiModel::SaiModel::environmentalGraspMatrixAtGeometricCenter,
			 py::arg("resultant_in_world_frame") = false,
			 py::arg("contact_forces_in_local_frames") = false)
		.def("display_joints", &SaiModel::SaiModel::displayJoints)
		.def("display_links", &SaiModel::SaiModel::displayLinks)
		.def("link_dependency", &SaiModel::SaiModel::linkDependency,
			 py::arg("link_name"), py::arg("update") = false)
		.def("link_dependency_vector", &SaiModel::SaiModel::linkDependencyVector,
			 py::arg("link_name"), py::arg("update") = false)
		.def("jdot_qdot", &SaiModel::SaiModel::jDotQDot, py::arg("link_name"),
			 py::arg("pos_in_link") = Eigen::Vector3d::Zero(),
			 py::arg("update_kinematics") = false)
		.def("get_centroidal_inertia_matrix",
			 &SaiModel::SaiModel::getCentroidalInertiaMatrix)
		.def("get_centroidal_momentum_matrix",
			 &SaiModel::SaiModel::getCentroidalMomentumMatrix)
		.def("get_point_inertia_matrix", &SaiModel::SaiModel::getPointInertiaMatrix,
			 py::arg("link_name"), py::arg("pos_in_link"))
		.def("get_jacobian_derivative", &SaiModel::SaiModel::getJacobianDerivative,
			 py::arg("link_name"), py::arg("pos_in_link"), py::arg("update") = false)
		.def("get_mass_matrix_derivative",
			 &SaiModel::SaiModel::getMassMatrixDerivative, py::arg("update") = false)
		.def("get_dynamic_bias_derivative_wrt_q",
			 &SaiModel::SaiModel::getDynamicBiasDerivativeWrtQ,
			 py::arg("gravity_opt") = true)
		.def("get_dynamic_bias_derivative_wrt_dq",
			 &SaiModel::SaiModel::getDynamicBiasDerivativeWrtDq,
			 py::arg("gravity_opt") = true)
		.def("get_gravity_derivative", &SaiModel::SaiModel::getGravityDerivative)
		.def("add_load",
			 [](SaiModel::SaiModel& model, const std::string& link_name, double mass,
				const Eigen::Vector3d& com_pos, const Eigen::Matrix3d& inertia,
				const Eigen::Matrix4d& link_transform,
				const std::string& body_name) {
				 model.addLoad(link_name, mass, com_pos, inertia,
							   matrixToAffine(link_transform), body_name);
			 },
			 py::arg("link_name"), py::arg("mass"), py::arg("com_pos"),
			 py::arg("inertia"),
			 py::arg("link_transform") = Eigen::Matrix4d::Identity(),
			 py::arg("body_name") = "")
		.def("remove_load", &SaiModel::SaiModel::removeLoad,
			 py::arg("body_name"))
			.def("com_acceleration", &SaiModel::SaiModel::comAcceleration)
			.def("add_muscle_system", &SaiModel::SaiModel::addMuscleSystem,
				 py::arg("muscle_xml"), py::arg("name"))
			.def("get_num_muscles", &SaiModel::SaiModel::getNumMuscles)
			.def("compute_muscle_capacity_matrix",
				 &SaiModel::SaiModel::computeMuscleCapacityMatrix)
			.def("compute_muscle_jacobian",
				 &SaiModel::SaiModel::computeMuscleJacobian,
				 py::arg("floating") = true)
			.def("compute_muscle_jacobian_derivative",
				 &SaiModel::SaiModel::computeMuscleJacobianDerivative,
				 py::arg("floating") = true)
			.def("compute_muscle_jacobian_inverse",
				 &SaiModel::SaiModel::computeMuscleJacobianInverse,
				 py::arg("W"))
			.def("compute_muscle_jacobian_inverse_derivative",
				 &SaiModel::SaiModel::computeMuscleJacobianInverseDerivative,
				 py::arg("W"));
	
		m.def("compute_pseudo_inverse", &SaiModel::computePseudoInverse,
			  py::arg("matrix"), py::arg("svd_epsilon") = 1e-6);
		m.def("matrix_range_basis", &SaiModel::matrixRangeBasis, py::arg("matrix"),
			  py::arg("svd_epsilon") = 1e-6);
	m.def("orientation_error", py::overload_cast<const Eigen::Matrix3d&,
												 const Eigen::Matrix3d&>(
								 &SaiModel::orientationError),
		  py::arg("desired_orientation"), py::arg("current_orientation"));
	m.def("orientation_error_quat",
		  py::overload_cast<const Eigen::Quaterniond&,
							const Eigen::Quaterniond&>(&SaiModel::orientationError),
		  py::arg("desired_orientation"), py::arg("current_orientation"));
	m.def("cross_product_operator", &SaiModel::crossProductOperator, py::arg("v"));
	m.def(
		"grasp_matrix_at_geometric_center",
		&SaiModel::graspMatrixAtGeometricCenter, py::arg("contact_locations"),
		py::arg("contact_types"));
}
