#include "rbdl/rbdl.h"

#include "ModelAD.h"
#include "ConstraintsAD.h"

struct CartPendulum {
    CartPendulum () {
        using namespace RigidBodyDynamics;
        using namespace RigidBodyDynamics::Math;

        ClearLogOutput();

        double cart_w = 0.5;
        double cart_d = 0.2;
        double cart_h = 0.2;
        double cart_m = 100.0;
        Vector3d cart_com (0., 0., 0.);

        double pend_l = 0.5;
        double pend_r = 0.1;
        double pend_m = 1.0;
        Vector3d pend_com (0., 0., pend_l);

        body_cart = Body (cart_m, cart_com, Matrix3d (
                    1. / 12. * cart_m * (cart_h*cart_h + cart_d * cart_d), 0., 0.,
                    0., 1. / 12. * cart_m * (cart_w*cart_w + cart_h * cart_h), 0.,
                    0., 0., 1. / 12. * cart_m * (cart_w*cart_w + cart_d * cart_d)
                    ));
        body_pendulum = Body (pend_m, pend_com, Matrix3d(
                    2. / 5. * pend_m * pend_r, 0., 0.,
                    0., 2. / 5. * pend_m * pend_r, 0.,
                    0., 0., 2. / 5. * pend_m * pend_r
                ));

        joint_cart = Joint (SpatialVector (0., 0., 0., 1., 0., 0.));
        joint_pendulum = Joint (SpatialVector (0., 1., 0., 0., 0., 0.));

        id_cart = model.AddBody (0, Xtrans(Vector3d (0., 0., 0.)), joint_cart, body_cart, "cart");
        id_pendulum = model.AddBody (id_cart, Xtrans(Vector3d (0., 0., 0.)), joint_pendulum, body_pendulum, "pendulum");

        q = VectorNd::Constant ((size_t) model.dof_count, 0.);
        qdot = VectorNd::Constant ((size_t) model.dof_count, 0.);
        qddot = VectorNd::Constant ((size_t) model.dof_count, 0.);
        tau = VectorNd::Constant ((size_t) model.dof_count, 0.);

        body_point = RigidBodyDynamics::Math::Vector3d (0., 0., pend_l);

        ClearLogOutput();
    }

    RigidBodyDynamics::Model model;

    unsigned int id_cart, id_pendulum;
    RigidBodyDynamics::Body body_cart, body_pendulum;
    RigidBodyDynamics::Joint joint_cart, joint_pendulum;

    RigidBodyDynamics::Math::VectorNd q;
    RigidBodyDynamics::Math::VectorNd qdot;
    RigidBodyDynamics::Math::VectorNd qddot;
    RigidBodyDynamics::Math::VectorNd tau;

    RigidBodyDynamics::Math::Vector3d body_point;
};

struct CartPendulumContact {
    CartPendulumContact () {
        using namespace RigidBodyDynamics;
        using namespace RigidBodyDynamics::Math;

        ClearLogOutput();

        double cart_w = 0.5;
        double cart_d = 0.2;
        double cart_h = 0.2;
        double cart_m = 100.0;
        Vector3d cart_com (0., 0., 0.);

        double pend_l = 0.5;
        double pend_r = 0.1;
        double pend_m = 1.0;
        Vector3d pend_com (0., 0., pend_l);

        body_cart = Body (cart_m, cart_com, Matrix3d (
                    1. / 12. * cart_m * (cart_h*cart_h + cart_d * cart_d), 0., 0.,
                    0., 1. / 12. * cart_m * (cart_w*cart_w + cart_h * cart_h), 0.,
                    0., 0., 1. / 12. * cart_m * (cart_w*cart_w + cart_d * cart_d)
                    ));
        body_pendulum = Body (pend_m, pend_com, Matrix3d(
                    2. / 5. * pend_m * pend_r, 0., 0.,
                    0., 2. / 5. * pend_m * pend_r, 0.,
                    0., 0., 2. / 5. * pend_m * pend_r
                ));

        // NOTE we use a free-flyer like joint for the cart and restrain it to
        //      the x-axis using a constraint set!
        joint_cart = Joint (
            SpatialVector (0., 0., 0., 1., 0., 0.),
            SpatialVector (0., 0., 0., 0., 0., 1.)//,
            // TODO add with 6DOF contact in CS
            // SpatialVector (0., 1., 0., 0., 0., 0.)
        );
        joint_pendulum = Joint (SpatialVector (0., 1., 0., 0., 0., 0.));

        id_cart = model.AddBody (0, Xtrans(Vector3d (0., 0., 0.)), joint_cart, body_cart, "cart");
        id_pendulum = model.AddBody (id_cart, Xtrans(Vector3d (0., 0., 0.)), joint_pendulum, body_pendulum, "pendulum");

        q = VectorNd::Constant ((size_t) model.dof_count, 0.);
        qdot = VectorNd::Constant ((size_t) model.dof_count, 0.);
        qddot = VectorNd::Constant ((size_t) model.dof_count, 0.);
        tau = VectorNd::Constant ((size_t) model.dof_count, 0.);

        body_point = RigidBodyDynamics::Math::Vector3d (0., 0., pend_l);

        ClearLogOutput();
    }

    RigidBodyDynamics::Model model;

    unsigned int id_cart, id_pendulum;
    RigidBodyDynamics::Body body_cart, body_pendulum;
    RigidBodyDynamics::Joint joint_cart, joint_pendulum;

    RigidBodyDynamics::Math::VectorNd q;
    RigidBodyDynamics::Math::VectorNd qdot;
    RigidBodyDynamics::Math::VectorNd qddot;
    RigidBodyDynamics::Math::VectorNd tau;

    RigidBodyDynamics::Math::Vector3d body_point;
};

struct MultiPendulumWithBranches {

  unsigned int nbranches;

  RigidBodyDynamics::Model model;

  // create multi-pendulum with given number of degrees
  MultiPendulumWithBranches()
    : nbranches(0)
  {};

  void create_model(
    const unsigned int& dof_count,
    const unsigned int& branch_count
  )
  {
    using namespace RigidBodyDynamics;
    using namespace RigidBodyDynamics::Math;

    // does number of branches exceed number of DoFs?
    if (branch_count > dof_count)
    {
      std::cerr << "no of branches exceeds number of DoFs" << std::endl;
      abort();
    }

    // coordinate axis convention
    //  x: points to the right
    //  y: points up
    //  z: out of the page
    model.gravity = Vector3d(0.,-9.81,0.);

    Joint joint (JointTypeRevoluteZ);
    Body body (
      1.0,
      Vector3d (0., 0., 0.),
      Matrix3d(
        1., 0., 0.,
        0., 1., 0.,
        0., 0., 1.
      )
    );
    SpatialTransform X_T = Xtrans (Vector3d (0., 1., 0.));

    unsigned int branches = 0;
    for (unsigned int i = 0; i < dof_count; i++) {
      if ((branch_count > 0) && (dof_count/branch_count)) {
        if ((i%(dof_count/branch_count)) && (branches < branch_count)) {
          model.AddBody (0, X_T, joint, body);
          branches++;
        } else {
          model.AppendBody (X_T, joint, body);
        }
      } else {
        model.AppendBody (X_T, joint, body);
      }
    }

    // set number of branches
    nbranches = branches;
  }
};

struct Arm2DofX {
    Arm2DofX() {
        using namespace RigidBodyDynamics;
        using namespace RigidBodyDynamics::Math;

        ClearLogOutput();

        double arm_l = 0.5;
        double arm_r = 0.1;
        double arm_m = 3.5;
        Matrix3d arm_I = .5 * arm_m * Matrix3d(
                    .5 * (arm_r * arm_r + arm_l * arm_l / 3.), 0., 0.,
                    arm_r * arm_r, 0., 0.,
                    .5 * (arm_r * arm_r + arm_l * arm_l / 3.), 0., 0.);
        Vector3d arm_com (0., .5 * arm_l, 0.);

        body_distal    = Body(arm_m, arm_com, arm_I);
        body_proximal  = Body(arm_m, arm_com, arm_I);

        joint_distal   = Joint(SpatialVector(1., 0., 0., 0., 0., 0.));
        joint_proximal = Joint(SpatialVector(1., 0., 0., 0., 0., 0.));

        id_distal      = model.AddBody(0, Xtrans(Vector3d(0., 0., 0.)),
                                       joint_distal, body_distal, "distal");

        id_proximal    = model.AddBody(id_distal, Xtrans(Vector3d(0., 0., 0.)),
                                       joint_proximal, body_proximal,
                                       "proximal");

        q     = VectorNd::Constant((size_t) model.dof_count, 0.);
        qdot  = VectorNd::Constant((size_t) model.dof_count, 0.);
        qddot = VectorNd::Constant((size_t) model.dof_count, 0.);
        tau   = VectorNd::Constant((size_t) model.dof_count, 0.);

        body_point = RigidBodyDynamics::Math::Vector3d (0., arm_l, 0.);

        ClearLogOutput();
    }

    RigidBodyDynamics::Model model;

    unsigned int id_distal, id_proximal;
    RigidBodyDynamics::Body  body_distal, body_proximal;
    RigidBodyDynamics::Joint joint_distal, joint_proximal;

    RigidBodyDynamics::Math::VectorNd q;
    RigidBodyDynamics::Math::VectorNd qdot;
    RigidBodyDynamics::Math::VectorNd qddot;
    RigidBodyDynamics::Math::VectorNd tau;

    RigidBodyDynamics::Math::Vector3d body_point;
};

struct Arm2DofZ {
    Arm2DofZ() {
        using namespace RigidBodyDynamics;
        using namespace RigidBodyDynamics::Math;

        ClearLogOutput();

        double arm_l = 0.5;
        double arm_r = 0.1;
        double arm_m = 3.5;
        Matrix3d arm_I = .5 * arm_m * Matrix3d(
                    .5 * (arm_r * arm_r + arm_l * arm_l / 3.), 0., 0.,
                    0., arm_r * arm_r, 0.,
                    0., 0., .5 * (arm_r * arm_r + arm_l * arm_l / 3.));
        Vector3d arm_com (0., .5 * arm_l, 0.);

        body_distal    = Body(arm_m, arm_com, arm_I);
        body_proximal  = Body(arm_m, arm_com, arm_I);

        joint_distal   = Joint(SpatialVector(0., 0., 1., 0., 0., 0.));
        joint_proximal = Joint(SpatialVector(0., 0., 1., 0., 0., 0.));

        id_distal      = model.AddBody(0, Xtrans(Vector3d(0., 0., 0.)),
                                       joint_distal, body_distal, "distal");

        id_proximal    = model.AddBody(id_distal, Xtrans(Vector3d(0., 0., 0.)),
                                       joint_proximal, body_proximal,
                                       "proximal");

        q     = VectorNd::Constant((size_t) model.dof_count, 0.);
        qdot  = VectorNd::Constant((size_t) model.dof_count, 0.);
        qddot = VectorNd::Constant((size_t) model.dof_count, 0.);
        tau   = VectorNd::Constant((size_t) model.dof_count, 0.);

        body_point = RigidBodyDynamics::Math::Vector3d (0., arm_l, 0.);

        ClearLogOutput();
    }

    RigidBodyDynamics::Model model;

    unsigned int id_distal, id_proximal;
    RigidBodyDynamics::Body  body_distal, body_proximal;
    RigidBodyDynamics::Joint joint_distal, joint_proximal;

    RigidBodyDynamics::Math::VectorNd q;
    RigidBodyDynamics::Math::VectorNd qdot;
    RigidBodyDynamics::Math::VectorNd qddot;
    RigidBodyDynamics::Math::VectorNd tau;

    RigidBodyDynamics::Math::Vector3d body_point;
};

struct Arm3DofXZYp {
    Arm3DofXZYp() {
        using namespace RigidBodyDynamics;
        using namespace RigidBodyDynamics::Math;

        ClearLogOutput();

        double arm_l = 0.35;
        double arm_r = 0.15;
        double arm_m = 4.0;
        Matrix3d arm_I = .5 * arm_m * Matrix3d(
                    .5 * (arm_r * arm_r + arm_l * arm_l / 3.), 0., 0.,
                    0., arm_r * arm_r, 0.,
                    0., 0., .5 * (arm_r * arm_r + arm_l * arm_l / 3.));
        Vector3d arm_com (0., .5 * arm_l, 0.);

        body_distal    = Body(arm_m, arm_com, arm_I);
        body_proximal  = Body(arm_m, arm_com, arm_I);
        body_slider  = Body(arm_m, arm_com, arm_I);

        joint_distal   = Joint(SpatialVector(1., 0., 0., 0., 0., 0.));
        joint_proximal = Joint(SpatialVector(0., 0., 1., 0., 0., 0.));
        joint_slider   = Joint(SpatialVector(0., 0., 0., 0., 1., 0.));

        id_distal      = model.AddBody(0, Xtrans(Vector3d(0., 0., 0.)),
                                       joint_distal, body_distal, "distal");

        id_proximal    = model.AddBody(id_distal, Xtrans(Vector3d(0., 0., 0.)),
                                       joint_proximal, body_proximal,
                                       "proximal");

        id_slider      = model.AddBody(id_proximal, Xtrans(Vector3d(0., 0., 0.)),
                                       joint_slider, body_slider,
                                       "slider");

        q     = VectorNd::Constant((size_t) model.dof_count, 0.);
        qdot  = VectorNd::Constant((size_t) model.dof_count, 0.);
        qddot = VectorNd::Constant((size_t) model.dof_count, 0.);
        tau   = VectorNd::Constant((size_t) model.dof_count, 0.);

        body_point = RigidBodyDynamics::Math::Vector3d (0., arm_l, 0.);

        ClearLogOutput();
    }

    RigidBodyDynamics::Model model;
    ADModel ad_model;

    unsigned int id_distal, id_proximal, id_slider;
    RigidBodyDynamics::Body  body_distal, body_proximal, body_slider;
    RigidBodyDynamics::Joint joint_distal, joint_proximal, joint_slider;

    RigidBodyDynamics::Math::VectorNd q;
    RigidBodyDynamics::Math::VectorNd qdot;
    RigidBodyDynamics::Math::VectorNd qddot;
    RigidBodyDynamics::Math::VectorNd tau;

    RigidBodyDynamics::Math::Vector3d body_point;
};

struct Arm3DofXZZp {
    Arm3DofXZZp() {
        using namespace RigidBodyDynamics;
        using namespace RigidBodyDynamics::Math;

        ClearLogOutput();

        double arm_l = 0.35;
        double arm_r = 0.15;
        double arm_m = 4.0;
        Matrix3d arm_I = .5 * arm_m * Matrix3d(
                    .5 * (arm_r * arm_r + arm_l * arm_l / 3.), 0., 0.,
                    arm_r * arm_r, 0., 0.,
                    .5 * (arm_r * arm_r + arm_l * arm_l / 3.), 0., 0.);
        Vector3d arm_com (0., .5 * arm_l, 0.);

        Matrix3d arm_slider_I = .5 * arm_m * Matrix3d(
                    .5 * (arm_r * arm_r + arm_l * arm_l / 3.), 0., 0.,
                    0., .5 * (arm_r * arm_r + arm_l * arm_l / 3.), 0.,
                    0., 0., arm_r * arm_r);
        Vector3d arm_slider_com(0., 0., .5 * arm_l);


        body_distal    = Body(arm_m, arm_com, arm_I);
        body_proximal  = Body(arm_m, arm_com, arm_I);
        body_slider  = Body(arm_m, arm_slider_com, arm_slider_I);

        joint_distal   = Joint(SpatialVector(1., 0., 0., 0., 0., 0.));
        joint_proximal = Joint(SpatialVector(0., 0., 1., 0., 0., 0.));
        joint_slider   = Joint(SpatialVector(0., 0., 0., 0., 0., 1.));

        id_distal      = model.AddBody(0, Xtrans(Vector3d(0., 0., 0.)),
                                       joint_distal, body_distal, "distal");

        id_proximal    = model.AddBody(id_distal, Xtrans(Vector3d(0., 0., 0.)),
                                       joint_proximal, body_proximal,
                                       "proximal");

        id_slider      = model.AddBody(id_proximal, Xtrans(Vector3d(0., 0., 0.)),
                                       joint_slider, body_slider,
                                       "slider");

        q     = VectorNd::Constant((size_t) model.dof_count, 0.);
        qdot  = VectorNd::Constant((size_t) model.dof_count, 0.);
        qddot = VectorNd::Constant((size_t) model.dof_count, 0.);
        tau   = VectorNd::Constant((size_t) model.dof_count, 0.);

        ad_model = ADModel(model);

        body_point = RigidBodyDynamics::Math::Vector3d (0., arm_l, 0.);

        ClearLogOutput();
    }

    RigidBodyDynamics::Model model;
    ADModel ad_model;

    unsigned int id_distal, id_proximal, id_slider;
    RigidBodyDynamics::Body  body_distal, body_proximal, body_slider;
    RigidBodyDynamics::Joint joint_distal, joint_proximal, joint_slider;

    RigidBodyDynamics::Math::VectorNd q;
    RigidBodyDynamics::Math::VectorNd qdot;
    RigidBodyDynamics::Math::VectorNd qddot;
    RigidBodyDynamics::Math::VectorNd tau;

    RigidBodyDynamics::Math::Vector3d body_point;
};

struct FixedBase6DoF {
    FixedBase6DoF () {
        using namespace RigidBodyDynamics;
        using namespace RigidBodyDynamics::Math;

        ClearLogOutput();

        model.gravity = Vector3d  (0., -9.81, 0.);

        /* 3 DoF (rot.) joint at base
         * 3 DoF (rot.) joint child origin
         *
         *          X Contact point (ref child)
         *          |
         *    Base  |
         *   / body |
         *  O-------*
         *           \
         *             Child body
         */

        // base body (3 DoF)
        base_rot_z = Body (
                0.,
                Vector3d (0., 0., 0.),
                Vector3d (0., 0., 0.)
                );
        joint_base_rot_z = Joint ( SpatialVector (0., 0., 1., 0., 0., 0.));
        base_rot_z_id = model.AddBody (0, Xtrans (Vector3d (0., 0., 0.)), joint_base_rot_z, base_rot_z);

        base_rot_y = Body (
                0.,
                Vector3d (0., 0., 0.),
                Vector3d (0., 0., 0.)
                );
        joint_base_rot_y = Joint ( SpatialVector (0., 1., 0., 0., 0., 0.));
        base_rot_y_id = model.AppendBody (Xtrans (Vector3d (0., 0., 0.)), joint_base_rot_y, base_rot_y);

        base_rot_x = Body (
                1.,
                Vector3d (0.5, 0., 0.),
                Vector3d (1., 1., 1.)
                );
        joint_base_rot_x = Joint ( SpatialVector (1., 0., 0., 0., 0., 0.));
        base_rot_x_id = model.AddBody (base_rot_y_id, Xtrans (Vector3d (0., 0., 0.)), joint_base_rot_x, base_rot_x);

        // child body (3 DoF)
        child_rot_z = Body (
                0.,
                Vector3d (0., 0., 0.),
                Vector3d (0., 0., 0.)
                );
        joint_child_rot_z = Joint ( SpatialVector (0., 0., 1., 0., 0., 0.));
        child_rot_z_id = model.AddBody (base_rot_x_id, Xtrans (Vector3d (1., 0., 0.)), joint_child_rot_z, child_rot_z);

        child_rot_y = Body (
                0.,
                Vector3d (0., 0., 0.),
                Vector3d (0., 0., 0.)
                );
        joint_child_rot_y = Joint ( SpatialVector (0., 1., 0., 0., 0., 0.));
        child_rot_y_id = model.AddBody (child_rot_z_id, Xtrans (Vector3d (0., 0., 0.)), joint_child_rot_y, child_rot_y);

        child_rot_x = Body (
                1.,
                Vector3d (0., 0.5, 0.),
                Vector3d (1., 1., 1.)
                );
        joint_child_rot_x = Joint ( SpatialVector (1., 0., 0., 0., 0., 0.));
        child_rot_x_id = model.AddBody (child_rot_y_id, Xtrans (Vector3d (0., 0., 0.)), joint_child_rot_x, child_rot_x);

        Q = VectorNd::Constant (model.mBodies.size() - 1, 0.);
        QDot = VectorNd::Constant (model.mBodies.size() - 1, 0.);
        QDDot = VectorNd::Constant (model.mBodies.size() - 1, 0.);
        Tau = VectorNd::Constant (model.mBodies.size() - 1, 0.);

        contact_body_id = child_rot_x_id;
        contact_point = Vector3d  (0.5, 0.5, 0.);
        contact_normal = Vector3d  (0., 1., 0.);

        ClearLogOutput();
        ad_model = ADModel(model);
    }

    ~FixedBase6DoF () {
    }

    RigidBodyDynamics::Model model;
    ADModel ad_model;

    unsigned int base_rot_z_id, base_rot_y_id, base_rot_x_id,
        child_rot_z_id, child_rot_y_id, child_rot_x_id,
        base_body_id;

    RigidBodyDynamics::Body base_rot_z, base_rot_y, base_rot_x,
        child_rot_z, child_rot_y, child_rot_x;

    RigidBodyDynamics::Joint joint_base_rot_z, joint_base_rot_y, joint_base_rot_x,
        joint_child_rot_z, joint_child_rot_y, joint_child_rot_x;

    RigidBodyDynamics::Math::VectorNd Q;
    RigidBodyDynamics::Math::VectorNd QDot;
    RigidBodyDynamics::Math::VectorNd QDDot;
    RigidBodyDynamics::Math::VectorNd Tau;

    unsigned int contact_body_id;
    RigidBodyDynamics::Math::Vector3d contact_point;
    RigidBodyDynamics::Math::Vector3d contact_normal;

    RigidBodyDynamics::ConstraintSet constraint_set;
};

struct FixedBase6DoF9DoF {
    FixedBase6DoF9DoF () {
        ClearLogOutput();

        model.gravity = RigidBodyDynamics::Math::Vector3d  (0., -9.81, 0.);

        /* 3 DoF (rot.) joint at base
         * 3 DoF (rot.) joint child origin
         *
         *          X Contact point (ref child)
         *          |
         *    Base  |
         *   / body |
         *  O-------*
         *           \
         *             Child body
         */

        // base body (3 DoF)
        base = RigidBodyDynamics::Body (
                1.,
                RigidBodyDynamics::Math::Vector3d (0.5, 0., 0.),
                RigidBodyDynamics::Math::Vector3d (1., 1., 1.)
                );
        joint_rotzyx = RigidBodyDynamics::Joint (
                RigidBodyDynamics::Math::SpatialVector (0., 0., 1., 0., 0., 0.),
                RigidBodyDynamics::Math::SpatialVector (0., 1., 0., 0., 0., 0.),
                RigidBodyDynamics::Math::SpatialVector (1., 0., 0., 0., 0., 0.)
                );
        base_id = model.AddBody (
                0,
                RigidBodyDynamics::Math::Xtrans (
                    RigidBodyDynamics::Math::Vector3d (0., 0., 0.)
                    ),
                joint_rotzyx, base
                );

        // child body 1 (3 DoF)
        child = RigidBodyDynamics::Body (
                1.,
                RigidBodyDynamics::Math::Vector3d (0., 0.5, 0.),
                RigidBodyDynamics::Math::Vector3d (1., 1., 1.)
                );
        child_id = model.AddBody (
                base_id,
                RigidBodyDynamics::Math::Xtrans (
                    RigidBodyDynamics::Math::Vector3d (0., 0., 0.)
                    ),
                joint_rotzyx, child
                );

        // child body (3 DoF)
        child_2 = RigidBodyDynamics::Body (
                1.,
                RigidBodyDynamics::Math::Vector3d (0., 0.5, 0.),
                RigidBodyDynamics::Math::Vector3d (1., 1., 1.)
                );
        child_2_id = model.AddBody (
                child_id,
                RigidBodyDynamics::Math::Xtrans (
                    RigidBodyDynamics::Math::Vector3d (0., 0., 0.)
                    ),
                joint_rotzyx,
                child_2
                );

        Q = RigidBodyDynamics::Math::VectorNd::Constant (model.mBodies.size() - 1, 0.);
        QDot = RigidBodyDynamics::Math::VectorNd::Constant (model.mBodies.size() - 1, 0.);
        QDDot = RigidBodyDynamics::Math::VectorNd::Constant (model.mBodies.size() - 1, 0.);
        Tau = RigidBodyDynamics::Math::VectorNd::Constant (model.mBodies.size() - 1, 0.);

        contact_body_id = child_id;
        contact_point = RigidBodyDynamics::Math::Vector3d  (0.5, 0.5, 0.);
        contact_normal =RigidBodyDynamics::Math::Vector3d  (0., 1., 0.);

        ClearLogOutput();
    }

    RigidBodyDynamics::Model model;
    ADModel ad_model;

    unsigned int base_id, child_id, child_2_id;

    RigidBodyDynamics::Body base, child, child_2;

    RigidBodyDynamics::Joint joint_rotzyx;

    RigidBodyDynamics::Math::VectorNd Q;
    RigidBodyDynamics::Math::VectorNd QDot;
    RigidBodyDynamics::Math::VectorNd QDDot;
    RigidBodyDynamics::Math::VectorNd Tau;

    unsigned int contact_body_id;
    RigidBodyDynamics::Math::Vector3d contact_point;
    RigidBodyDynamics::Math::Vector3d contact_normal;
    RigidBodyDynamics::ConstraintSet constraint_set;
};
