// The following function is a replacement for the function SparseSolveLx in src/rbdl_mathutils.cc . It can take Vectors and ColXpr as inputs.
// Following the website http://eigen.tuxfamily.org/dox/TopicFunctionTakingEigenTypes.html
// this is realised with a little hack, passing a constant x and casting the constness away.
// The function still has to be tested for performance.
template <typename Derived>
void SparseSolveLxtemplated (Model &model, Math::MatrixNd &L, Eigen::MatrixBase<Derived> const &x_aux) {
  Eigen::MatrixBase<Derived> &x = const_cast <Eigen::MatrixBase<Derived> & > (x_aux);

  for (unsigned int i = 1; i <= model.qdot_size; i++) {
        unsigned int j = model.lambda_q[i];
        while (j != 0) {
            x[i - 1] = x[i - 1] - L(i - 1,j - 1) * x[j - 1];
            j = model.lambda_q[j];
        }
        x[i - 1] = x[i - 1] / L(i - 1,i - 1);
    }
}

// The following function is a replacement for the function SparseSolveLTx in src/rbdl_mathutils.cc . It can take Vectors and ColXpr as inputs for x.
// Following the website http://eigen.tuxfamily.org/dox/TopicFunctionTakingEigenTypes.html
// this is realised with a little hack, passing a constant x and casting the constness away.
// The function still has to be tested for performance.
template <typename Derived>
void SparseSolveLTxtemplated (Model &model, Math::MatrixNd &L, Eigen::MatrixBase<Derived> const &x_aux) {
  Eigen::MatrixBase<Derived> &x = const_cast <Eigen::MatrixBase<Derived> & > (x_aux);

  for (int i = model.qdot_size; i > 0; i--) {
        x[i - 1] = x[i - 1] / L(i - 1,i - 1);
        unsigned int j = model.lambda_q[i];
        while (j != 0) {
            x[j - 1] = x[j - 1] - L(i - 1,j - 1) * x[i - 1];
            j = model.lambda_q[j];
        }
    }
}




/*This functin shall compute the Forward Dynamics by solving M qddot = tau - N(q,qdot)*/
void ForwardDynamicsCholesky (
                  Model &model,
                  const VectorNd &q,
                  const VectorNd &qdot,
                  const VectorNd &tau,
                  VectorNd &qddot,
                  std::vector<SpatialVector>* f_ext
                  ){

  VectorNd zero_vector (VectorNd::Zero (tau.size()));

  //Here we geht the zero_vector by calling Inverse dynamics with qddot = 0
  InverseDynamics(model,q,qdot,zero_vector,qddot,f_ext);
  // cout << "qddot: " << qddot << endl;
  // cout << "zero_vector: " << zero_vector << endl;
  // cout << "tau: " << tau << endl;

  MatrixNd M ( MatrixNd::Zero(qddot.size(),qddot.size()));
  CompositeRigidBodyAlgorithm (model, q, M);  //bool update kinematics?? Works without it
  //  cout << "Mass Matrix: " << M << endl;

  SparseFactorizeLTL(model,M);

  qddot = tau - qddot;
  SparseSolveLTxtemplated(model,M,qddot);
  SparseSolveLxtemplated(model,M,qddot);

};



void ad_ForwardDynamicsCholesky (
                 Model& model,
                 ADModel& ad_model,
                 const VectorNd& q,
                 const MatrixNd& q_dirs,
                 const VectorNd& qdot,
                 const MatrixNd& qdot_dirs,
                 const VectorNd& tau,
                 const MatrixNd& tau_dirs,
                 VectorNd& qddot,
                 MatrixNd& ad_qddot,
                 std::vector<SpatialVector>* f_ext) {

  unsigned int ndirs = q_dirs.cols();
  ad_model.resize_directions(ndirs);


  VectorNd zero_vector (VectorNd::Zero (q.size()));
  MatrixNd zero_vector_dirs (MatrixNd::Zero(q.size(),ndirs));

  /*Here we get the coriolis term by calling Inverse dynamics with accelerations = 0.
    The result is stored in qddot. */

  ad_InverseDynamics(model,
             ad_model,
             q,
             q_dirs,
             qdot,
             qdot_dirs,
             zero_vector,
             zero_vector_dirs,
             qddot,
             ad_qddot,
             f_ext);

  MatrixNd M (MatrixNd::Zero(qddot.size(),qddot.size()));
  std::vector<MatrixNd> ad_M (q_dirs.cols(),MatrixNd::Zero(model.q_size,model.q_size));

  ad_CompositeRigidBodyAlgorithm(model,
                 ad_model,
                 q,
                 q_dirs,
                 M,
                 ad_M);

  //    CompositeRigidBodyAlgorithm (model, q, M);  //bool update kinematics?? Works without it

  //  cout << "Mass Matrix: " << M << endl;

  SparseFactorizeLTL(model,M);

  qddot = tau - qddot;
  SparseSolveLTx(model,M,qddot);
  SparseSolveLx(model,M,qddot);


  // Now there is a templated version of SparseSolveLTx and SparseSolveLx and we don't need the local copy anymore. This still has to be tested for performance of the local copy against the templated version.
  //VectorNd local_column_copy;

  for (unsigned int j = 0; j < ndirs; j++)
    {
      ad_qddot.col(j)=tau_dirs.col(j) - ad_qddot.col(j)-ad_M[j]*qddot;
      //local_column_copy=ad_qddot.col(j);
      SparseSolveLTxtemplated(model,M,ad_qddot.col(j));
      SparseSolveLxtemplated(model,M,ad_qddot.col(j));
      //ad_qddot.col(j)=local_column_copy;
    }


};


