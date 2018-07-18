#ifndef SOLVER_HPP_INCLUDED
#define SOLVER_HPP_INCLUDED

void SOLVER(int Gi, int Gj, int Gk,
            gsl_matrix * JACOBIAN, gsl_matrix * C_M_WGO, gsl_matrix * QR_N,
            gsl_vector * FN_V, gsl_vector * C_N_WGO, gsl_vector * C_V_WGO, gsl_vector * Q_V, gsl_vector * X_V, gsl_vector * residual);

#endif // SOLVER_HPP_INCLUDED
