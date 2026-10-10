#ifndef MTQFSW_ADCS_EKF_RMM_H
#define MTQFSW_ADCS_EKF_RMM_H
#include "adcs/config.h"

struct RmmEkf {
   double x[6];
   double P[6][6];
   double Phi[6][6];
   double Q[6][6];
   double Qmin[6][6];
   double R[3][3];
   double alpha;
   long initialized;
};

void RmmEkfInit(struct RmmEkf *ekf, const Ekf_config_t *ekf_config);

int RmmEkfStep(
   struct RmmEkf *ekf,
   const double magneticCommand[3],
   const double inertia[3][3],
   const double magneticField[3],
   const double angularRateMeasurement[3],
   double dt
);

#endif