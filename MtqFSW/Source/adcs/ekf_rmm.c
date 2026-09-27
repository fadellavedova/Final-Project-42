#include "adcs/ekf_rmm.h"
#include "abMATH.h"

#include <math.h>
#include <string.h>

#define EKF_STATE_SIZE 6
#define EKF_MEASUREMENT_SIZE 3

void RmmEkfInit(
   struct RmmEkf *ekf,
   const double initialState[6],
   const double initialCovariance[6][6],
   const double processNoise[6][6],
   const double minimumProcessNoise[6][6],
   const double measurementNoise[3][3],
   double alpha)
{
   memcpy(ekf->x,initialState,sizeof(ekf->x));
   memcpy(ekf->P,initialCovariance,sizeof(ekf->P));
   memcpy(ekf->Q,processNoise,sizeof(ekf->Q));
   memcpy(ekf->Qmin,minimumProcessNoise,sizeof(ekf->Qmin));
   memcpy(ekf->R,measurementNoise,sizeof(ekf->R));
   ekf->alpha = alpha;
   ekf->initialized = 1;
}

int RmmEkfStep(
   struct RmmEkf *ekf,
   const double magneticCommand[3],
   const double inertia[3][3],
   const double magneticField[3],
   const double angularRateMeasurement[3],
   double dt
)
{
   const double velocityGain = 10.0;
   const double epsilon = 0.001;
   double inertiaInverse[3][3];
   double angularMomentum[3];
   double commandTorque[3];
   double residualTorque[3];
   double gyroscopicTorque[3];
   double netTorque[3];
   double angularAcceleration[3];
   double angularRateSkew[3][3];
   double angularMomentumSkew[3][3];
   double magneticFieldSkew[3][3];
   double temporary3a[3][3];
   double temporary3b[3][3];
   double linearizedDynamics[6][6] = {{0.0}};
   double scaledDynamics[6][6];
   double gamma[6][6];
   double predictedState[6];
   double predictedCovariance[6][6];
   double temporary6a[6][6];
   double temporary6b[6][6];
   double innovationCovariance[3][3];
   double innovationCovarianceInverse[3][3];
   double kalmanGain[6][3];
   double innovation[3];
   double correctedInnovation[3];
   double stateCorrection[6] = {0.0};
   double covarianceFactor[6][6] = {{0.0}};
   double correctedCovariance[6][6];
   double fieldNormSquared = 0.0;
   long row;
   long column;
   long inner;

   if (!ekf || !ekf->initialized || dt <= 0.0) return 0;
   if (!InvertMatrix3(inertia,inertiaInverse)) return 0;

   for (row = 0; row < 3; row++) {
      fieldNormSquared += magneticField[row]*magneticField[row];
   }
   if (fieldNormSquared < 1.0E-30) return 0;

   Matrix3Vector(inertia,ekf->x,angularMomentum);
   CrossProduct(magneticCommand,magneticField,commandTorque);
   CrossProduct(&ekf->x[3],magneticField,residualTorque);
   CrossProduct(ekf->x,angularMomentum,gyroscopicTorque);
   for (row = 0; row < 3; row++) {
      netTorque[row] = commandTorque[row] + residualTorque[row]
         - gyroscopicTorque[row];
   }
   Matrix3Vector(inertiaInverse,netTorque,angularAcceleration);

   Skew(ekf->x,angularRateSkew);
   Skew(angularMomentum,angularMomentumSkew);
   Skew(magneticField,magneticFieldSkew);
   Matrix3Multiply(angularRateSkew,inertia,temporary3a);
   Matrix3Multiply(magneticFieldSkew,inertia,temporary3b);
   Matrix3Multiply(temporary3b,magneticFieldSkew,temporary3a);
   for (row = 0; row < 3; row++) {
      for (column = 0; column < 3; column++) {
         double rateTerm = 0.0;
         double magneticTerm = temporary3a[row][column]
            *epsilon*velocityGain/fieldNormSquared;
         for (inner = 0; inner < 3; inner++) {
            rateTerm += angularRateSkew[row][inner]*inertia[inner][column];
         }
         temporary3b[row][column] = angularMomentumSkew[row][column]
            - rateTerm - magneticTerm;
      }
   }
   Matrix3Multiply(inertiaInverse,temporary3b,temporary3a);
   Matrix3Multiply(inertiaInverse,magneticFieldSkew,temporary3b);
   for (row = 0; row < 3; row++) {
      for (column = 0; column < 3; column++) {
         linearizedDynamics[row][column] = temporary3a[row][column];
         linearizedDynamics[row][column + 3] = -temporary3b[row][column];
      }
   }

   for (row = 0; row < 6; row++) {
      predictedState[row] = ekf->x[row];
      if (row < 3) predictedState[row] += angularAcceleration[row]*dt;
      for (column = 0; column < 6; column++) {
         scaledDynamics[row][column] = linearizedDynamics[row][column]*dt;
      }
   }
   MatrixExponential6(scaledDynamics,ekf->Phi);
   for (row = 0; row < 6; row++) {
      for (column = 0; column < 6; column++) {
         gamma[row][column] = dt*ekf->Phi[row][column];
      }
   }

   Matrix6Multiply(ekf->Phi,ekf->P,temporary6a);
   Matrix6MultiplyTranspose(temporary6a,ekf->Phi,predictedCovariance);
   Matrix6Multiply(gamma,ekf->Q,temporary6a);
   Matrix6MultiplyTranspose(temporary6a,gamma,temporary6b);
   for (row = 0; row < 6; row++) {
      for (column = 0; column < 6; column++) {
         predictedCovariance[row][column] += temporary6b[row][column];
      }
   }

   for (row = 0; row < 3; row++) {
      for (column = 0; column < 3; column++) {
         innovationCovariance[row][column] =
            predictedCovariance[row][column] + ekf->R[row][column];
      }
   }
   if (!InvertMatrix3(innovationCovariance,innovationCovarianceInverse)) return 0;

   for (row = 0; row < 6; row++) {
      for (column = 0; column < 3; column++) {
         kalmanGain[row][column] = 0.0;
         for (inner = 0; inner < 3; inner++) {
            kalmanGain[row][column] += predictedCovariance[row][inner]
               *innovationCovarianceInverse[inner][column];
         }
      }
   }

   for (row = 0; row < 3; row++) {
      innovation[row] = angularRateMeasurement[row] - predictedState[row];
   }
   for (row = 0; row < 6; row++) {
      for (column = 0; column < 3; column++) {
         stateCorrection[row] += kalmanGain[row][column]*innovation[column];
      }
      ekf->x[row] = predictedState[row] + stateCorrection[row];
   }

   for (row = 0; row < 6; row++) {
      covarianceFactor[row][row] = 1.0;
      for (column = 0; column < 3; column++) {
         covarianceFactor[row][column] -= kalmanGain[row][column];
      }
   }
   Matrix6Multiply(covarianceFactor,predictedCovariance,correctedCovariance);
   memcpy(ekf->P,correctedCovariance,sizeof(ekf->P));

   for (row = 0; row < 3; row++) {
      correctedInnovation[row] = angularRateMeasurement[row] - ekf->x[row];
   }
   for (row = 0; row < 6; row++) {
      stateCorrection[row] = 0.0;
      for (column = 0; column < 3; column++) {
         stateCorrection[row] += kalmanGain[row][column]
            *correctedInnovation[column];
      }
   }
   for (row = 0; row < 6; row++) {
      for (column = 0; column < 6; column++) {
         ekf->Q[row][column] = (1.0 - ekf->alpha)*ekf->Q[row][column]
            + ekf->alpha*(ekf->Qmin[row][column]
            + stateCorrection[row]*stateCorrection[column]);
      }
   }

   return 1;
}