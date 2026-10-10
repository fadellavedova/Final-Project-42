#include "adcs/magsun.h"
#include "adcs/ekf_rmm.h"
#include <math.h>
#include "abMATH.h"
#include "adcs/config.h"
#include "adcs/state_machine.h"
#include "mathkit.h"

// Constants
const double pi = 3.14159265e0;
const double DEGTORAD = pi / 180.;
const double RADTODEG = 180. / pi;
const double RADSECTORPM = 30. / pi;
const double RPMTORADSEC = pi / 30.;
const double EARTHGRAVCO = 3.986004418e14; // in METER^3/SECOND^2 (es mu=GM)
const double RE = 6.378135e06;             // in METERS
const double RE2 = RE * RE;
const double EARTHOMEGA = 7.2921151467e-5; // in RADIANS/SECOND
const double EARTH_J2 = 1.082616e-3;
const double EARTH_J3 = -2.538810e-6;
const double EARTH_J4 = -1.655970e-6;
const double  EYE[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

double earthradius = RE;

/* ********************************************************************
   UBA
/* ********************************************************************


/* ********************************************************************
   Sensor processing
   CSS, Mag, Gyros. HS
********************************************************************/

void cssProcessing(struct AcType *AC)
{
   /* Get ADCS configuration */
   AcConfig_t *config = GetAcConfig();
   struct AcCssType *Css;
   double AtA[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
   double Atb[3] = {0.0, 0.0, 0.0};
   double AtAi[3][3];
   double A[2][3], b[2];
   long Ic, i, j;
   long Nvalid = 0;
   double InvalidSVB[3] = {0.0, 0.0, 0.0}; /* Safe vector if SunValid == FALSE */

   if (AC->Ncss == 0)
   {
      /* AC->svb populated by true S->svb in 42sensors.c */
   }
   else
   {
      for (Ic = 0; Ic < AC->Ncss; Ic++)
      {
         Css = &AC->CSS[Ic];
         if (Css->Illum > config->eclipse_threshold)
         {
            Nvalid++;
            /* Normal equations, assuming Nvalid will end up > 2 */
            for (i = 0; i < 3; i++)
            {
               Atb[i] += Css->Axis[i] * Css->Illum / Css->Scale;

               for (j = 0; j < 3; j++)
               {
                  AtA[i][j] += Css->Axis[i] * Css->Axis[j];
               }
            }
            /* In case Nvalid ends up == 2 */
            for (i = 0; i < 3; i++)
            {
               A[0][i] = A[1][i];
               A[1][i] = Css->Axis[i];
            }
            b[0] = b[1];
            b[1] = Css->Illum / Css->Scale;
         }
      }
      if (Nvalid > 2)
      {
         AC->SunValid = TRUE;
         InvertMatrix3(AtA, AtAi);
         Matrix3Vector(AtAi, Atb, AC->svb);
         Normalize3Vector(AC->svb);
      }
      else if (Nvalid == 2)
      {
         AC->SunValid = TRUE;
         for (i = 0; i < 3; i++)
            AC->svb[i] = b[0] * A[0][i] + b[1] * A[1][i];
         Normalize3Vector(AC->svb);
      }
      else if (Nvalid == 1)
      {
         AC->SunValid = TRUE;
         for (i = 0; i < 3; i++)
            AC->svb[i] = Atb[i];
         Normalize3Vector(AC->svb);
      }
      else
      {
         AC->SunValid = FALSE;
         for (i = 0; i < 3; i++)
            AC->svb[i] = InvalidSVB[i];
      }
   }
}

void esProcessing(struct AcType *AC) // FUNCION: tiene que calcular el vecornadir a partir de los datos Roll y Pitch YA CALCULADOS en el modelo de sensor
{
   struct AcEarthSensorType *ES = &AC->ES;
   double Roll, Pitch;

   if (AC->Nst == 0)
   {
      /* AC->qbn populated by true S->B[0].qn in 42sensors.c */
      /* AC->CBN populated by true S->B[0].CBN in 42sensors.c */
      /* AC->CLN populated by true S->B[0].CLN in 42sensors.c */
   }
   else
   {
      /* AC->qbn populated by true S->B[0].qn in 42sensors.c */
      /* AC->CBN populated by true S->B[0].CBN in 42sensors.c */
      /* AC->CLN populated by true S->B[0].CLN in 42sensors.c */
   }

   Roll = atan2(AC->CLN[2][1], AC->CLN[2][2]);
   Pitch = -asin(AC->CLN[2][0]);

   ES->Roll = Roll;
   ES->Pitch = Pitch;
}

/* ********************************************************************
   Attitude Est and Control Algorithms
   Detumbling:
      TAM, Gyros and MTQ
   Sun Pointing:
      TAM, Gyros, CSS and MTQ
   Nadir Pointing:
      TAM, Gyros, CSS, Horizon Sensor and MTQ
  ********************************************************************/
// double framefactor=-1;
/* En algún lado tenemos que implementar el ekf_rmm.h
      en matlabtr está:
        phik = 0;
        rmm_cov_diag = [0;0;0];
        if RMM_ESTIMATE == 1
            %[x_pred, sigma, phikm1, Q] = ekf_rmm(x_pred, (mag_mom - mom_res), iner, earth_field_b, sigma, dw, Q, Q_min, alpha, R, tstep);
            [ekf_rmm, phik] = update(ekf_rmm, mag_mom - mom_res, earth_field_b, dw);
            [x_pred, sigma, Q] = get_estimates(ekf_rmm);
            rmm_cov_diag = [sigma(4,4); sigma(5,5); sigma(6,6)];
        end

   Se puede tener el objeto filtro dentro del AC y llamarlo desde ahí? Lo inicializamos en la primera iteracion con un static int que despues vaya a 0
*/

int adcsDetumbling(struct AcType *AC, struct AcConfig *config)
{
   // Detumbling Control Algorithm

   double b_norm = Norm3Vector(AC->bvb);
   double (*J)[3] = AC->MOI;
   double u[3];
   double m[3];
   double kw = config->kw_detumb;
   double eps = config->eps_detumb;
   double aux[3];

   if (b_norm > 1e-12)
   {
      Matrix3Vector_trunc(J, AC->wbn, u, 0, 3);
      for (int i = 0; i < 3; i++)
      {
         u[i] = -kw * eps * u[i]; // Control law
      }
      CrossProduct(AC->bvb, u, m); // Compute the magnetic moment command
      Matrix3Vector_trunc(J, &(AC->ekf.x), aux, 0, 3);
      for (size_t i = 0; i < 3; i++)
      {
         m[i] = (m[i] / (b_norm * b_norm)) ; // Normalize the magnetic moment command
         AC->Mcmd[i] = m[i];
      }
   }
   else
   {

      return 1;
   }

   return 0; // Return 0 to indicate successful execution
}

int adcsNadirPointing(struct AcType *AC, struct AcConfig *config)
{
   // Nadir Pointing Control Algorithm

   double b_norm = Norm3Vector(AC->bvb);
   AC -> ekf.x;

   if (b_norm > 1e-12)
   {

      double (*J)[3] = AC->MOI;
      double inverse_J[3][3];
      InvertMatrix3(AC->MOI, inverse_J);
      double (*invJ)[3] = inverse_J;
      double d_component[3];
      double p_component[3];
      double u[3];
      double m[3];
      double qs[3];
      double qs0;
      double angle;
      double qs_norm;
      static int init_sp_loop = 1; // Sign of the scalar part of the quaternion representing the sun pointing error
      int signqs0 = 1;
      double kw = config->kw_nadirpointing;
      double kp = config->kp_nadirpointing;
      double eps = config->eps_nadirpointing;
      double aux[3];

      if (Norm3Vector(AC->svb) > 0)
      {                                                            // Poner acá la condicion de nadir
         CrossProduct(AC->svb, config->nadir_pointing_vector, qs); // Compute the nadir vector error
         qs_norm = Norm3Vector(qs);
         angle = asin(qs_norm);
         qs0 = cos(angle / 2.);
         if (qs_norm > 1e-12)
         {
            for (int i = 0; i < 3; i++)
            {
               qs[i] = (qs[i] / qs_norm) * sin(angle / 2.); // Normalize the nadir vector error
            }
         }
         else
         {
            for (int i = 0; i < 3; i++)
            {
               qs[i] = 0.0;
            }
         }
         Matrix3Vector(invJ, qs, p_component); // Compute the proportional component
         if (init_sp_loop)
         {
            signqs0 = (qs0 >= 0) ? 1.0 : -1.0; // Determine the sign of the scalar part of the quaternion
            init_sp_loop = 0;
         }
      }
      else
      {
         for (int i = 0; i < 3; i++)
         {
            p_component[i] = 0.0;
         }
         qs0 = 1.0; // Set the scalar part of the quaternion to 1 (no error)
      }

      Matrix3Vector_trunc(J, &(AC->ekf.x), d_component, 0, 3); // Compute the derivative component
      for (int i = 0; i < 3; i++)
      {
         u[i] = -kw * eps * d_component[i] - kp * eps * eps * signqs0 * p_component[i]; // Control law
      }
      CrossProduct(AC->svb, u, m); // Compute the magnetic moment command for sun pointing
      Matrix3Vector_trunc(EYE, &(AC->ekf.x), aux, 3, 6);
      for (size_t i = 0; i < 3; i++)
      {
         m[i] = (m[i] / (b_norm * b_norm)) - aux[i]; // Normalize the magnetic moment command
         AC->Mcmd[i] = m[i];
      }
   }
   else
   {
      return 1;
   }

   return 0;
}

int adcsSunPointing(struct AcType *AC, struct AcConfig *config)
{
   // Sun Pointing Control Algorithm

   double b_norm = Norm3Vector(AC->bvb);
   if (b_norm > 1e-12)
   {

      double (*J)[3] = AC->MOI;
      double invJ[3][3];
      InvertMatrix3(AC->MOI, invJ);
      double d_component[3];
      double p_component[3];
      double u[3];
      double m[3];
      double qs[3];
      double qs0;
      double angle;
      double qs_norm;
      static int init_sp_loop = 1; // Sign of the scalar part of the quaternion representing the sun pointing error
      int signqs0 = 1;
      double kw = config->kw_sunpointing;
      double kp = config->kp_sunpointing;
      double eps = config->eps_sunpointing;
      double aux[3];

      if (AC->SunValid)
      {
         CrossProduct(AC->svb, config->sun_pointing_vector, qs); // Compute the sun vector error
         qs_norm = Norm3Vector(qs);
         double dot = VoV(AC -> svb, config -> sun_pointing_vector);
         angle = atan2(qs_norm, dot);
         qs0 = cos(angle / 2.);
         if (qs_norm > 1e-12)
         {
            for (int i = 0; i < 3; i++)
            {
               qs[i] = (qs[i] / qs_norm) * sin(angle / 2.); // Normalize the sun vector error
            }
         }
         else
         {
            for (int i = 0; i < 3; i++)
            {
               qs[i] = 0.0; // If the sun vector error is negligible, set it to zero
               qs0 = 1.0;
            }
         }
         Matrix3Vector(invJ, qs, p_component);

         if (init_sp_loop)
         {
            signqs0 = (qs0 >= 0) ? 1.0 : -1.0; // Determine the sign of the scalar part of the quaternion
            init_sp_loop = 0;
         }
      }
      else
      {
         for (int i = 0; i < 3; i++)
         {
            p_component[i] = 0.0; // Set the proportional component to zero
         }
         qs0 = 1.0; // Set the scalar part of the quaternion to 1 (no error)
      }

      Matrix3Vector(J, AC->wbn, d_component); // Compute the control moment based on angular velocity and inertia
      for (int i = 0; i < 3; i++)
      {
         u[i] = -kw * eps * d_component[i] - kp * eps * eps * signqs0 * p_component[i]; // Control law
      }
      printf("Sun pointing u = [%.17g, %.17g, %.17g]\n", u[0], u[1], u[2]);
      CrossProduct(AC->svb, u, m); // Compute the magnetic moment command for sun pointing
      Matrix3Vector_trunc(EYE, &(AC->ekf.x), aux, 3, 6);
      for (size_t i = 0; i < 3; i++)
      {
         m[i] = (m[i] / (b_norm * b_norm)) - aux[i];
         if (!isfinite(m[i]))
         {
            fprintf(stderr, "Sun pointing nonfinite axis=%zu b_norm=%.17g svb=%.17g u=%.17g d=%.17g p=%.17g ekf=%.17g Mcmd=%.17g\n",
                    i, b_norm, AC->svb[i], u[i], d_component[i], p_component[i], aux[i], m[i]);
         }
         AC->Mcmd[i] = m[i];
      }
   }
   else
   {
      return 1;
   }

   return 0;
}

int adcsUBA(struct AcType *AC)
{

   /* Get EKF configuration */
   Ekf_config_t *config_ekf = GetEkfConfig();
   if (AC->ekf.initialized != 1)
   {
      RmmEkfInit(&(AC->ekf), config_ekf);
   }

   int retval = 0.;

   /* Get ADCS configuration */
   AcConfig_t *config = GetAcConfig();

   UpdateMode(AC, config);
   enum AcMode currentMode = GetCurrentMode(AC);

   RmmEkfStep(&(AC -> ekf), &(AC -> Mcmd), &(AC -> MOI), &(AC -> bvb), &(AC -> wbn), AC -> DT);


   printf("%d", currentMode);
   switch (currentMode)
   {
   case AC_MODE_DETUMBLE:

      retval = adcsDetumbling(AC, config);

      break;

   case AC_MODE_SUN_POINTING:

      retval = adcsSunPointing(AC, config); // Hacer algo despues con el retval

      break;
   case AC_MODE_NADIR_POINTING:

      retval = adcsNadirPointing(AC, config);

      break;
   }
}
