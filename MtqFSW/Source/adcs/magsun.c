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
const double EYE[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

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
         if (Css->Illum / Css->Scale > config->eclipse_threshold)
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

void magProcessing(struct AcType *AC) {

   if (AC->Nmag != 3) {
         printf("Wrong processing function, Nmag != 3.\n");
         return;
   }

   struct AcMagnetometerType *M;
   AcConfig_t *config = GetAcConfig();
   double bvs[3];
   double bvb[3];

   for(int i = 0; i < AC->Nmag; i++) {
      M = &AC->MAG[i];
      bvs[i] = M->Field - config->mag_bias[i];
   }
   Matrix3Vector(config->mag2b_matrix, bvs, bvb);
   for(int i = 0; i < 3; i++) {
      AC->bvb[i] = bvb[i];
   }
   return;
}

void gyroProcessing(struct AcType *AC) {

   if (AC->Ngyro != 3) {
         printf("Wrong processing function, Ngyro != 3.\n");
         return;
   }

   struct AcGyroType *G;
   AcConfig_t *config = GetAcConfig();
   double wsn[3];
   double wbn[3];

   for(int i = 0; i < AC->Ngyro; i++) {
      G = &AC->Gyro[i];
      wsn[i] = G->Rate - config->gyro_bias[i];
   }
   Matrix3Vector(config->gyro2b_matrix, wsn, wbn);
   for(int i = 0; i < 3; i++) {
      AC->wbn[i] = wbn[i];
   }
   return;
}


void esProcessing(struct AcType *AC) // FUNCION: tiene que calcular el vecornadir a partir de los datos Roll y Pitch YA CALCULADOS en el modelo de sensor
{
   AcConfig_t *config = GetAcConfig();
   struct AcEarthSensorType *ES;
   double Roll, Pitch;
   int Nvalid = 0;
   double invalidNV[3] = {0.0, 0.0, 0.0};
   double nvs[3] = {0.0, 0.0, 0.0};
   double nvb[3] = {0.0, 0.0, 0.0};

   for(long i = 0; i < AC->Nehs; i++){
      ES = &AC->EHS[i];
      if(ES->Valid) {
         Nvalid++;
         nvs[0] += sin(ES->Pitch) * cos(ES->Roll);
         nvs[1] += sin(ES->Pitch) * sin(ES->Roll);
         nvs[2] += cos(ES->Pitch);         
      } 
   }

   if (AC->Nehs == 0 || Nvalid == 0)
   {
      for(int i = 0; i < 3; i++){
         AC->nvb[i] = invalidNV[i];
      }
      AC->EarthValid = FALSE;
      return;
   }
   if (Nvalid > 0) {
      AC->EarthValid = TRUE;
      Matrix3Vector(config->ehs2b_matrix, nvs, nvb);
      for(int i = 0; i < 3; i++) {
         AC->nvb[i] = nvb[i]/Nvalid;
      }
   }
}

void mtbProcessing(struct AcType *AC)
{
      if(AC->Nmtb != 3) {
         printf("Wrong processing function, Nmtb != 3.\n");
         return;
      }

      AcConfig_t *config = GetAcConfig();
      struct AcMtbType *M;
      long Im;
      double m[3];
      double max_m = 0.0;

      Matrix3Vector(config->b2mtq_matrix, AC->Mcmd, m); //TF a los ejes de los mtq

      for (size_t i = 0; i < 3; i++) {
         if (fabs(m[i]) > max_m) {
            max_m = fabs(m[i]);
         }
      }
      if (max_m > config->max_dipole) {
         double scale = config->max_dipole / max_m;
         for (size_t i = 0; i < 3; i++) {
            m[i] *= scale;
         }
      }
      for(Im=0;Im<AC->Nmtb;Im++) {
         M = &AC->MTB[Im];
         M->Mcmd = m[Im];
      }

      // Renovamos el Mcmd del AC
      Matrix3Vector(config->mtq2b_matrix, m, AC->Mcmd);

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

int adcsDetumbling(struct AcType *AC)
{
   // Detumbling Control Algorithm
   AcConfig_t *config = GetAcConfig();

   double b_norm = Norm3Vector(AC->bvb);
   double (*J)[3] = config->J;
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
      Matrix3Vector_trunc(J, (AC->ekf)->x, aux, 0, 3);
      for (size_t i = 0; i < 3; i++)
      {
         m[i] = (m[i] / (b_norm * b_norm)); // Normalize the magnetic moment command
      }

      for (size_t i = 0; i < 3; i++) {
         AC->Mcmd[i] = m[i];
      }
   }
   else
   {

      return 1;
   }

   return 0; // Return 0 to indicate successful execution
}

int adcsNadirPointing(struct AcType *AC, int *saw_horizon)
{
   // Nadir Pointing Control Algorithm
   AcConfig_t *config = GetAcConfig();

   double b_norm = Norm3Vector(AC->bvb);
   if (b_norm < 1e-9) {
      return 1;
   }

   double (*J)[3] = config->J;
   double (*inverse_J)[3] = config->inv_J;
   double d_component[3];
   double p_component[3];
   double u[3];
   double m[3];
   double qs[3];
   double qs0;
   double angle;
   double qs_norm;
   static double nvb[3];
   static int init_sp_loop = 1; // Sign of the scalar part of the quaternion representing the sun pointing error
   static int signqs0 = 1;
   double kw = config->kw_nadirpointing_nohorizon;
   double kp = config->kp_nadirpointing_nohorizon;
   double eps = config->eps_nadirpointing_nohorizon;
   double aux[3];

   // Espera a ver por primera vez el horizonte (capaz esto se puede hacer en sun pointing tmb(?))
   if (AC->EarthValid) {
      *saw_horizon = 1;
      for(int i = 0; i < 3; i++) nvb[i] = AC->nvb[i];
      kw = config->kw_nadirpointing;
      kp = config->kp_nadirpointing;
      eps = config->eps_nadirpointing;
   }

   if (*saw_horizon)
   {  
      CrossProduct(nvb, config->nadir_pointing_vector, qs); // Compute the nadir vector error
      qs_norm = Norm3Vector(qs);
      double dot = DotProduct(nvb, config->nadir_pointing_vector);
      angle = atan2(qs_norm, dot);
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
      Matrix3Vector(inverse_J, qs, p_component); // Compute the proportional component
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
      qs0 = 1.0;
   }

   Matrix3Vector_trunc(J, (AC->ekf)->x, d_component, 0, 3); // Compute the derivative component
   for (int i = 0; i < 3; i++)
   {
      u[i] = -kw * eps * d_component[i] - kp * eps * eps * signqs0 * p_component[i]; // Control law
   }
   CrossProduct(AC->bvb, u, m); // Compute the magnetic moment command for nadir pointing
   Matrix3Vector_trunc(EYE, (AC->ekf)->x, aux, 3, 6);
   for (size_t i = 0; i < 3; i++)
   {
      m[i] = (m[i] / (b_norm * b_norm)); // Normalize the magnetic moment command
      if(config->rmm_compensate) {
         m[i] -= aux[i];
      }
   }

   for (size_t i = 0; i < 3; i++) {
      AC->Mcmd[i] = m[i];
   }
   return 0;
}

int adcsSunPointing(struct AcType *AC)
{
   // Sun Pointing Control Algorithm
   AcConfig_t *config = GetAcConfig();

   double b_norm = Norm3Vector(AC->bvb);
   if (b_norm < 1e-9) {
      return 1;
   }

   double (*J)[3] = config->J;
   double (*inverse_J)[3] = config->inv_J;
   double d_component[3];
   double p_component[3];
   double u[3];
   double m[3];
   double qs[3];
   double qs0;
   double angle;
   double qs_norm;
   static int init_sp_loop = 1; // Sign of the scalar part of the quaternion representing the sun pointing error
   static int signqs0 = 1;
   double kw = config->kw_sunpointing;
   double kp = config->kp_sunpointing;
   double eps = config->eps_sunpointing;
   double aux[3];

   if (AC->SunValid)
   {
      CrossProduct(AC->svb, config->sun_pointing_vector, qs); // Compute the sun vector error
      qs_norm = Norm3Vector(qs);
      double dot = DotProduct(AC->svb, config->sun_pointing_vector);
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
      Matrix3Vector(inverse_J, qs, p_component);

      if (init_sp_loop)
      {
         signqs0 = (qs0 >= 0) ? 1.0 : -1.0; // Determine the sign of the scalar part of the quaternion
         init_sp_loop = 0;
      }
   }
   else
   {  
      kw = config->kw_sunpointing_eclipse;
      eps = config->eps_sunpointing_eclipse;

      for (int i = 0; i < 3; i++)
      {
         p_component[i] = 0.0; 
      }
      qs0 = 1.0; 
   }

   Matrix3Vector_trunc(J, (AC->ekf)->x, d_component, 0, 3);// Compute the control moment based on angular velocity and inertia
   for (int i = 0; i < 3; i++)
   {
      u[i] = -kw * eps * d_component[i] - kp * eps * eps * signqs0 * p_component[i]; // Control law  u = u - obj.eps*obj.eps*obj.k_p*inv(obj.iner)*dq*signq4;
   }
   Matrix3Vector_trunc(EYE, (AC->ekf)->x, aux, 3, 6);
   CrossProduct(AC->bvb, u, m); 
   
   for (size_t i = 0; i < 3; i++)
   {
      m[i] = (m[i] / (b_norm * b_norm));
      
      if(config->rmm_compensate) {
         m[i] -= aux[i];
      }
   }

   for (size_t i = 0; i < 3; i++) {
      AC->Mcmd[i] = m[i];
   }

   return 0;
}

int adcsUBA(struct AcType *AC)
{
   /* Get EKF configuration */

   if ((AC->ekf)->initialized != 1)
   {
      RmmEkfInit(AC->ekf);
   }

   int retval = 0.;
   static int saw_horizon = 0; // Para guardar primer avistamiento el horizonte en nad pointing

   /* Get ADCS configuration */
   AcConfig_t *config_ac = GetAcConfig();

   UpdateMode(AC, config_ac);
   enum AcMode currentMode = GetCurrentMode(AC);
   RmmEkfStep(AC, config_ac);
   switch (currentMode)
   {
   case AC_MODE_DETUMBLE:
      saw_horizon = 0;
      retval = adcsDetumbling(AC);

      break;

   case AC_MODE_SUN_POINTING:
      saw_horizon = 0;
      retval = adcsSunPointing(AC); // Hacer algo despues con el retval

      break;
   case AC_MODE_NADIR_POINTING:

      retval = adcsNadirPointing(AC, &saw_horizon);

      break;
   }
   return retval;
}
