#include "42.h"
#include <math.h>
#include "abMATH.h"

// Constants
const double pi=3.14159265e0;
const double DEGTORAD = pi/180.;
const double RADTODEG = 180./pi;
const double RADSECTORPM = 30./pi;
const double RPMTORADSEC = pi/30.;
const double EARTHGRAVCO = 3.986004418e14; // in METER^3/SECOND^2 (es mu=GM)
const double RE = 6.378135e06;             // in METERS
const double RE2 = RE*RE;
const double EARTHOMEGA = 7.2921151467e-5; // in RADIANS/SECOND
const double EARTH_J2 =  1.082616e-3;
const double EARTH_J3 = -2.538810e-6;
const double EARTH_J4 = -1.655970e-6;

double earthradius = RE;


  /* ********************************************************************
   UBA
/* ********************************************************************

/* ********************************************************************
   ADCS based on TAM, CSS, Gyros, Horizon Sensor and MTQ
   Detumbling: 
      TAM, Gyros and MTQ
   Sun Pointing: 
      TAM, Gyros, CSS and MTQ
   Nadir Pointing: 
      TAM, Gyros, CSS, Horizon Sensor and MTQ
  ********************************************************************/
  // double framefactor=-1;
  int adcsMagSunUBA(struct AcType *AC)
  {
  
   struct AcCfsCtrlType *C;
   C = &AC->CfsCtrl;

   int retval = 0.;
   static int secondspointed = 0, secondssun = 0, spinseconds = 0, timepostdetumbling = 0;
   double qe[4]={0.,0.,0.,1.};

   timepostdetumbling += AC->DT;

   // Matriz Gamma promedio móvil
   static double GAV[3][3] = {{0,0,0}, {0,0,0}, {0,0,0}};

   double mbvb = MAGV(AC->bvb);        // Norm of B
   static double kq = 0., kw = 10.;    // Initial Gains
   int coriolis = 1;                   // == 1 => coriolis compensation
   int iwearth = 1;                    // >= 0 => M3, axis of rotation, x=0 y=1 z=2
   static double earthvel = 0.;        // != 0 => M3,                   earthvel * wo
   double kw0 = 20.0;                  // Derivative Gain with Eclipse
   double kw1 = 5.0;                   // Derivative Gain without Eclipse
   double kq0 = 0.075;                 // Proportional Gain with Eclipse
   double kq1 = 0.0075;                // Proportional Gain without Eclipse
   double eps = 0.001;                 // Averaging Parameter
   double thnw = 0.001;                // End of Detumbling Velocity Condition [rad/sec] 
   double nm=0.;                       // Auxiliary Variable
   double threclipse = 0.5;            // Eclipse detection threshold (between 0 and 6, 6 never detects daylight, 0 never detects eclipse)
   double J[3] = {AC->MOI[0][0], AC->MOI[1][1], AC->MOI[2][2]};   // Inertia Matrix Diagonal
   double nw = sqrt(AC->wbn[0]*AC->wbn[0]+AC->wbn[1]*AC->wbn[1]+AC->wbn[2]*AC->wbn[2]);   // Norm of angular velocity in b
   int iwheel = 2; // x=0 to z=2 uses one wheel on this axis, otherwise does not use it (must be consistent with SC_xx.txt)
   if (AC->Nwhl==0) iwheel = -2;
   else if (AC->Whl[0].Axis[0]==1) iwheel = 0;
   else if (AC->Whl[0].Axis[1]==1) iwheel = 1;
   else if (AC->Whl[0].Axis[2]==1) iwheel = 2;
   else iwheel = -1;
   if (mbvb > 0.) {
      double m[3], u[3];
      for (int i = 0; i < 3; i++) {
         /* Angular Velocity Feedback */
         u[i] = - kw * eps * J[i] * AC->wbn[i] / mbvb / mbvb;
         if (i==iwearth && kq!=0.) {
            //if (earthvel==0) spinseconds=0; 
            //else if (spinseconds<100*60) spinseconds++;
            u[i] = u[i] + kw * eps * J[i] * 0.001082 * earthvel / mbvb / mbvb;
         }
         /* Sun Pointing Feedback */
         double sumillum = 0.;
         for (int ss=0; ss<AC->Ncss; ss++) sumillum += AC->CSS[ss].Illum;
         if (sumillum>threclipse) {
            /* Sun with AC->svb */
            secondssun++;
            double svb[3] = {AC->svb[0], AC->svb[1], AC->svb[2]};
            double qc[4] = {-svb[2], 0., svb[0], svb[1]};
            
            static int secondschange = 0;
            static double signq = 1.;
            if (svb[1]<0. && secondschange>100*60) {
               signq = -1.;
               secondschange=0;
            } if (svb[1]>0. && secondschange>100*60) {
               signq = 1.;
               secondschange=0;
            } else {
               secondschange++;
            }
            qc[0] = signq*qc[0];
            qc[1] = signq*qc[1];
            qc[2] = signq*qc[2];
            
            if (nw < thnw && kq==0.) {kq = kq0; kw = kw0;}
            if (kq!=0.) {
               if (secondssun>2*100*60) {
                  kq = kq1;
                  kw = kw1;
               } else {
                  kq = kq0;
                  kw = kw0;
               }
            }
            if (nw < thnw && fabs(svb[1])>0.98) {
               secondspointed++;
            } else if (fabs(svb[1])<0.25) {
               secondspointed=0.;
            }
            u[i] = u[i] - kq * eps * eps * (1/J[i]) * qc[i] / mbvb / mbvb;
            if (i==iwheel && kq!=0.) {
               AC->Tcmd[i] = - kw * eps * J[i] * AC->wbn[i] - kq * eps * eps * (1/J[i]) * qc[i];
               if (i==iwearth) AC->Tcmd[i] += kw * eps * J[i] * 0.001082 * earthvel;
               AC->Tcmd[i] = 10.0*AC->Tcmd[i];
               u[i] = 0.; // Descargar momento de ruedas
            }
            qe[0]=qc[0];qe[1]=qc[1];qe[2]=qc[2];qe[3]=qc[3];
         } else {
            if(kq!=0.) {
               if (i==iwearth) u[i] = u[i] + kw * eps * J[i] * 0.001082 * earthvel / mbvb / mbvb; 
               u[i] = 0.5*u[i];
               if (i==iwheel) {
                  AC->Tcmd[i] = - kw * eps * J[i] * AC->wbn[i];
                  if (i==iwearth) AC->Tcmd[i] += kw * eps * J[i] * 0.001082 * earthvel;
                  AC->Tcmd[i] = 10.0*AC->Tcmd[i];
                  u[i] = 0.;
               }
            }
            secondspointed=0.;
            secondssun=0;
         }
      }
      /* m = b x u */
      m[0] = -(u[1]*AC->bvb[2]-u[2]*AC->bvb[1]);
      m[1] = -(u[2]*AC->bvb[0]-u[0]*AC->bvb[2]);
      m[2] = -(u[0]*AC->bvb[1]-u[1]*AC->bvb[0]);
      /* Direct Control Allocation */
      int dnm = 0;
      nm = sqrt(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);
      if (fabs(m[0])*(1-dnm)+nm*dnm>AC->MTB[0].Mmax) {
         m[0]= m[0]*AC->MTB[0].Mmax/(nm*dnm+(1-dnm)*fabs(m[0]));
         m[1]= m[1]*AC->MTB[0].Mmax/(nm*dnm+(1-dnm)*fabs(m[0]));
         m[2]= m[2]*AC->MTB[0].Mmax/(nm*dnm+(1-dnm)*fabs(m[0]));
      }
      nm = sqrt(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);
      if (fabs(m[1])*(1-dnm)+nm*dnm>AC->MTB[1].Mmax) {
         m[0]= m[0]*AC->MTB[1].Mmax/(nm*dnm+(1-dnm)*fabs(m[1]));
         m[1]= m[1]*AC->MTB[1].Mmax/(nm*dnm+(1-dnm)*fabs(m[1]));
         m[2]= m[2]*AC->MTB[1].Mmax/(nm*dnm+(1-dnm)*fabs(m[1]));
      }
      nm = sqrt(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);
      if (fabs(m[2])*(1-dnm)+nm*dnm>AC->MTB[2].Mmax) {
         m[0]= m[0]*AC->MTB[2].Mmax/(nm*dnm+(1-dnm)*fabs(m[2]));
         m[1]= m[1]*AC->MTB[2].Mmax/(nm*dnm+(1-dnm)*fabs(m[2]));
         m[2]= m[2]*AC->MTB[2].Mmax/(nm*dnm+(1-dnm)*fabs(m[2]));
      }
      /* .. Momentum Management */
      C->Kunl = 0.5E5;
      if (iwheel<2) C->Kunl = 2.5E5;
      double um[3]={0.,0.,0.}, mm[3]={0.,0.,0.}, tm[3]={0.,0.,0.};
      if (fabs(m[0])<AC->MTB[0].Mmax/2. && fabs(m[1])<AC->MTB[1].Mmax/2. && fabs(m[2])<AC->MTB[2].Mmax/2.) {
         if (iwheel>=0 && iwheel<3) {
            um[iwheel] -= C->Kunl*AC->Whl[0].H;
            mm[0] = -(um[1]*AC->bvb[2]-um[2]*AC->bvb[1]);
            mm[1] = -(um[2]*AC->bvb[0]-um[0]*AC->bvb[2]);
            mm[2] = -(um[0]*AC->bvb[1]-um[1]*AC->bvb[0]);
            tm[0] = (mm[1]*AC->bvb[2]-mm[2]*AC->bvb[1]);
            tm[1] = (mm[2]*AC->bvb[0]-mm[0]*AC->bvb[2]);
            tm[2] = (um[0]*AC->bvb[1]-mm[1]*AC->bvb[0]);
            AC->Tcmd[iwheel] -= tm[iwheel];
         }

         double vii[3] = {0.,0.,0.}, MM[3][3], M[3];
         static int numgamma = 0;
         numgamma ++;
         for (int ii=0;ii<3;ii++) {
            vii[ii]     = 1;
            M[0]        = - ( vii[1]*AC->bvb[2] - vii[2]*AC->bvb[1] );
            M[1]        = - ( vii[2]*AC->bvb[0] - vii[0]*AC->bvb[2] );
            M[2]        = - ( vii[0]*AC->bvb[1] - vii[1]*AC->bvb[0] );
            MM[0][ii]   =   ( M[1]*AC->bvb[2] - M[2]*AC->bvb[1] ) / mbvb / mbvb;
            MM[1][ii]   =   ( M[2]*AC->bvb[0] - M[0]*AC->bvb[2] ) / mbvb / mbvb;
            MM[2][ii]   =   ( M[0]*AC->bvb[1] - M[1]*AC->bvb[0] ) / mbvb / mbvb;
            GAV[0][ii]  = 0.0001 * MM[0][ii] + 0.9999 * GAV[0][ii];
            GAV[1][ii]  = 0.0001 * MM[1][ii] + 0.9999 * GAV[1][ii];
            GAV[2][ii]  = 0.0001 * MM[2][ii] + 0.9999 * GAV[2][ii];
            vii[ii]     = 0;
         }
         
         if (coriolis==1 && kq!=0. && iwheel>0) {
            double cor[3] = { (J[1]-J[2])*AC->wbn[1]*AC->wbn[2]-AC->wbn[1]*AC->Whl[iwheel].H, (J[2]-J[0])*AC->wbn[2]*AC->wbn[0]+AC->wbn[0]*AC->Whl[iwheel].H, (J[0]-J[1])*AC->wbn[0]*AC->wbn[1]};
            //double cor[3] = { -AC->wbn[1]*AC->Whl[iwheel].H , AC->wbn[0]*AC->Whl[iwheel].H, 0.};
            AC->Tcmd[iwheel] -= cor[iwheel];
            cor[iwheel] = 0;
            mm[0] += (cor[1]*AC->bvb[2]-cor[2]*AC->bvb[1])/0.5;
            mm[1] += (cor[2]*AC->bvb[0]-cor[0]*AC->bvb[2])/0.9;
            mm[2] += (cor[0]*AC->bvb[1]-cor[1]*AC->bvb[0])/0.5;
         }
      }
      for (int i = 0; i < 3; i++) {
         AC->Mcmd[i] = m[i] + mm[i];
      }

      /*
      if (secondssun<2500*2 && secondspointed>100*2) {
         retval = 1;
         printf("\n RETVAL=1 ss = %d, sp = %d ", secondssun, secondspointed);
      }
      */
      
      // Adaptation of drag model
      if (iwheel>=0 && kq!=0.) {
         double chat0 = AC->Whl[0].Tmax/AC->Whl[0].Hmax/AC->Whl[0].Hmax;
         static double chat = 0.00023/0.00177/0.00177;
         double pchat = 1E8;
         double Gamma[3] = {0.55,0.9,0.55}; 
         chat = chat + pchat*AC->Whl[0].H*fabs(AC->Whl[0].H)*J[iwheel]*AC->wbn[iwheel]/Gamma[iwheel];
         if (chat>1.5*chat0) chat=1.15*chat0;
         if (chat<0.5*chat0) chat=0.85*chat0;
         AC->Tcmd[iwheel] -= chat*AC->Whl[0].H*fabs(AC->Whl[0].H);
         // Saturation of wheel
         double satw = 0.1;
         if ( AC->Tcmd[iwheel] >  AC->Whl[0].Tmax*satw ) AC->Tcmd[iwheel] =  AC->Whl[0].Tmax*satw;
         if ( AC->Tcmd[iwheel] < -AC->Whl[0].Tmax*satw ) AC->Tcmd[iwheel] = -AC->Whl[0].Tmax*satw;
      }

   }

   // Mode 3 Test
   // To be sent by command, this is just an example for dawn/dusk
   if ( secondssun > 30000*2 || timepostdetumbling > 10000) {
      //earthvel = -1.0;
   }
   retval = 0.;

   static int first=1;
   FILE *FilePtr;
   if (first) {
      first=0;
      FilePtr = fopen("MTQSAT/mission.m", "w");
      if (FilePtr != NULL) {
         fprintf(FilePtr, "vm=[%f %f %f %f %f %f %f %f %f %f %f %f %f %f];\n", qe[0], qe[1], qe[2], qe[3], AC->Tcmd[iwheel], GAV[0][0], GAV[0][1], GAV[0][2], GAV[1][0], GAV[1][1], GAV[1][2], GAV[2][0], GAV[2][1], GAV[2][2]);
         fflush(FilePtr);
         fclose(FilePtr);
         printf("DEBUG: Successfully wrote to MTQSAT/mission.m\n");
      } else {
         printf("ERROR: Could not open MTQSAT/mission.m for writing (errno: %d)\n", errno);
         perror("fopen");
      }
   } else {
      FilePtr = fopen("MTQSAT/mission.m", "a");
      if (FilePtr != NULL) {
         fprintf(FilePtr, "vm=[vm;%f %f %f %f %f %f %f %f %f %f %f %f %f %f];\n", qe[0], qe[1], qe[2], qe[3], AC->Tcmd[iwheel], GAV[0][0], GAV[0][1], GAV[0][2], GAV[1][0], GAV[1][1], GAV[1][2], GAV[2][0], GAV[2][1], GAV[2][2]);
         fflush(FilePtr);
         fclose(FilePtr);
      }
   }

   return retval;
  
  }


  int adcsRwTriadTLEUBA(struct AcType *AC)
  {

   int sunwheel = 0;

   int retval = 0.;
   static int secondsmis = 0, secondsnmis = 0;

   secondsnmis++;
   if (secondsnmis>21000) {
      secondsmis++;
   }
   if (secondsmis>2500) {
      secondsmis = 0;
      secondsnmis = 0;
   }

   struct AcCfsCtrlType *C;
   double L1[3],L2[3],L3[3];
   double HxB[3];
   double AngErr;
   long i,j;

   C = &AC->CfsCtrl;

   if (C->Init) {
      C->Init = 0;
      for(i=0;i<3;i++) FindPDGains(AC->MOI[i][i],0.1,0.7,&C->Kr[i],&C->Kp[i]);
      C->Kunl = 1.0E6;
   }

/* .. Commanded Attitude */
   double qln2[4]={0.,0.,0.,1.}, qbr2[4]={0.,0.,0.,1.};
   if (AC->GPS[0].Valid) {
      CopyUnitV(AC->PosN,L3);
      VxV(AC->PosN,AC->VelN,L2);
      UNITV(L2);
      UNITV(L3);
      //for(i=0;i<3;i++) {
      //   L2[i] = -L2[i];
      //   L3[i] = -L3[i];
      //}
      VxV(L2,L3,L1);
      UNITV(L1);
      double phi=0.*3.1416/180.;
      if (secondsmis==0) phi=0.*3.1416/180.;
      for(i=0;i<3;i++) {
         AC->CLN[0][i] = L1[i]*cos(phi)-L3[i]*sin(phi);
         AC->CLN[1][i] = L2[i];
         AC->CLN[2][i] = L3[i]*cos(phi)+L1[i]*sin(phi);
      }
      C2Q(AC->CLN,AC->qln);
      double phi2=15.*3.1416/180.;
      if (secondsmis==0) phi2=0.*3.1416/180.;
      for(i=0;i<3;i++) {
         AC->CLN[0][i] = L1[i]*cos(phi2)-L3[i]*sin(phi2);
         AC->CLN[1][i] = L2[i];
         AC->CLN[2][i] = L3[i]*cos(phi2)+L1[i]*sin(phi2);
      }
      C2Q(AC->CLN,qln2);
      AC->wln[1] = MAGV(AC->VelN)/MAGV(AC->PosN);
   } else {
      for(i=0;i<3;i++) {
         for(j=0;j<3;j++) {
            AC->CLN[i][j] = 0.0;
         }
         AC->CLN[i][i] = 1.0;
         AC->qln[i] = 0.0;
         AC->wln[i] = 0.0;
      }
      AC->qln[3] = 1.0;
   }
            
/* .. Attitude Control Reference */
   QxQT(AC->qbn,AC->qln,AC->qbr);
   RECTIFYQ(AC->qbr);

   QxQT(AC->qbn,qln2,qbr2);
   RECTIFYQ(qbr2);

   static double sqc[4]={0.,0.,0.,0.};

   double mbvb = MAGV(AC->bvb);        // Norm of B
   static double kq = 0., kw = 10.;    // Initial Gains
   int coriolis = 0;                   // == 1 => coriolis compensation
   int iwearth = 1;                    // >= 0 => axis of rotation, x=0 y=1 z=2
   static double earthvel = 1.0;       // != 0 => earthvel * wo
   static double earthang = 0.0;       // != 0 => shift
   double kw0 = 20.0;                  // Derivative Gain with Eclipse
   double kw1 = 5.0;                   // Derivative Gain without Eclipse
   double kq0 = 0.075;                 // Proportional Gain with Eclipse
   double kq1 = 0.0075;                // Proportional Gain without Eclipse
   double eps = 0.001;                 // Averaging Parameter
   if (secondsmis>0) eps = 0.002;
   double thnw = 0.001;                // End of Detumbling Velocity Condition [rad/sec] 
   double nm=0.;                       // Auxiliary Variable
   double threclipse = 0.5;            // Eclipse detection threshold (between 0 and 6, 6 never detects daylight, 0 never detects eclipse)
   double J[3] = {AC->MOI[0][0], AC->MOI[1][1], AC->MOI[2][2]};   // Inertia Matrix Diagonal
   double nw = sqrt(AC->wbn[0]*AC->wbn[0]+AC->wbn[1]*AC->wbn[1]+AC->wbn[2]*AC->wbn[2]);   // Norm of angular velocity in b
   int iwheel = 2; // x=0 to z=2 uses one wheel on this axis, otherwise does not use it (must be consistent with SC_xx.txt)
   if (AC->Nwhl==0) iwheel = -2;
   else if (AC->Whl[0].Axis[0]==1) iwheel = 0;
   else if (AC->Whl[0].Axis[1]==1) iwheel = 1;
   else if (AC->Whl[0].Axis[2]==1) iwheel = 2;
   else iwheel = -1;
   if (mbvb > 0.) {
      double m[3], u[3];
      for (int i = 0; i < 3; i++) {
         /* Angular Velocity Feedback */
         u[i] = - kw * eps * J[i] * AC->wbn[i] / mbvb / mbvb;
         if (i==iwearth && kq!=0.) u[i] = u[i] + kw * eps * J[i] * 0.001082 * earthvel / mbvb / mbvb;
         if (i==iwearth && kq!=0.) u[i] = u[i] + kq * eps * eps / J[i] * earthang / mbvb / mbvb;
         double svb[3] = {AC->svb[0], AC->svb[1], AC->svb[2]};
         //double qc[4] = {AC->qbr[0], AC->qbr[1], AC->qbr[2], AC->qbr[3]};
         double qc[4] = {qbr2[0], qbr2[1], qbr2[2], qbr2[3]};
         double sqc3=1.;if (qc[3]<0.) sqc3=-1.;
         qc[0]=qc[0]*sqc3;
         qc[1]=qc[1]*sqc3;
         qc[2]=qc[2]*sqc3;
         qc[3]=qc[3]*sqc3;
         kq = kq0;
         kw = kw0;
         double ki = kq*0.;
         sqc[i] = 0.001*qc[i] + 0.999*sqc[i];
         u[i] = u[i] - kq * eps * eps * (1/J[i]) * qc[i] / mbvb / mbvb;
         u[i] = u[i] - ki * eps * eps * (1/J[i]) * sqc[i] / mbvb / mbvb;
         if (i==iwheel && kq!=0.) {
            AC->Tcmd[i] = - kw * eps * J[i] * AC->wbn[i] - kq * eps * eps * (1/J[i]) * (qc[i]+sqc[i]*ki/kq);
            if (i==iwearth) AC->Tcmd[i] += kw * eps * J[i] * AC->wln[1] * earthvel;
            AC->Tcmd[i] = 10.0*AC->Tcmd[i];
            //u[i] = 0.; // Descargar momento de ruedas
         }
      }
      /* m = b x u */
      m[0] = -(u[1]*AC->bvb[2]-u[2]*AC->bvb[1]);
      m[1] = -(u[2]*AC->bvb[0]-u[0]*AC->bvb[2]);
      m[2] = -(u[0]*AC->bvb[1]-u[1]*AC->bvb[0]);
      /* Direct Control Allocation */
      if (fabs(m[0])>AC->MTB[0].Mmax) {
         m[0]= m[0]*AC->MTB[0].Mmax/fabs(m[0]);
         m[1]= m[1]*AC->MTB[0].Mmax/fabs(m[0]);
         m[2]= m[2]*AC->MTB[0].Mmax/fabs(m[0]);
      }
      if (fabs(m[1])>AC->MTB[1].Mmax) {
         m[0]= m[0]*AC->MTB[1].Mmax/fabs(m[1]);
         m[1]= m[1]*AC->MTB[1].Mmax/fabs(m[1]);
         m[2]= m[2]*AC->MTB[1].Mmax/fabs(m[1]);
      }
      if (fabs(m[2])>AC->MTB[2].Mmax) {
         m[0]= m[0]*AC->MTB[2].Mmax/fabs(m[2]);
         m[1]= m[1]*AC->MTB[2].Mmax/fabs(m[2]);
         m[2]= m[2]*AC->MTB[2].Mmax/fabs(m[2]);
      }
      /* .. Momentum Management */
      C->Kunl = 0.5E5;
      if (iwheel<2) C->Kunl = 2.5E5;
      double um[3]={0.,0.,0.}, mm[3]={0.,0.,0.}, tm[3]={0.,0.,0.};
      if (fabs(m[0])<AC->MTB[0].Mmax/2. && fabs(m[1])<AC->MTB[1].Mmax/2. && fabs(m[2])<AC->MTB[2].Mmax/2.) {
         if (iwheel>=0 && iwheel<3) {
            um[iwheel] -= C->Kunl*AC->Whl[0].H;
            mm[0] = -(um[1]*AC->bvb[2]-um[2]*AC->bvb[1]);
            mm[1] = -(um[2]*AC->bvb[0]-um[0]*AC->bvb[2]);
            mm[2] = -(um[0]*AC->bvb[1]-um[1]*AC->bvb[0]);
            tm[0] = (mm[1]*AC->bvb[2]-mm[2]*AC->bvb[1]);
            tm[1] = (mm[2]*AC->bvb[0]-mm[0]*AC->bvb[2]);
            tm[2] = (um[0]*AC->bvb[1]-mm[1]*AC->bvb[0]);
            AC->Tcmd[iwheel] -= tm[iwheel];
         }
         if (coriolis==1 && kq!=0.) {
            double cor[3] = { (J[1]-J[2])*AC->wbn[1]*AC->wbn[2], (J[2]-J[0])*AC->wbn[2]*AC->wbn[0], (J[0]-J[1])*AC->wbn[0]*AC->wbn[1]};
            AC->Tcmd[iwheel] -= cor[iwheel];
            cor[iwheel]=0;
            mm[0] += (cor[1]*AC->bvb[2]-cor[2]*AC->bvb[1]);
            mm[1] += (cor[2]*AC->bvb[0]-cor[0]*AC->bvb[2]);
            mm[2] += (cor[0]*AC->bvb[1]-cor[1]*AC->bvb[0]);
         }
      }
      for (int i = 0; i < 3; i++) {
         AC->Mcmd[i] = m[i] + mm[i];
      }
      
      // Adaptation of drag model
      if (iwheel>=0 && kq!=0.) {
         double chat0 = AC->Whl[0].Tmax/AC->Whl[0].Hmax/AC->Whl[0].Hmax;
         static double chat = 0.00023/0.00177/0.00177;
         double pchat = 1E8;
         chat = chat + pchat*AC->Whl[0].H*fabs(AC->Whl[0].H)*J[iwheel]*AC->wbn[iwheel];
         if (chat>1.15*chat0) chat=1.15*chat0;
         if (chat<0.85*chat0) chat=0.85*chat0;
         AC->Tcmd[iwheel] -= chat*AC->Whl[0].H*fabs(AC->Whl[0].H);
         // Saturation of wheel
         double satw = 0.1;
         if ( AC->Tcmd[iwheel] >  AC->Whl[0].Tmax*satw ) AC->Tcmd[iwheel] =  AC->Whl[0].Tmax*satw;
         if ( AC->Tcmd[iwheel] < -AC->Whl[0].Tmax*satw ) AC->Tcmd[iwheel] = -AC->Whl[0].Tmax*satw;
      }

   }

   static int first=1;
   FILE *FilePtr;
   if (first) {
            first=0;
            FilePtr = fopen("TSAT/mission.m", "w");
            fprintf(FilePtr, "vm=[%f %f %f %f %f];\n", AC->qbr[0], AC->qbr[1], AC->qbr[2], AC->qbr[3], AC->Tcmd[iwheel]);
            fclose(FilePtr);
   } else {
            FilePtr = fopen("TSAT/mission.m", "a");
            fprintf(FilePtr, "vm=[vm;%f %f %f %f %f];\n", AC->qbr[0], AC->qbr[1], AC->qbr[2], AC->qbr[3], AC->Tcmd[iwheel]);
            fclose(FilePtr);
   }

   return retval;
 
  }