#ifndef MTQFSW_ADCS_ABMATH_H
#define MTQFSW_ADCS_ABMATH_H

#ifdef __cplusplus
extern "C" {
#endif

double det4(double M[4][4]);
double det5(double M[5][5]);
double det4(double M[4][4]);
void ker45(double X[4][5],double wx[5]);
void ker56(double X[5][6],double wx[6]);
void ker67(double M[6][7],double kerM[7]);
void pinv34(double M[3][4], double piM[4][3]);
void pinv37(double M[3][7], double piM[7][3]);
void pinv45(double M[4][5], double piM[5][4]);
void pinv46(double M[4][6], double piM[6][4]);
void pinv47(double M[4][7], double piM[7][4]);
void pinv410(double M[4][10], double piM[10][4]);
void triad(double v1m[4], double v2m[4], double v1t[4], double v2t[4], double Cmt[4][4]);
void vectorialproduct(double v1[4],double v2[4],double r[4]);
void quaterniontomatrix(double q0, double q1, double q2, double q3, double C[4][4]);
void matrixtoquaternion(double C[4][4], double *pq0, double *pq1, double *pq2, double *pq3);
void rotation(double C[4][4],double vin[4],double vout[4]);
void normalize(double v[4]);
void normalizes(double v[4], int sender);
void q_mult(double q0, double q1, double q2, double q3, double q0prima, double q1prima, double q2prima, double q3prima, double *pq0prima2, double *pq1prima2, double *pq2prma2, double *pq3prima2);
void q_a_rad(double q0, double q1, double q2, double q3, double *protacion, double *pascension, double *pdeclinacion);
void rad_a_q(double *pq0, double *pq1, double *pq2, double *pq3, double rotacion, double ascension, double declinacion);
void getCd(double vpoint[4], double Cd[4][4]);
void pVec(double v1[], double v2[], double v3[]);
void CrossProduct(const double left[3], const double right[3], double result[3]);
void Skew(const double vector[3], double matrix[3][3]);
void Matrix3Vector(const double matrix[3][3], const double vector[3], double result[3]);
void Matrix3Multiply(const double left[3][3], const double right[3][3], double result[3][3]);
int InvertMatrix3(const double matrix[3][3], double inverse[3][3]);
double Norm3Vector(const double vector[3]);
void Normalize3Vector(double vector[3]);
void Matrix6Multiply(const double left[6][6], const double right[6][6], double result[6][6]);
void Matrix6MultiplyTranspose(const double left[6][6], const double right[6][6], double result[6][6]);
void MatrixExponential6(const double matrix[6][6], double exponential[6][6]);

#ifdef __cplusplus
}
#endif

#endif