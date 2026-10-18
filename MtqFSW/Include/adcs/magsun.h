#ifndef MTQFSW_ADCS_MAGSUN_H
#define MTQFSW_ADCS_MAGSUN_H

#include "42.h"

int adcsPropatMagController(struct AcType *AC);
int adcsMagSunUBA(struct AcType *AC);
int adcsRwTriadTLEUBA(struct AcType *AC);

int adcsUBA(struct AcType *AC);
void cssProcessing(struct AcType *AC);
void esProcessing(struct AcType *AC);


#endif
