#ifndef MTQFSW_ADCS_MAGSUN_H
#define MTQFSW_ADCS_MAGSUN_H

#include "42.h"

enum MtqMode {
	MTQ_MODE_DETUMBLE = 0,
	MTQ_MODE_SUN_POINTING = 1,
	MTQ_MODE_NADIR_POINTING = 2
};

int adcsPropatMagController(struct AcType *AC);
int adcsMagSunUBA(struct AcType *AC);
int adcsRwTriadTLEUBA(struct AcType *AC);

int adcsUBA(struct AcType *AC);


#endif
