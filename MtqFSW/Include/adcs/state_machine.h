#ifndef STATE_MACHINE_H
#define STATE_MACHINE_H

#include "42.h"
#include "adcs/config.h"

typedef enum AcMode {
    AC_MODE_DETUMBLE = 1,
    AC_MODE_SUN_POINTING = 2,
    AC_MODE_NADIR_POINTING = 3
} AcMode_t;

int DetumbleDone(struct AcType *AC, AcConfig_t *config);
int ExitNominal(struct AcType *AC, AcConfig_t *config);
void InitAcMode(struct AcType *AC);

AcMode_t GetCurrentMode(struct AcType *AC);
void UpdateMode(struct AcType *AC, AcConfig_t *config);


#endif 