#include "adcs/state_machine.h"
#define MAX_MODE_ENTRIES 64

AcMode_t mode[MAX_MODE_ENTRIES];
int mode_set[MAX_MODE_ENTRIES];

int DetumbleDone(struct AcType *AC, AcConfig_t *config)
{
    double wmag = sqrt(AC->wbn[0]*AC->wbn[0] + // MODIFICAR cuando se tenga estado del EKF
                      AC->wbn[1]*AC->wbn[1] +
                      AC->wbn[2]*AC->wbn[2]);

    return wmag < config->detumble_angular_velocity_threshold_low;

}

int ExitNominal(struct AcType *AC, AcConfig_t *config)
{
    double wmag = sqrt(AC->wbn[0]*AC->wbn[0] + // MODIFICAR cuando se tenga estado del EKF
                      AC->wbn[1]*AC->wbn[1] +
                      AC->wbn[2]*AC->wbn[2]);

    return wmag < config->detumble_angular_velocity_threshold_high;

}

void UpdateMode(struct AcType *AC, AcConfig_t *config)
{
    static AcMode_t mode[MAX_MODE_ENTRIES];
    static int initialized = 0;
    if (!initialized) {
        for (int i = 0; i < MAX_MODE_ENTRIES; i++) {
            mode[i] = AC_MODE_DETUMBLE;
        }
        initialized = 1;
    }

    enum AcMode currentMode = mode[AC->ID];

    if (AC->ReqMode != currentMode && AC->ReqMode > 0) {
        currentMode = (enum AcMode)AC->ReqMode;
        mode[AC->ID] = currentMode;
        AC->ReqMode = 0;
        return;
    }

    // Normal FSW-controlled mode transitions
    switch(currentMode) {
        case AC_MODE_DETUMBLE:
            if (DetumbleDone(AC, config)) {
                mode[AC->ID] = AC_MODE_SUN_POINTING;
            }
            break;
        case AC_MODE_SUN_POINTING:
            if (ExitNominal(AC, config)) {
                mode[AC->ID] = AC_MODE_DETUMBLE;
            }
            break;
        case AC_MODE_NADIR_POINTING:
            if (ExitNominal(AC, config)) {
                mode[AC->ID] = AC_MODE_DETUMBLE;
            }
            break;
    }






}

enum AcMode GetCurrentMode(struct AcType *AC)
{
    return mode[AC->ID];
}

