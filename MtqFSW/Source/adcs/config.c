#include "adcs/config.h" 

AcConfig_t *GetAcConfig() {
    return &ac_config;
}

Ekf_config_t *GetEkfConfig() {
    return &ekf_config;
}