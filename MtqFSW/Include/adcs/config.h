#ifndef CONFIG_H
#define CONFIG_H

typedef struct AcConfig
{
    double detumble_angular_velocity_threshold_low;
    double detumble_angular_velocity_threshold_high;

    double kw_detumb;
    double eps_detumb;

    double kw_sunpointing;
    double kp_sunpointing;
    double eps_sunpointing;

    double kw_nadirpointing;
    double kp_nadirpointing;
    double eps_nadirpointing;

    double nadir_pointing_vector[3]; // Desired nadir pointing vector in body frame
    double sun_pointing_vector[3];   // Desired sun pointing vector in body frame

    float eclipse_threshold; // Threshold for eclipse detection
    float magnetic_field_min;
} AcConfig_t;

// Poner todas las configs del sistema de control (parametros, thresholds, etc.) capaz se puede hacer otro tipo de config para el EKF(?)

static AcConfig_t ac_config = {
    .detumble_angular_velocity_threshold_low = 0.01,
    .detumble_angular_velocity_threshold_high = 0.6,

    .kw_detumb = 0.33,
    .eps_detumb = 0.01,

    .kw_sunpointing = 0.33, // Añadir diferencia con ganancias en eclipse
    .kp_sunpointing = 0.000025,
    .eps_sunpointing = 0.01,

    .kw_nadirpointing = 0.33, // Añadir diferencia con ganancias cuando no hay horizonte
    .kp_nadirpointing = 0.000025,
    .eps_nadirpointing = 0.01,

    .nadir_pointing_vector = {0.0, 0.0, -1.0}, // +Z to nadir
    .sun_pointing_vector = {1.0, 0.0, 0.0},    // -Z to sun

    .eclipse_threshold = 0.5,
    .magnetic_field_min = 1e-6

};

typedef struct Ekf_config
{
    double initialState[6];
    double initialCovariance[6][6];
    double processNoise[6][6];
    double minimumProcessNoise[6][6];
    double measurementNoise[3][3];
    double alpha;
}Ekf_config_t;

static Ekf_config_t ekf_config = {
    .initialState = {0},
    .initialCovariance = {0},
    .processNoise = {{1e-8, 0, 0, 0, 0, 0}, {0, 1e-8, 0, 0, 0, 0}, {0, 0, 1e-8, 0, 0, 0}, {0, 0, 0, 1e-6, 0, 0}, {0, 0, 0, 0, 1e-6, 0}, {0, 0, 0, 0, 0, 1e-6}},
    .minimumProcessNoise = {{1e-8, 0, 0, 0, 0, 0}, {0, 1e-8, 0, 0, 0, 0}, {0, 0, 1e-8, 0, 0, 0}, {0, 0, 0, 1e-6, 0, 0}, {0, 0, 0, 0, 1e-6, 0}, {0, 0, 0, 0, 0, 1e-6}},
    .measurementNoise = {{0.01, 0, 0}, {0, 0.01, 0}, {0, 0, 0.01}},
    .alpha = 0.001f
};

AcConfig_t *GetAcConfig();

Ekf_config_t *GetEkfConfig();

#endif