#ifndef CONFIG_H
#define CONFIG_H

typedef struct AcConfig
{
    double detumble_w_thld_low;
    double detumble_w_thld_high;

    //Physics constants
    double J[3][3];
    double inv_J[3][3];

    //Gains
    double kw_detumb;
    double eps_detumb;

    double kw_sunpointing;
    double kp_sunpointing;
    double eps_sunpointing;

    double kw_sunpointing_eclipse;
    double eps_sunpointing_eclipse;

    double kw_nadirpointing;
    double kp_nadirpointing;
    double eps_nadirpointing;

    double kw_nadirpointing_nohorizon;
    double kp_nadirpointing_nohorizon;
    double eps_nadirpointing_nohorizon;

    // References
    double nadir_pointing_vector[3]; // Desired nadir pointing vector in body frame
    double sun_pointing_vector[3];   // Desired sun pointing vector in body frame

    // Popurri
    float eclipse_threshold; // Threshold for eclipse detection
    float magnetic_field_min;
    int rmm_compensate;

    float max_dipole;
    double b2mtq_matrix[3][3];
    double mtq2b_matrix[3][3];

    // EHS
    double ehs2b_matrix[3][3];

    //Mag
    double mag_bias[3];
    double mag2b_matrix[3][3];

    //Gyro
    double gyro_bias[3];
    double gyro2b_matrix[3][3];

} AcConfig_t;

// Poner todas las configs del sistema de control (parametros, thresholds, etc.) capaz se puede hacer otro tipo de config para el EKF(?)

static AcConfig_t ac_config = {
    .detumble_w_thld_low = 0.005,
    .detumble_w_thld_high = 0.1,

    .kw_detumb = 0.3,
    .eps_detumb = 0.01,

    .J = {{0.0021, -0.0001, 0.0001}, {-0.0001, 0.0018, -0.0001}, {0.0001, -0.0001, 0.002}},

    .inv_J = {{478.475,  25.3232, -22.6576}, { 25.3232, 558.443,  26.656}, {-22.6576,  26.656, 502.466}},

    .kw_sunpointing = 0.3, // Añadir diferencia con ganancias en eclipse
    .kp_sunpointing = 0.000025,
    .eps_sunpointing = 0.01,

    .kw_sunpointing_eclipse = 0.3, // Añadir diferencia con ganancias en eclipse
    .eps_sunpointing_eclipse = 0.04,

    .kw_nadirpointing = 0.3, // Añadir diferencia con ganancias cuando no hay horizonte
    .kp_nadirpointing = 0.000025,
    .eps_nadirpointing = 0.01,

    .kw_nadirpointing_nohorizon= 0.3,
    .kp_nadirpointing_nohorizon= 0.000025,
    .eps_nadirpointing_nohorizon= 0.03,

    .nadir_pointing_vector = {0.0, 0.0, 1.0}, // +Z to nadir
    .sun_pointing_vector = {0.0, 0.0, -1.0},    // -Z to sun

    .eclipse_threshold = 0.2,
    .magnetic_field_min = 1e-6,
    .rmm_compensate = 1,

    .max_dipole = 0.2,
    .b2mtq_matrix = {{1.0, 0.0, 0.0},
                    {0.0, 1.0, 0.0},
                    {0.0, 0.0, 1.0}},
    .mtq2b_matrix = {{1.0, 0.0, 0.0},
                    {0.0, 1.0, 0.0},
                    {0.0, 0.0, 1.0}},

    .ehs2b_matrix = {{1.0, 0.0, 0.0},
                  {0.0, 0.3746, -0.9271}, // basada en 68°
                  {0.0, 0.9271, 0.3746}},

    .mag_bias = {{0.0}, {0.0}, {0.0}},
    .mag2b_matrix = {{1.0, 0.0, 0.0},
                    {0.0, 1.0, 0.0},
                    {0.0, 0.0, 1.0}},

    .gyro_bias = {{0.0}, {0.0}, {0.0}},
    .gyro2b_matrix = {{1.0, 0.0, 0.0},
                    {0.0, 1.0, 0.0},
                    {0.0, 0.0, 1.0}},
};


typedef struct Ekf_config
{
    double initialState[6];
    double initialCovariance[6][6];
    double processNoise[6][6];
    double minimumProcessNoise[6][6];
    double measurementNoise[3][3];
    double minimumMeasNoise[3][3];
    double J[3][3];
    double alpha;
}Ekf_config_t;

static Ekf_config_t ekf_config = {
    .initialState = {},
    .initialCovariance = {{2e-7, 0, 0, 0, 0, 0}, {0, 2e-7, 0, 0, 0, 0}, {0, 0, 2e-7, 0, 0, 0}, {0, 0, 0, 2e-3, 0, 0}, {0, 0, 0, 0, 2e-3, 0}, {0, 0, 0, 0, 0, 2e-3}}, //[[eye(3)*0.00001^2 zeros(3)];[zeros(3) eye(3)*0.001^2]]*2000
    .processNoise = {{1e-7, 0, 0, 0, 0, 0}, {0, 1e-7, 0, 0, 0, 0}, {0, 0, 1e-7, 0, 0, 0}, {0, 0, 0, 1e-5, 0, 0}, {0, 0, 0, 0, 1e-5, 0}, {0, 0, 0, 0, 0, 1e-5}},
    .minimumProcessNoise = {{1e-7, 0, 0, 0, 0, 0}, {0, 1e-7, 0, 0, 0, 0}, {0, 0, 1e-7, 0, 0, 0}, {0, 0, 0, 1e-5, 0, 0}, {0, 0, 0, 0, 1e-5, 0}, {0, 0, 0, 0, 0, 1e-5}},
    .measurementNoise = {{0.01, 0, 0}, {0, 0.01, 0}, {0, 0, 0.01}},
    .minimumMeasNoise = {{0.01, 0, 0}, {0, 0.01, 0}, {0, 0, 0.01}},
    .alpha = 0.01f
};

AcConfig_t *GetAcConfig();

Ekf_config_t *GetEkfConfig();

#endif