#ifndef CONFIG_H
#define CONFIG_H

typedef struct AcConfig {
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
    double sun_pointing_vector[3]; // Desired sun pointing vector in body frame

    float eclipse_threshold; // Threshold for eclipse detection
    float magnetic_field_min;
} AcConfig_t;

static AcConfig_t ac_config = {
    .detumble_angular_velocity_threshold_low = 0.001, // Example value
    .detumble_angular_velocity_threshold_high = 0.1, // Example value
    
    .kw_detumb = 0.33, // Example value
    .eps_detumb = 0.01, // Example value

    .kw_sunpointing = 0.33, // Example value
    .kp_sunpointing = 0.000025, // Example value
    .eps_sunpointing = 0.01, // Example value

    .kw_nadirpointing = 0.33, // Example value
    .kp_nadirpointing = 0.000025, // Example value
    .eps_nadirpointing = 0.01, // Example value

    .nadir_pointing_vector = {0.0, 0.0, -1.0}, // Pointing towards nadir in body frame
    .sun_pointing_vector = {1.0, 0.0, 0.0}, // Pointing towards sun in body frame

    .eclipse_threshold = 0.5,                           /* Moderate eclipse sensitivity */
    .magnetic_field_min = 1e-6                      /* Minimum magnetic field strength for valid readings */


};

AcConfig_t* GetAcConfig();

#endif