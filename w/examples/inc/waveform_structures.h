// waveform complete structure definition
struct wav_raw {
  int weeksow;                  // 0
  short millisecond;            // 4
  char status_numcorr;          // 6
  char link_updw;               // 7
  char prn;                     // 8
  char max_pos;                 // 9
  char amplitude;               // 10
  char phase;                   // 11
  int range_model;              // 12
  int doppler_avion;            // 16
  int sampling_freq_int;        // 20
  short sampling_freq_frac;     // 24
  short d_freq;                 // 26
  short d_tao;                  // 28
  short sin_elevation;          // 30
  signed char data[128];    // 32
};


// waveform complete structure definition
struct wav_int {
  double              inertialR[3];        // Attitude parameters
  double              positionR[3];        // Position Receiver
  double              positionS[3];        // Position Spec Point
  double              positionT[3];        // Transmitter position
  double              velocityR[3];        // Position Velocity
  double              velocityT[3];        // Transmitter velocity
  double              ant_beam[3];         // Orientation of the antenna over the TSR system (Transmitter-Specular_point-Receiver)
  double              AltDelayData;
  double              ScattDelayData;
  double              AltDelayModel;
  double              ScattDelayModel;
  double              GeomDelayModel;
  double              AltDelayResidual;
  float               Anisofact;                // Anisotropy factor v mss
  float               AtmosphericDelay;         // Atmospheric delay
  float               AircInertialDelay;        // Delay due to the different position of the antennas
  float               Azimuth_T;
  float               Elevation_T;
  float               Height_ellipsoidal_T;
  float               Height_ellipsoidal_R;
  float               Height_local_T;
  float               Height_local_R;
  float               InstrumentalDelay;        // Constant Instrumental delay
  float               Latitude_S;
  float               Longitude_S;
  float               Latitude_R;
  float               Longitude_R;
  float               Latitude_T;
  float               Longitude_T;
  float               Mss[2];                   // mean squares slopes
  float               Omega;                    // Fetch parameter
  float               Peak_Sensitivity;         // delay with respect to mss
  float               Spec_Sensitivity;         // delay with respect to mss
  float               s;                        // salinity
  float               T;                        // Temperature
  float               undulation;               // undulation
  float               WavAmpMax;                // SNRv
  float               WavPeakPos;               // Peak delay wrtith SP delay
  float               WavSpecPos;               // Specular delay
//  float               SigmaPeakPos;
//  float               SigmaSpecPos;
  float               WindDirection;            // Wind Azimuth from North
  float               WindMagnitude;            // Magnitude  Wind velocity
  float               MssEstimate;
  float               DeltaHEstimate;
  float               SigmaMssEstimate;
  float               SigmaDeltaHEstimate;
  float               PostFitResidual;
  float               BiasEstimate;
  float               data[64];
  float               chisq;                    // Normalized Chi2 of the Wav adj.
  int                 doppler_avion;
  int                 range_model;
  int                 sampling_freq_int;
  int                 weeksow;
  int                 num_coherent;
  int                 num_uncoherent;
  short               count_waves;
  short               d_freq;
  short               d_tao;
  short               millisecond;
  char                model_mss;
  char                numcorr;
  char                link_updw;
  char                Polarization;              // Polarization L or R
  char                prn;
  char                complete;
  char                valid;
  char                dummy[1];                 //used to set the size of the struct to a multiple of 8
};
