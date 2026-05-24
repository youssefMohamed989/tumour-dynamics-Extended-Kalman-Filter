/*
 * TCS3200 RGB Colour Sensor Readout
 * Instrument-free GPC3 quantification channel
 *
 * Wiring:
 *   S0 -> D4   (frequency scaling)
 *   S1 -> D5   (frequency scaling, set to 20% with S0=HIGH, S1=LOW)
 *   S2 -> D6   (colour filter select)
 *   S3 -> D7   (colour filter select)
 *   OUT -> D8  (frequency output)
 *   OE -> GND  (output enable, active low)
 *   VCC -> 5V
 *
 * Output (Serial, 9600 baud):
 *   R:<value> G:<value> B:<value> GPC3_est:<value_ng_per_mL>
 *
 * Calibration (from manuscript Figure 5E, multi-day replicated):
 *   Red_intensity = -4.54 * [GPC3 ng/mL] + 116.68
 *   Rearranged: [GPC3] = (116.68 - Red_intensity) / 4.54
 *   LOD = 0.31 ng/mL
 *   Linear range: 0.25 - 5 ng/mL
 *
 * Note: Apply fresh single-use calibration against a known standard on
 * each assay day. The equation above is from the replicated multi-day
 * calibration; single-run R2 may vary (see manuscript Figure 5E).
 */

// Pin assignments
const int S0_PIN = 4;
const int S1_PIN = 5;
const int S2_PIN = 6;
const int S3_PIN = 7;
const int OUT_PIN = 8;

// Calibration parameters (multi-day replicated; update from daily standard if required)
const float SLOPE     = -4.54f;   // AU per ng/mL
const float INTERCEPT = 116.68f;  // AU
const float LOD       = 0.31f;    // ng/mL

// Pulse counting window in milliseconds
const unsigned long COUNT_WINDOW_MS = 100UL;

// Number of readings to average per reported value
const int N_AVERAGE = 5;

// Forward declarations
long readChannel(int s2State, int s3State);
float rawToGPC3(float redIntensity);

void setup() {
    Serial.begin(9600);

    pinMode(S0_PIN, OUTPUT);
    pinMode(S1_PIN, OUTPUT);
    pinMode(S2_PIN, OUTPUT);
    pinMode(S3_PIN, OUTPUT);
    pinMode(OUT_PIN, INPUT);

    // Set frequency scaling to 20% (S0=HIGH, S1=LOW)
    digitalWrite(S0_PIN, HIGH);
    digitalWrite(S1_PIN, LOW);

    Serial.println("# TCS3200 GPC3 Readout - initialised");
    Serial.println("# Format: R:<counts> G:<counts> B:<counts> GPC3_est:<ng_per_mL>");
    Serial.println("# LOD = 0.31 ng/mL. Values below LOD reported as <LOD.");
}

void loop() {
    long sumR = 0, sumG = 0, sumB = 0;

    for (int i = 0; i < N_AVERAGE; i++) {
        sumR += readChannel(LOW,  LOW);   // Red filter: S2=LOW, S3=LOW
        sumG += readChannel(HIGH, HIGH);  // Green filter: S2=HIGH, S3=HIGH
        sumB += readChannel(LOW,  HIGH);  // Blue filter: S2=LOW, S3=HIGH
        delay(10);
    }

    float avgR = (float)sumR / N_AVERAGE;
    float avgG = (float)sumG / N_AVERAGE;
    float avgB = (float)sumB / N_AVERAGE;

    float gpc3 = rawToGPC3(avgR);

    Serial.print("R:");
    Serial.print((long)avgR);
    Serial.print(" G:");
    Serial.print((long)avgG);
    Serial.print(" B:");
    Serial.print((long)avgB);
    Serial.print(" GPC3_est:");

    if (gpc3 < LOD) {
        Serial.println("<LOD");
    } else {
        Serial.println(gpc3, 3);
    }

    delay(1000);
}

/*
 * readChannel: Set filter pins and count output pulses over COUNT_WINDOW_MS.
 * Returns pulse count (proportional to light intensity for selected colour).
 */
long readChannel(int s2State, int s3State) {
    digitalWrite(S2_PIN, s2State);
    digitalWrite(S3_PIN, s3State);
    delay(5); // Allow filter to stabilise

    long count = 0;
    unsigned long t_start = millis();
    while ((millis() - t_start) < COUNT_WINDOW_MS) {
        if (pulseIn(OUT_PIN, LOW, 1000UL) > 0) {
            count++;
        }
    }
    return count;
}

/*
 * rawToGPC3: Apply inverse of the calibration equation.
 *   Red_intensity = SLOPE * [GPC3] + INTERCEPT
 *   [GPC3] = (Red_intensity - INTERCEPT) / SLOPE
 * Returns estimated GPC3 in ng/mL. Negative values clamped to 0.
 */
float rawToGPC3(float redIntensity) {
    float gpc3 = (redIntensity - INTERCEPT) / SLOPE;
    if (gpc3 < 0.0f) gpc3 = 0.0f;
    return gpc3;
}
