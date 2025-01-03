---

# mimo-4-dual-polarized-antennas

This repository "Quaternion-Based MIMO Systems and BER Simulations" contains scripts and simulations for Bit Error Rate (BER) evaluations in various communication systems, including MIMO configurations, leveraging quaternion mathematics and advanced modulation schemes.

---

## Directory: **SPL Paper**

### Files and Features:
1. **Combine_constellation**:
   - Simulates **100k random bits**.
   - Evaluates two case scenarios:
     - **QPSK + QPSK** 1x1 system.
     - **QPSK + 16-PSK** 1x1 system.
   - BER evaluation performed at **18 SNR bins** ranging from '-4' to '30' dB with a step size of '2'.

2. **SPL_attempt**:
   - Simulates **100k random bits**.
   - BER evaluation conducted at **18 SNR bins** ranging from '-4' to '30' dB with a step size of '2'.

---

## Directory: **1x1 QODs with phase_rot**

### Description:
- Simulates BER for a **1x1 MIMO system**.
- Utilizes quaternion mathematics under a **rotated QPSK modulation scheme**.

---

## Directory: **2x1--RR**

### Description:
- Simulates BER for a **2x1 MIMO system**.
- Implements quaternion mathematics with a **rotated QPSK modulation scheme**.

---

## Directory: **QODs 2x1**

### Description:
- Conducts a complete simulation and analysis of a **2x1 system**.
- Employs quaternion mathematics to model the system.
- Utilizes a **tweaked QPSK constellation scheme** for modulation.

---

## Directory: **Corrected**

### Description:
- Contains the **combined overall analysis**.
- Features improved simulation methodologies for enhanced accuracy.

---

## Getting Started

### Prerequisites:
- MATLAB or Python (depending on the implementation).
- Basic understanding of MIMO systems, quaternion mathematics, and modulation schemes.

### Usage:
1. Clone the repository:  
   '''bash
   git clone https://github.com/SaadIqbalEE/mimo-4-dual-polarized-antennas.git
   cd mimo-4-dual-polarized-antennas
   '''
2. Navigate to the desired directory and run the simulation files.

---

## Contributing
Contributions are welcome! If you find bugs or have suggestions for improvement, please open an issue or submit a pull request.

---

## License
This project is licensed under the MIT License. See the 'LICENSE' file for more details.

---

