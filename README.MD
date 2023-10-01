# Underwater Backscatter Link Budget Tool

Acoustic backscatter is a recent communication technology that enables ultra-low power sensing and connectivity in underwater environments. This repository includes an interactive graphical tool (Matlab App) for quickly evaluating the performance of acoustic backscatter in real-world scenarios. The tool enables users to customize a backscatter system with practical underwater transducers and channel properties, then uses an analytical framework to predict the system's power-up and communication range. Researchers and IoT designers with no acoustics background can leverage our tool to make informed design choices and optimize their system's performance. 

A demo video showing how to use an early version of the tool can be found here: https://youtu.be/CJBYtOHUjC8

## Installing the Tool

### Installer for Windows
You can download and install a standalone executable of the tool as desktop app. The latest version of the installer can be found in the [latest release](https://github.com/signalkinetics/linkbudget/releases).
### Other
If you have MATLAB 2023a or newer you can clone this repository and run the tool by opening Backscatter.mlapp from MATLAB.

## Using the Tool
The GUI is arranged into panels to configure the different components and parameters of an underwater acoustic backscatter system. The panels are arranged from left-to-right and top-to-bottom as a typical workflow for setting up your system.

### Physical Setup Panel
You start by selecting and configuring the underwater transducers to be used in your system, mainly the transmitter (Tx), receiver (Rx), and backscatter node. 
For each component, you can click the dropdown menu to select a transducer of your choice:
1. **Spherical**: Assumes an omnidirectional transducer with a uniform gain (G=0) across all frequencies. Select this option if you want to model an ideal hydrophone or thin-shell transducer.
1. **Custom**: Select this option to define uniform properties for your transducer across frequencies. Upon selecting custom, you can click on the "Config." button that appears next to the dropdown menu to input parameters of your transducer.
	You can either specify the:
	1. Transducer's **gain** (G) in dB.
	1. Transducer's **directivity index** (DI) in dB and **efficiency** in %.
	1. Transducer's **transmit voltage response** (TVR) in dB. Note that you have to specify the transducers complex electrical impedance (Z=R+jX) to use TVR in the tool.
1. A commercial transducer: The tool includes the frequency dependent paramters of three popular underwater communication transducers (B-Tech BT-2RCL, Benthowave BII-7511, and Geospectrum M27-931) from their datasheets. Select any of these transducers and click on "View" to see its gain, TVR and impedance versus frequency.

You can also import measured or datasheet parameters of other transducers using csv files. To import additional transducers into the tool, save their properties as a CSV file and place that file into the "Transducers" folder of the tool. The "Transducers" folder is under "{INSTALL_DIR}\application" for the standalone version of the tool.
The data in the CSV file should be formated as columns with following format: 
| Frequency [Hz] |	TVR [dB] |	 R [ohm] |	 X [ohm] |
| --- | --- | --- | --- |
| 1000 | 150 | 0.1 | 50 |
| 1100 | 152 | 0.5 | 53 | 
| ... | ... | ... | ... |

The first row is reserved for headers.

After selecting the desired transducer for Tx, Rx, and the backscatter node. You can configure the parameter of the underwater channel (acoustic propagation medium) using the "Config. Channel" button.
You can enter a numerical value for the following parameters for the medium properties:
1. The **spreading factor** *k* which defines how the acoustic waves spread in the channel. Use k=2 for spherical spreading (ex. a vertical channel in the ocean (i.e. surface to ocean floor)), k=1 for cylindrical spreading (ex. a shallow river with hard bottom), and k=1.5 for practical spreading (ex. shallow water, underwater channels, etc.).
1. Water **density** in kg/m^3 (assumed uniform throughout the channel).
1. **Speed of sound** in m/s (assumed uniform throughout the channel).
1. Average **depth** of the channel in meters. Used in calculating the attenuation (absorption) coefficient of acoustic waves underwater.
1. Water **temperature** in degrees Celsius. The temperature is assumed uniform throughout the channel and used to calculte acoustic attenuation.
1. Water **salinity** in parts per thousand (ppt). A default value of 35 is used for the average salinity of the ocean. Use a value lower than 0.5 for freshwaters in rivers. The water salinity has a strong effect on its acoustic attenuation especially at low frequencies (below 200 kHz).
1. The **pH** value of the water indicating its acidity which is used in calculating attenuation. A default pH of 8 is used for average ocean acidity (alkaline).  

The underwater ambient noise strongly influences the acheivable SNR in the underwater channel. You can define the following parameters for approximating the noise at Rx:
1. **Wind speed** at the sea surface in m/s.
1. **Shipping** activity factor that approximates the shipping noise.

### Analysis Panel
Use this panel to specify the backscatter parameters of interest for plotting. You can choose to plot one of the following parameters on the y-axis:
1. **Harv Power/SNR**: Plots the harvested power at the backscatter node in the downlink plot and the acheivable SNR at the receiver (Rx) in the uplink plot. 
1. **Max Range**: Plots the maximum operating range (distance) for powering up the backscatter node in the downlink plot and the maximum operating range for backscatter communication in the uplink plot. Note that the maximum range for communication and power are usually different and they depend on the power and SNR thresholds set in the "Data & Power Settings" panel and the secondary parameters in the "Analysis Settings" panel.
1. **Min P_in**: Plots the minimum input power at the transmitter that is need for power-up (downlink plot) and communication (uplink plot).

You can also choose one of the following parameters to have as the dependent parameter varied on the x-axis. The values of the other parameters are kept constant, and they are determined from the "Analysis Settings" panel.
1. **Frequency**: The carrier frequency transmitted from Tx.
1. **Range**: The distance between Tx/Rx and the backscatter node.
1. **P_in**: The input electrical power to Tx.

Note that the units of the dependent parameter is determined from the "Analysis Settings" panel.

### Analysis Settings Panel
This panel defines the secondary operation parameters used for the analysis. The analyzed (primary) operation parameters are selected in the "Analysis Panel" and are plotted on the x-axis of the downlink/uplink plots. The secondary parameters are used to evaluate the plots.

### Data & Power Settings Panel
Use this panel to define:
1. The **data bandwidth** (communication rate) in (Hz/kHz/MHz). The bandwidth is the rate (frequency) of switching the backscatter node between the reflection and absorption state to transmit uplink data. This parameter only affects the (communication range/uplink) plot.
1. The **power threshold** in (W,dBW,dBm). The power threshold defines the minimum required power for operating the backscatter node electronics (microcontroller, FPGA, sensors, etc). The threshold is shown as a dashed red line in the downlink plot when Harv power/SNR are plotted. The power threshold also determines the maximum power-up range and the minimum input power for powering up the backscatter electronics at a given range.
1. The **SNR threshold** for successful communication in dB. The SNR threshold is shown as a dashed red line in the uplink plot when Harv power/SNR are plotted. The SNR threshold also determines the maximum communication range and the minimum input power for successful backscatter communication at a given range.

## Known Limitations
1. The theoeretical model assumes a uniform underwater acoustic channel between Tx/Rx and the backscatter node. Other channel effects such as thermoclines and multipath are not currently supported.
1. The theoretical model assumes the far-field approximation for acoustic propagation. Non-physical results might be generated for ranges smaller than the sum of the Rayleigh distance of Tx/Rx and the backscatter node. 
1. Perfect electrical impedance matching is assumed between the backscatter node and the connected electronics. The effect of electrical impedance mismatch might be included in future versions.


## Cite This Work

**Authors**:  [Waleed Akbar](https://signal-kinetics.media.mit.edu/people/waleed-akbar/), [Ahmed Allam](https://ahmed-allam.com/), [Fadel Adib](http://www.mit.edu/~fadel/)

**Paper**: https://dl.acm.org/doi/10.1145/3603269.3610836

Cite the paper using BibTex as follows:

```
@inproceedings{10.1145/3603269.3610836,
author = {Allam, Ahmed and Akbar, Waleed and Adib, Fadel},
title = {Demo: Underwater Backscatter Link Budget Tool},
year = {2023},
isbn = {9798400702365},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
url = {https://doi.org/10.1145/3603269.3610836},
doi = {10.1145/3603269.3610836},
abstract = {Acoustic backscatter is a recent communication technology that enables ultra-low power sensing and connectivity in underwater environments. This demo presents an interactive graphical tool (Matlab App) for quickly evaluating the performance of acoustic backscatter in real-world scenarios. Our tool enables users to customize a backscatter system with practical underwater transducers and channel properties, then uses an analytical framework to predict the system's power-up and communication range. Researchers and IoT designers with no acoustics background can leverage our tool to make informed design choices and optimize their system's performance. A demo video of the tool can be found here: https://youtu.be/CJBYtOHUjC8},
booktitle = {Proceedings of the ACM SIGCOMM 2023 Conference},
pages = {1191–1192},
numpages = {2},
keywords = {underwater backscatter networking, link budget analysis, subsea IoT},
location = {New York, NY, USA},
series = {ACM SIGCOMM '23}
}

  
```