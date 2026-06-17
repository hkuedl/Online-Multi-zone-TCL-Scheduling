# Online Scheduling of Thermostatically Controlled Loads: A Stepwise Reformulation Approach

_This work presents a computationally efficient method for online scheduling of thermostatically controlled loads (TCLs) in multi-zone buildings under uncertain real-time prices. We reformulate the original scheduling problem (multi-stage, uncertain, and coupled) into a single-stage, deterministic, and decoupled problem type, which is solving-efficient and accurate for online applications._

Codes for Submitted Paper "Online Scheduling of Thermostatically Controlled Loads: A Stepwise Reformulation Approach".

Authors: Xueyuan Cui, Liudong Chen, Yi Wang, and Bolun Xu.

## Experiments

To access the required data, please go to [Google Drive](https://drive.google.com/drive/folders/1n2SEkaIN_YUOnrSdE14ptHkDB4ud57J_?usp=drive_link) for direct use. The data of electricity prices is included in ```Data_10days.mat``` (real data) and ```Data_quan_Final_10.mat``` (forecasting data), respectively. The raw data of prices can be collected at [NYISO](https://www.nyiso.com/energy-market-operational-data). The data of RC model parameters and disturbances are included in ```R&C.csv``` and ```90zone_15min.csv```, respectively. 

To reproduce the proposed method, please run the MATLAB codes in ```Codes```: ```Main_proposed.m```, where the sensitivity analysis on the clustering number is also presented. To reproduce the comparisons, please run the MATLAB codes in ```Codes```: ```Main_DE.m```, ```Main_SO.m```, ```Main_GT.m```, ```Main_MMPC.m```.

Furthermore, to reproduce the performance analysis between ASDP and SDDP, please run the MATLAB codes in ```Codes```: ```Com_main.m``` and ```Com_SDDP.m```. We acknowledge the open-source toolbox named [FAST](https://stanford.edu/~lcambier/cgi-bin/fast/tuto.php) to support our experiment on SDDP.

## Citation
```
```
