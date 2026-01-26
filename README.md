# Efficient_Response_to_Realtime_Price

_This work presents a computationally efficient method for online scheduling of thermostatically controlled loads (TCLs) in multi-zone buildings under uncertain real-time prices. We reformulate the original scheduling problem (multi-stage, uncertain, and decoupled) into a single-stage, deterministic, and decoupled problem type, which is solving-efficient and accurate for online applications._

Codes for Submitted Paper "A Stepwise Reformulation Approach for Online Thermostatically Controlled Loads Scheduling".

Authors: Xueyuan Cui, Liudong Chen, Yi Wang, and Bolun Xu.

## Experiments

To access the required data, please go to ```Data``` folder for direct use. The data of electricity prices is included in ```Data_10days.mat``` (real data) and ```Data_quan_Final_10.mat``` (forecasting data), respectively. The raw data of prices can be collected at [NYISO](https://www.nyiso.com/energy-market-operational-data). The data of RC model parameters and disturbances are included in ```R&C.csv``` and ```90zone_15min.csv```, respectively. 

To reproduce the proposed method, please run the MATLAB codes in ```Codes```: ```Main_proposed.m```, where the sensitivity analysis on the clustering number is also presented. To reproduce the comparisons, please run the MATLAB codes in ```Codes```: ```Main_DE.m```, ```Main_SO.m```, ```Main_GT.m```.

Furthermore, to reproduce the performance analysis between ASDP and SDDP, please run the MATLAB codes in ```Codes```: ```Com_main.m``` and ```Com_SDDP.m```. We acknowledge the open-source toolbox named [FAST](https://stanford.edu/~lcambier/cgi-bin/fast/tuto.php) to support our experiment on SDDP.

## Citation
```
```
