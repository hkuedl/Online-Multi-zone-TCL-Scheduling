# Efficient_Response_to_Realtime_Price

_This work presents a computationally efficient method for online scheduling of thermostatically controlled loads (TCLs) in multi-zone buildings under uncertain real-time prices. We reformulate the original scheduling problem (multi-stage, uncertain, and decoupled) into a single-stage, deterministic, and decoupled problem type, which is solving-efficient and accurate for online applications._

Codes for Submitted Paper "A Stepwise Reformulation Approach for Online Thermostatically Controlled Loads Scheduling".

Authors: Xueyuan Cui, Liudong Chen, Yi Wang, and Bolun Xu.

## Experiments

To access the required data, please go to ```Data``` folder for direct use. The data of electricity prices is included in ```Data_10days.mat``` (real data) and ```Data_quan_Final_10.mat``` (forecasting data), respectively. The raw data of prices can be collected at [NYISO](https://www.nyiso.com/energy-market-operational-data). The data of RC model parameters and disturbances are included in ```.mat``` 

To reproduce the proposed method, please run the MATLAB codes in ```Codes```: ```Com_main.m```. To reproduce the comparisons, please run the MATLAB codes in ```Codes```: ```Com_main.m```

between the proposed method and SDDP, please run the MATLAB codes in ```Codes```: ```Com_main.m``` and ```matlab Com_SDDP.m```. We acknowledge the open-source toolbox named [FAST](https://stanford.edu/~lcambier/cgi-bin/fast/tuto.php) to support our experiment on SDDP.

To reproduce the experiments of real-time demand response with different comfort functions, please run the MATLAB code in ```Codes```: ```Response.m```.

The figures are generated with the results from MATLAB codes and the Python code in ```Codes```: ```Figures.py```.

## Citation
```
```
