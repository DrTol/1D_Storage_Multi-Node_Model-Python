[![ko-fi](https://ko-fi.com/img/githubbutton_sm.svg)](https://ko-fi.com/F2F01JB1KE)

# 1D_Storage_Multi-Node_Model-Python
A Python implementation of the one-dimensional multi-node method for simulating stratified sensible thermal storage tanks. The model resolves the axial temperature distribution along the tank height using a multi-node approach, enabling detailed analysis of charging, discharging, idle operation, and stratification effects.

This repository provides both the original reference script (`MultiNodeModel.py`) and the new improved implementation (`storage1d.py` & `ExampleScript.py`).

## Table of Contents
- [Description](README.md#Description)
- [Previous Limitations](README.md#Previous-Limitations)
- [Recent Changes](README.md#Recent-Changes)
- [Example Script](README.md#Example-Script)
- [License](README.md#License)
- [Acknowledgement](README.md#Acknowledgement)
- [How2Cite](README.md#How2Cite)
- [References](README.md#References)

## Description
Thermal energy storage is a key component in modern district heating (DH), renewable energy systems, and building energy applications. The multi-node method is a widely used 1D stratified tank modeling approach:

1. The tank is discretized into several vertical nodes (layers).
2. Heat and mass balance equations are applied per node.
3. Effects included:
* Charging (hot fluid injection at a chosen inlet node).
* Discharging (fluid withdrawal at a chosen outlet).
* Idle operation with thermal losses and axial conduction.
* Stratification behavior due to buoyancy and flow direction.

## Previous Limitations
The file `MultiNodeModel.py` (kept for reference) was the first implementation with limitations:
* Use of global variables (not modular).
* Susceptibility to NumPy deprecation warnings.
* Limited operation modes (charging/discharging separately defined).
* No execution-time reporting.
* Simplistic input/output interface.

For new projects, please use the improved scripts (`storage1d.py` + `ExampleScript.py`).

## Recent Changes

The new version (2025 update) introduces:
1. Modular structure (storage1d.py) with a clean function simulate_stratified_tank().
2. Dataclass-based configuration (TankConfig) for clarity and validation.
3. Single signed flow input (m_flow): (i) Positive = charging (ii) Negative = discharging and (iii) Zero = idle
4. Idle mode fully supported.
5. Execution time reported in outputs.
6. Input validation and robust error messages.
7. No global variables; function is self-contained.

## Example Script
The file [ExampleScript.py](https://github.com/DrTol/1D_Storage_Multi-Node_Model-Python/blob/main/ExampleScript.py) demonstrates how to use the model in four realistic scenarios:
1. Full charging — start with a cold tank and charge with hot water from the bottom.
2. Discharging — start from a charged state and discharge from the top outlet.
3. Idle — start with a stratified profile and simulate losses and axial conduction.
4. Combined operation — charge → idle → discharge with varying flow rates.

The script produces:
1. Top/bottom temperature curves
2. Temperature field plots (time vs. height)
3. Mass flow schedule plots
4. Vertical temperature profiles vs. tank height at selected times

This makes it suitable for research, teaching, and system evaluation.

Please see [MultiNodeModel.py](https://github.com/DrTol/1D_Storage_Multi-Node_Model-Python/blob/main/MultiNodeModel.py) for the older version. 

## License
You are free to use, modify and distribute the code as long as **authorship is properly acknowledged**.

## Acknowledgement
Above all, I give thanks to **Allah, The Creator (C.C.)**, and honor His name **Al-‘Alīm (The All-Knowing)**.

This repository is lovingly **dedicated to my parents** who have passed away, in remembrance of their guidance and support.

I would also like to thank ChatGPT (by OpenAI).

We would like to acknowledge all of the open-source minds in general for their willing of share (as apps or comments/answers in forums), which has encouraged our department to publish our tools developed.

## How2Cite
1. Tol, Hİ. pressure_loss_calculator-Python. DOI: 10.5281/zenodo.4167751. GitHub Repository 2020; https://github.com/DrTol/1D_Storage_Multi-Node_Model-Python

## References
- Unrau, Cody. Numerical investigation of one-dimensional storage tank models and the development of analytical modelling techniques. M.Sc. Thesis. McMaster University, Hamilton, Canada. 
- Kleinbach, Eberhard Markus. Performance study of one-dimensional models for stratified thermal storage tank. M.Sc. Thesis. University of Wisconsin-Madison, Madison, Wisconsin.
