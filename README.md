# rs_runoff

## Requirements
* Linux or unix environment
* gfortran

## How to run on the wsl?
Firstly, install wsl and a linux distribution (Ubuntu 22.04 is recommended) on windows environment.
Open the terminal, and type

```bash
sudo apt update
sudo apt istall gfortran git
```

Download the source code by,
```bash
git clone https://github.com/ka-yamanoi/rs_runoff.git
```
Then the directory named 'rs_runoff' is generated on the current dirrectory.

After that you can compile the progam by, 
```bash
cd ./rs_runoff
gfortran ./ls_runoff.f90
```
Then, the executive file "a.out" is generated.

After that, type
```bash
./a.out
```
then the simulation programm will be executed.

To set the simulating conditions and/or parameters, edit the following files.
| File name | Description |
| ---- | ---- |
| ./conditions/parameters.dat | Phisical parameters for rainfall and sediment runoff simulation |
| ./conditions/rainfall.dat | Rainfall conditions. |
| ./watershed_data/unitchannels_02.dat | Topographic condition of unit channels|
| ./watershed_data/unit_slopes.dat | Topographic condition of unit slopes |

The output files are followings.
| File name | Description |
| ---- | ---- |
|Q.dat| Temporal change of water discharge|
|Qb.dat| Temporal change of bedload discharge|
|Dzb.dat| Temporal change of river bed deformation|
|Bw.dat| Widths of the unit channels|
|Car.dat| Catchment area for each unit channel downstream ends|
|Cp.dat| Temporal change of gradient of unit channels (radians)|
|H.dat| Temporal change of water depth|
|Zb.dat| Temporal change of river bed elevation|
