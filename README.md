# rs_runoff

## What is this?

rs_runoff is the source code and example input/output files of rainfall and sediment runoff simulation developed by (Egashira and Matsuki, 2000) and coded by Kazuki Yamanoi for mainly education purposes.
This is simplified and limited version of SiMHiS (Yamanoi and Fujita, 2014;2016)  which excluded a function of sediment production, sediment supply, treatment of non-uniform sediment and suspended sediment transport.
These codes can be used for conducting distributed rainfall runoff simulation, sestimating a sediment transport in the unit-channel system.

## Requirements
* Linux or unix environment

## How to run on the WSL (Windows Subsystem for Linux).
Firstly, install wsl and a linux distribution (Ubuntu 22.04 is recommended) on windows environment.
Open the terminal, and type

```bash
sudo apt update
sudo apt install gfortran git
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

# Input files

To set the simulating conditions and/or parameters, edit the following files.
| File name | Description |
| ---- | ---- |
| ./conditions/parameters.dat | Phisical parameters for rainfall and sediment runoff simulation |
| ./conditions/rainfall.dat | Rainfall conditions. |
| ./watershed_data/unitchannels_02.dat | Topographic condition of unit channels|
| ./watershed_data/unit_slopes.dat | Topographic condition of unit slopes |

## unitchannels_02.dat
The unit channel data looks like
```
7	!	number	of	unit	channels	
cat	next_stream	prev_str01	prev_str02	length(m)	source_elev(m)	outlet_elev(m)
1	-1	2	3	2764.285714	1220	950
2	1	0	0	1028.571429	1520	1220
3	2	4	7	257.1428571	1260	1220
4	3	5	6	1414.285714	1620	1260
5	4	0	0	578.5714286	1860	1620
6	4	0	0	621.4285714	1870	1620
7	3	0	0	2421.428571	2030	1260
```
In table style,
|7|!	number	of	unit	channels| | | | |
| ---- | ---- | ---- | ---- | ----| ---- |
|cat |next_stream	prev_str01|prev_str02|length(m)|source_elev(m)|outlet_elev(m)|
|1|-1|2|3|2764.285714|1220|950|
|2|1|0|0|1028.571429|1520|1220|
|3|2|4|7|257.1428571|1260|1220|
|4|3|5|6|1414.285714|1620|1260|
|5|4|0|0|578.5714286|1860|1620|
|6|4|0|0|621.4285714|1870|1620|
|7|3|0|0|2421.428571|2030|1260|

The frst row indicates the number of unit channels.
The second row is text, do not edit this row.
The column cat indicates ID of the unit channel.
The column next_stream indicates the ID of downstream unit channel.
The columns prev_str01 and prev_str02 indicate the IDs of upstream unit channels, respectively
The column length indicate the length of the unit channels in meter.
The columns source_elev and outlet_elev indicate the elevation at the upstream and downstream end of unit channels, respectively.

## unit_slopes.dat
The unit slope data looks like
```
cat	area(m2)	us_sl(degrees)
1	1103418.367	70.70236774
2	1376632.653	62.33568424
3	588214.2857	59.73447584
4	509234.6939	63.88648774
5	14693.87755	64.53665494
6	147857.1429	49.78253248
7	1265969.388	45.15544123
8	202040.8163	73.74620996
9	247500	53.58944813
10	337040.8163	55.09835624
11	103775.5102	78.69039254
12	179540.8163	70.91399654
13	537704.0816	79.47309275
14	895867.3469	72.01971999
```
The number of unit slopes shuld be twice as the number of unit channels. 

The column cat indicates ID of the unit slope.
The column area and us_sl indicate the area and slopes of unit slope, respectively. 

## parameters.dat
The parameter file looks like,
```
  0.7     !Cmns !Manning number of slopes (m^(-1/3)*s)
  0.05    !Cmn !Manning number of channels (m^(-1/3)*s)
  0.2     !Da !depth of A layer
  0.1     !Db !depth of B layer
  0.001   !Rka !Coefficient of permeability (m/s)
  0.000001 !Rkb !Coefficient for a layer (m/s)
  0.1     ! Hsini !Initial water depth of ground water at downstream end of each slope (m)
  0.0     ! Hcini !Initial water depth of ground water at downstream end of channels (m)
  1.    !Dt !Time step (s)
  10.      !Q0 !Specific water discharge at downstream end for deciding channel width (m^3/s)
  5.       !Parameter for regeme fomula (m^(-1/2)*s^(1/2))
  0.05  !Dm1 !Mean diameter of river bed material (only for uniform bed) (m)
  300.   !data_time !Time interval of output for files (s)
  1800.   !out_time !Time interval of output for display (s)
  -10.    !dzbmin
```

The values at the left indicate the parameter values, and the sentences after "!" indicates the description of each parameters.

## rainfall.dat
The rainfall file looks like
```
24 !number of data
3600.d0 !time interval
T(sec) prec(mm)
3600 0.
7200 0.
10800 0.
14400 10.
18000 20.
21600 40.
25200 80.
28800 20.
32400 10.
36000 0.
39600 0.
43200 0.
46800 0.
50400 0.
54000 0.
57600 0.
61200 0.
64800 0.
68400 0.
72000 0.
75600 0.
79200 0.
82800 0.
86400 0.
```
The value in the first row indicates the number of data (n). 
The value in the second row indicates the time interval (sec.) of each values in the data (Train).
The third row is the text row, which is not used for the simulation.
The forth to (N+3)th rows are the rainfall data, first column is time index and second column is rainfall intensity in mm/h.

The first rainfall intensity is used for the calculation from T=0 to T=Train. 
In the case of this sample file, the value is used from T=0 to T=3600 [sec.).

## Output files

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

