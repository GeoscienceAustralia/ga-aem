[IMPORT ARCHIVE]
FILEHEADER	0
RECORDFORM	FIXED
DATA	0, 12, NORMAL, 1, 0, , 
CHAN	uniqueid, LONG, NORMAL, 12, 0, Label=uniqueid;DESC=Inversion sequence number;
DATA	12, 6, NORMAL, 1, 0, , 
CHAN	Project_GA, LONG, NORMAL, 6, 0, Label=Project_GA;NULL=-9999;LONGNAME=GA Project Number;
DATA	18, 9, NORMAL, 1, 0, , 
CHAN	Date, LONG, NORMAL, 9, 0, Label=Date;NULL=-9999999;LONGNAME=Date (yyyymmdd);
DATA	27, 4, NORMAL, 1, 0, , 
CHAN	Flight, LONG, NORMAL, 4, 0, Label=Flight;NULL=-99;LONGNAME=Flight Number;
LINENUMBER	31, 10, NORMAL, 1, 0, , 
DATA	41, 8, NORMAL, 1, 0, , 
CHAN	Fiducial, DOUBLE, NORMAL, 8, 1, Label=Fiducial;NULL=-999999.9;LONGNAME=Fiducial Number;
DATA	49, 13, NORMAL, 1, 0, , 
CHAN	Longitude, DOUBLE, NORMAL, 13, 7, Label=Longitude;NULL=-999.9999999;UNITS=deg;LONGNAME=Longitude;extra1=DATUM;extra2=GDA94;
DATA	62, 12, NORMAL, 1, 0, , 
CHAN	Latitude, DOUBLE, NORMAL, 12, 7, Label=Latitude;NULL=-99.9999999;UNITS=deg;LONGNAME=Latitude;extra1=DATUM;extra2=GDA94;
DATA	74, 10, NORMAL, 1, 0, , 
CHAN	Easting, DOUBLE, NORMAL, 10, 2, Label=Easting;NULL=-99999.99;UNITS=m;LONGNAME=Easting;extra1=PROJECTION;extra2=GDA94/MGA54;
DATA	84, 11, NORMAL, 1, 0, , 
CHAN	Northing, DOUBLE, NORMAL, 11, 2, Label=Northing;NULL=-999999.99;UNITS=m;LONGNAME=Northing;extra1=PROJECTION;extra2=GDA94/MGA54;
DATA	95, 8, NORMAL, 1, 0, , 
CHAN	Lidar, DOUBLE, NORMAL, 8, 2, Label=Lidar;NULL=-999.99;UNITS=m;LONGNAME=Final Lidar altimeter;
DATA	103, 8, NORMAL, 1, 0, , 
CHAN	Radalt, DOUBLE, NORMAL, 8, 2, Label=Radalt;NULL=-999.99;UNITS=m;LONGNAME=Final Radar altimeter;
DATA	111, 8, NORMAL, 1, 0, , 
CHAN	DTM, DOUBLE, NORMAL, 8, 2, Label=DTM;NULL=-999.99;UNITS=m;LONGNAME=Final Ground Elevation (AHD);extra1=DATUM;extra2=AHD;
DATA	119, 9, NORMAL, 1, 0, , 
CHAN	input_tx_height, DOUBLE, NORMAL, 9, 2, Label=input_tx_height;UNITS=m;DESC=Input Tx height above ground level;
DATA	128, 9, NORMAL, 1, 0, , 
CHAN	input_tx_roll, DOUBLE, NORMAL, 9, 2, Label=input_tx_roll;UNITS=degrees;DESC=Input Tx roll - left side up + ve;
DATA	137, 9, NORMAL, 1, 0, , 
CHAN	input_tx_pitch, DOUBLE, NORMAL, 9, 2, Label=input_tx_pitch;UNITS=degrees;DESC=Input Tx pitch - nose down + ve;
DATA	146, 9, NORMAL, 1, 0, , 
CHAN	input_tx_yaw, DOUBLE, NORMAL, 9, 2, Label=input_tx_yaw;UNITS=degrees;DESC=Input Tx yaw - turn left + ve;
DATA	155, 9, NORMAL, 1, 0, , 
CHAN	input_txrx_dx, DOUBLE, NORMAL, 9, 2, Label=input_txrx_dx;UNITS=m;DESC=Input Tx - Rx horizonatl inline separation;
DATA	164, 9, NORMAL, 1, 0, , 
CHAN	input_txrx_dy, DOUBLE, NORMAL, 9, 2, Label=input_txrx_dy;UNITS=m;DESC=Input Tx - Rx horizonatl transverse separation;
DATA	173, 9, NORMAL, 1, 0, , 
CHAN	input_txrx_dz, DOUBLE, NORMAL, 9, 2, Label=input_txrx_dz;UNITS=m;DESC=Input Tx - Rx vertical separation;
DATA	182, 9, NORMAL, 1, 0, , 
CHAN	input_rx_roll, DOUBLE, NORMAL, 9, 2, Label=input_rx_roll;UNITS=degrees;DESC=Input Rx roll - left side up + ve;
DATA	191, 9, NORMAL, 1, 0, , 
CHAN	input_rx_pitch, DOUBLE, NORMAL, 9, 2, Label=input_rx_pitch;UNITS=degrees;DESC=Input Rx pitch - nose down + ve;
DATA	200, 9, NORMAL, 1, 0, , 
CHAN	input_rx_yaw, DOUBLE, NORMAL, 9, 2, Label=input_rx_yaw;UNITS=degrees;DESC=Input Rx yaw - turn left + ve;
DATA	209, 4, NORMAL, 1, 0, , 
CHAN	ndata, LONG, NORMAL, 4, 0, Label=ndata;DESC=Number of data in inversion;
DATA	213, 4, NORMAL, 1, 0, , 
CHAN	nlayers, LONG, NORMAL, 4, 0, Label=nlayers;DESC=Number of layers ;
DATA	217, 15, EXP, 1, 0, , 
CHAN	conductivity{30}, DOUBLE, EXP, 15, 6, Label=conductivity;UNITS=S/m;DESC=Layer conductivity;
DATA	667, 9, NORMAL, 1, 0, , 
CHAN	thickness{30}, DOUBLE, NORMAL, 9, 2, Label=thickness;UNITS=m;DESC=Layer thickness;
DATA	937, 9, NORMAL, 1, 0, , 
CHAN	depth_top{30}, DOUBLE, NORMAL, 9, 2, Label=depth_top;UNITS=m;DESC=Depth to top of layer;
DATA	1207, 9, NORMAL, 1, 0, , 
CHAN	depth_bottom{30}, DOUBLE, NORMAL, 9, 2, Label=depth_bottom;UNITS=m;DESC=Depth to bottom of layer;
DATA	1477, 9, NORMAL, 1, 0, , 
CHAN	depth_bottom_negative{30}, DOUBLE, NORMAL, 9, 2, Label=depth_bottom_negative;UNITS=m;DESC=Negative of depth to bottom of layer;
DATA	1747, 9, NORMAL, 1, 0, , 
CHAN	elevation_interface{30}, DOUBLE, NORMAL, 9, 2, Label=elevation_interface;UNITS=m;DESC=Elevation of interface;
DATA	2017, 15, EXP, 1, 0, , 
CHAN	conductivity_sensitivity{30}, DOUBLE, EXP, 15, 6, Label=conductivity_sensitivity;DESC=Conductivity parameter sensitivity;
DATA	2467, 15, EXP, 1, 0, , 
CHAN	conductivity_uncertainty{30}, DOUBLE, EXP, 15, 6, Label=conductivity_uncertainty;UNITS=log10(S/m);DESC=Conductivity parameter uncertainty;
DATA	2917, 15, EXP, 1, 0, , 
CHAN	observed_EMSystem_1_ZS{15}, DOUBLE, EXP, 15, 6, Label=observed_EMSystem_1_ZS;UNITS=fT;DESC=Observed EMSystem 1 Z-component secondary field;
DATA	3142, 15, EXP, 1, 0, , 
CHAN	noise_EMSystem_1_ZS{15}, DOUBLE, EXP, 15, 6, Label=noise_EMSystem_1_ZS;UNITS=fT;DESC=Estimated noise EMSystem 1 Z-component secondary field;
DATA	3367, 15, EXP, 1, 0, , 
CHAN	predicted_EMSystem_1_ZS{15}, DOUBLE, EXP, 15, 6, Label=predicted_EMSystem_1_ZS;UNITS=fT;DESC=Predicted EMSystem 1 Z-component secondary field;
DATA	3592, 15, EXP, 1, 0, , 
CHAN	Alpha_RefCon, DOUBLE, EXP, 15, 6, Label=Alpha_RefCon;DESC=Conductivity reference model constraint alpha parameter;
DATA	3607, 15, EXP, 1, 0, , 
CHAN	Phi_RefCon, DOUBLE, EXP, 15, 6, Label=Phi_RefCon;DESC=Conductivity reference model constraint model norm;
DATA	3622, 15, EXP, 1, 0, , 
CHAN	Alpha_VCsmth, DOUBLE, EXP, 15, 6, Label=Alpha_VCsmth;DESC=Vertical conductivity smoothness constraint alpha parameter;
DATA	3637, 15, EXP, 1, 0, , 
CHAN	Phi_VCsmth, DOUBLE, EXP, 15, 6, Label=Phi_VCsmth;DESC=Vertical conductivity smoothness constraint model norm;
DATA	3652, 15, EXP, 1, 0, , 
CHAN	PhiD, DOUBLE, EXP, 15, 6, Label=PhiD;DESC=Normalised data misfit;
DATA	3667, 15, EXP, 1, 0, , 
CHAN	PhiM, DOUBLE, EXP, 15, 6, Label=PhiM;DESC=Combined model norm;
DATA	3682, 15, EXP, 1, 0, , 
CHAN	Lambda, DOUBLE, EXP, 15, 6, Label=Lambda;DESC=Lambda regularization parameter;
DATA	3697, 4, NORMAL, 1, 0, , 
CHAN	Iterations, LONG, NORMAL, 4, 0, Label=Iterations;DESC=Number of iterations;
