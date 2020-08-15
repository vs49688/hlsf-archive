
$modelname "eli_dead_hev.mdl"
$cd ".\"
$cdtexture ".\"
$scale 1.0
$cliptotextures


$bbox 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
$cbox 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
$eyeposition 0.000000 0.000000 60.000000

$bodygroup "torso"
{
	//studio "SCI_torso_reference"
	studio "HEV_torso_reference"
}

$bodygroup "legs"
{
	//studio "SCI_legs_reference"
	studio "HEV_legs_reference"
}

$bodygroup "heads"
{
	studio "HEV_luther_head"
}

// 2 bone controller(s)
$controller 0 "Bip01 Head" XR -60.000000 60.000000
$controller 4 "Mouth" ZR 0.000000 22.000000

// 18 hit box(es)
$hbox 3 "Bip01 Pelvis" -9.890000 -5.660000 -8.600000 2.940000 6.520000 6.790000
$hbox 6 "Bip01 L Thigh" 0.000000 -5.920000 -4.000000 19.370001 3.780000 3.730000
$hbox 6 "Bip01 L Calf" 0.000000 -5.690000 -3.620000 18.830000 3.070000 2.940000
$hbox 6 "Bip01 L Foot" 0.000000 -2.990000 -4.040000 6.420000 8.680000 1.610000
$hbox 7 "Bip01 R Thigh" 0.000000 -5.930000 -3.870000 19.270000 3.750000 3.790000
$hbox 7 "Bip01 R Calf" -0.330000 -5.730000 -2.980000 18.950001 2.960000 3.190000
$hbox 7 "Bip01 R Foot" 0.000000 -3.240000 -2.320000 6.210000 8.440000 3.280000
$hbox 3 "Bip01 Spine" 0.000000 -4.890000 -6.970000 8.370000 8.000000 6.730000
$hbox 2 "Bip01 Spine1" -1.000000 -5.990000 -8.220000 5.960000 8.000000 8.480000
$hbox 4 "Bip01 L Clavicle" -0.400000 -4.180000 0.000000 7.200000 4.750000 4.200000
$hbox 4 "Bip01 L UpperArm" -2.520000 -4.570000 -3.920000 9.860000 2.280000 3.430000
$hbox 4 "Bip01 L Forearm" 0.000000 -2.850000 -2.080000 11.900000 2.850000 2.750000
$hbox 4 "Bip01 L Hand" 0.000000 -1.320000 -2.290000 4.660000 2.600000 1.650000
$hbox 5 "Bip01 R Clavicle" 0.000000 -4.150000 -4.190000 7.530000 4.120000 4.240000
$hbox 5 "Bip01 R UpperArm" -2.690000 -4.650000 -3.420000 9.470000 2.280000 3.760000
$hbox 5 "Bip01 R Forearm" 0.000000 -2.730000 -2.560000 11.000000 3.140000 2.190000
$hbox 5 "Bip01 R Hand" 0.000000 -1.320000 -2.910000 10.380000 3.230000 4.000000
$hbox 1 "Bip01 Head" -0.660000 -4.250000 -5.240000 11.190000 7.970000 4.630000

// 128 animation sequence(s)
$sequence "lying_on_back" "HEV_lying_on_back" fps 30 
$sequence "lying_on_stomach" "HEV_lying_on_stomach" fps 30 
$sequence "lying_on_side" "HEV_lying_on_side" fps 30 