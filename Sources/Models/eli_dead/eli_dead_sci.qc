
$modelname "eli_dead_sci.mdl"
$cd ".\"
$cdtexture ".\"
$scale 1.0
$cliptotextures


$bbox 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
$cbox 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
$eyeposition 0.000000 0.000000 60.000000

$bodygroup "torso"
{
	studio "SCI_torso_reference"
	//studio "HEV_torso_reference"
}

$bodygroup "legs"
{
	studio "SCI_legs_reference"
	//studio "HEV_legs_reference"
}

$bodygroup "heads"
{
	studio "SCI_luther_head"
}

// 2 bone controller(s)
$controller 0 "Bip01 Head" XR -60.000000 60.000000
$controller 4 "Mouth" ZR 0.000000 22.000000

// 19 hit box(es)
$hbox 6 "Bip01 L Thigh" 0.000000 -5.690000 -3.550000 18.709999 3.610000 3.430000
$hbox 6 "Bip01 L Calf" 0.000000 -3.770000 -2.970000 16.959999 2.210000 2.410000
$hbox 6 "Bip01 L Foot" -0.070000 -2.290000 -3.140000 3.620000 8.170000 1.450000
$hbox 7 "Bip01 R Thigh" 0.000000 -5.960000 -3.770000 18.639999 3.580000 3.390000
$hbox 7 "Bip01 R Calf" 0.000000 -3.850000 -2.580000 16.969999 2.160000 2.980000
$hbox 7 "Bip01 R Foot" -0.110000 -2.290000 -1.380000 3.610000 8.170000 2.990000
$hbox 3 "Bip01 Pelvis" -6.670000 -4.400000 -6.460000 5.120000 5.170000 6.830000
$hbox 3 "Bip01 Spine" 0.000000 -3.290000 -6.830000 5.320000 5.000000 6.910000
$hbox 2 "Bip01 Spine1" 0.000000 -3.900000 -6.200000 5.700000 5.670000 5.910000
$hbox 2 "Bip01 Spine2" 0.000000 -3.650000 -5.970000 8.990000 4.930000 5.440000
$hbox 4 "Bip01 L UpperArm" 0.000000 -3.170000 -1.870000 5.470000 2.950000 1.580000
$hbox 4 "Bip01 L Forearm" 0.000000 -3.400000 -1.950000 10.950000 2.060000 2.460000
$hbox 4 "Bip01 L UpperArm" 0.000000 -2.880000 -2.280000 11.550000 2.100000 1.700000
$hbox 4 "Bip01 L Hand" 0.000000 -1.640000 -1.960000 6.120000 1.790000 1.770000
$hbox 5 "Bip01 R UpperArm" 0.000000 -3.080000 -1.440000 5.890000 2.970000 2.450000
$hbox 5 "Bip01 R Forearm" 0.000000 -3.350000 -2.710000 10.920000 2.400000 2.110000
$hbox 5 "Bip01 R UpperArm" 0.000000 -2.950000 -1.970000 11.250000 2.220000 2.450000
$hbox 5 "Bip01 R Hand" 0.000000 -1.850000 -1.770000 5.930000 1.590000 2.080000
$hbox 1 "Bip01 Head" -1.000000 -3.960000 -3.860000 8.020000 5.230000 4.070000

// 128 animation sequence(s)
$sequence "lying_on_back" "SCI_lying_on_back" fps 30 
$sequence "lying_on_stomach" "SCI_lying_on_stomach" fps 30 