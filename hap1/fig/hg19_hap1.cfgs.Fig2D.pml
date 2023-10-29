set connect_mode, 1
load hg19_hap1.cfg1.xyz.pdb, cj
load hg19_hap1.cfg1.eps.pdb, ep

hide everything
bg_color white

alter all, vdw = 0.63
color grey, cj
color cyan, cj and resn LAD
show sphere, cj and z<0

color black, ep
set_bond stick_radius, 0.2, ep
show stick, ep
hide everything, ep and chain 0+1 and z>0

set ambient, 0.45
set stick_quality, 3
set sphere_quality, 3
set light_count, 8
set specular, 0.11
set ray_opaque_background, on

reset
set_view (\
     0.965925813,   -0.044943456,    0.254886985,\
     0.000000000,    0.984807730,    0.173648179,\
    -0.258819044,   -0.167731255,    0.951251209,\
     0.000000000,    0.000000000,  -97.158004761,\
     0.007293701,    0.000661850,   -0.945179939,\
    66.073577881,  128.242492676,  -20.000000000 )

