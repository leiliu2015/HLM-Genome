set connect_mode, 1
load hg19_hap1.cfg1.xyz.pdb, cj
load hg19_hap1.cfg1.eps.pdb, ep

hide everything
bg_color white

set ribbon_trace_atoms, 1
set ribbon_sampling, 10
set ribbon_width, 6.0
show ribbon, cj
spectrum chain, carbon_cyan_lightmagenta_yellow_salmon_hydrogen_slate_orange, cj

color black, ep
set_bond stick_radius, 0.2, ep
show stick, ep

set ambient, 0.45
set stick_quality, 3
set sphere_quality, 3
set light_count, 8
set specular, 0.11
set ray_opaque_background, on

reset
turn y, 180

