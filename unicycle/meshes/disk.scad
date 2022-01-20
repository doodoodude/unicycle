$fn = 100;
$fs = 0.01; 
$fa = 0.01;
translate([0, 0, -0.007])
    cylinder(h = 0.006, r = 0.08, center = true);
translate([0, 0, -0.007])
    cylinder(h = 0.008, r = 0.045, center = true);
cylinder(h = 0.02, r = 0.03, center = true);
translate([0.064, 0, 0])
    cylinder(h = 0.02, r = 0.016, center = true);
translate([-0.064, 0, 0])
    cylinder(h = 0.02, r = 0.016, center = true);
translate([0, 0.064, 0])
    cylinder(h = 0.02, r = 0.016, center = true);
translate([0, -0.064, 0])
    cylinder(h = 0.02, r = 0.016, center = true);

