$fn = 100;
$fs = 0.01; 
$fa = 0.01;
difference(){
translate([0, 0, -0.007])
    cylinder(h = 0.006, r = 0.08, center = true);
translate([0.035, 0.035, 0])
    cylinder(h = 0.021, r = 0.016, center = true);   
translate([0.035, -0.035, 0])
    cylinder(h = 0.021, r = 0.012, center = true);   
translate([-0.035, -0.035, 0])
    cylinder(h = 0.021, r = 0.008, center = true);   
translate([-0.035, 0.035, 0])
    cylinder(h = 0.021, r = 0.004, center = true);   
}
cylinder(h = 0.02, r = 0.03, center = true);
translate([0.07, 0, 0])
    cylinder(h = 0.02, r = 0.01, center = true);
translate([-0.07, 0, 0])
    cylinder(h = 0.02, r = 0.01, center = true);
translate([0, 0.07, 0])
    cylinder(h = 0.02, r = 0.01, center = true);
translate([0, -0.07, 0])
    cylinder(h = 0.02, r = 0.01, center = true);

