$fn = 100;
rotate_extrude(convexity = 100)
    translate([0.04, 0, 0])
        circle(r = 0.01);
cylinder(h = 0.02, r = 0.02, center = true);
difference(){
cylinder(h = 0.015, r = 0.04, center = true);
translate([0.02, 0, 0])
    cylinder(h = 0.016, r = 0.01, center = true);
translate([-0.01, -0.01732, 0])
    cylinder(h = 0.016, r = 0.01, center = true);
translate([-0.01, 0.01732, 0])
    cylinder(h = 0.016, r = 0.01, center = true);
}
