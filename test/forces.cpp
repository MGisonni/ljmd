#include <gtest/gtest.h>
extern "C" { // since our lib is C code, but the test here is C++
#include <mdlib-util.h>
#include <mdlib.h>
}

//test forces.c for the function force, which is used to compute the forces on each atom of a system of 2 atoms
TEST(forces2, force) {
    mdsys_t sys;
    sys.natoms = 2;
    sys.rx = new double[2];
	sys.ry = new double[2];
	sys.rz = new double[2];
    sys.fx = new double[2];
    sys.fy = new double[2];
    sys.fz = new double[2];
    sys.rx[0] = 0.0;
    sys.rx[1] = 0.0;
    sys.ry[0] = 0.0;
    sys.ry[1] = 0.0;
    sys.rz[0] = 0.0;
    sys.rz[1] = 3.0;
    sys.fx[0] = 0.0;
    sys.fx[1] = 0.0;
    sys.fy[0] = 0.0;
    sys.fy[1] = 0.0;
    sys.fz[0] = 0.0;
    sys.fz[1] = 0.0;
    sys.mass = 39.948;
    sys.epsilon = 0.2379;
    sys.sigma = 3.405;
    sys.rcut = 8.5;
    sys.box = 17.1580;
    force(&sys);
    EXPECT_EQ(0.0, sys.fx[0]);
    EXPECT_EQ(0.0, sys.fx[1]);
    EXPECT_EQ(0.0, sys.fy[0]);
    EXPECT_EQ(0.0, sys.fy[1]);
    EXPECT_NEAR(-13.32787861802134, sys.fz[0], 10e-10); // Optimization creates a small error in last digit.
    EXPECT_NEAR(13.32787861802134, sys.fz[1], 10e-10);
    //EXPECT_DOUBLE_EQ(-13.32787861802134, sys.fz[0]);
    //EXPECT_DOUBLE_EQ(13.32787861802134, sys.fz[1]);
}
//test forces.c for the function force, which is used to compute the forces on each atom of a system of 3 atoms
