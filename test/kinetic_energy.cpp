#include <gtest/gtest.h>
extern "C" { // since our lib is C code, but the test here is C++
#include <mdlib-util.h>
#include <mdlib.h>
}

TEST(KineticEnergy, ComputeKE) {
	mdsys_t sys;
	sys.natoms = 2;
	sys.mass = 39.948; //argon
	sys.vx = new double[2];
	sys.vy = new double[2];
	sys.vz = new double[2];
	sys.vx[0] = 1.5;
	sys.vy[0] = 2.5;
	sys.vz[0] = 3.5;
	sys.vx[1] = 10.5;
	sys.vy[1] = 15.6;
	sys.vz[1] = 20.7;
	ekin(&sys);
	ASSERT_DOUBLE_EQ(sys.ekin, 38327260.75777286);
}

