#include <gtest/gtest.h>
extern "C" { // since our lib is C code, but the test here is C++
#include <mdlib-util.h>
#include <mdlib.h>
}

TEST(Utilities, AzzeroArray) {
	int n = 3;
	double *arr = new double[n];
	arr[0] = 1.0;
	arr[1] = 2.0;
	arr[2] = 3.0;
	azzero(arr, n);
	ASSERT_DOUBLE_EQ(arr[0], 0.0);
	ASSERT_DOUBLE_EQ(arr[1], 0.0);
	ASSERT_DOUBLE_EQ(arr[2], 0.0);
}

TEST(Utilities, WallClock) {
	double t1 = wallclock();
	sleep(1);
	double t2 = wallclock();
	ASSERT_TRUE(t2 > t1+0.9 && t2 < t2+1.1);
}

