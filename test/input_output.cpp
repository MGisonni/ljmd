#include <gtest/gtest.h>
extern "C" { // since our lib is C code, but the test here is C++
#include <mdlib-util.h>
#include <mdlib.h>
}

TEST(Input, GetLine) {
	char input_filename[] = "test_input_17.inp";
	char line[BLEN];
	FILE *fp;
	fp = fopen(input_filename, "w");
	ASSERT_NE(fp, nullptr);
	fprintf(fp,"%s","17               # natoms");
	fclose(fp);
	fp = fopen(input_filename, "r");
	get_a_line(fp,line);
	ASSERT_STREQ(line, "17");

	remove(input_filename);
}

TEST(Output, WriteToFile) {
	char erg_filename[] = "test_erg.out";
	char traj_filename[] = "test_traj.out";
	mdsys_t sys;
    sys.natoms = 2;
    sys.rx = new double[2];
    sys.ry = new double[2];
    sys.rz = new double[2];
    sys.rx[0] = 10.1;
	sys.ry[0] = 11.2;
    sys.rz[0] = 12.3;
    sys.rx[1] = 20.2;
    sys.ry[1] = 21.3;
    sys.rz[1] = 22.4;
    sys.nfi = 3;
	sys.temp = 1.5;
	sys.ekin = 2.6;
	sys.epot = 3.7;
	FILE *fp_erg = fopen(erg_filename, "w");
	ASSERT_NE(fp_erg, nullptr);
	FILE *fp_traj = fopen(traj_filename, "w");
    ASSERT_NE(fp_traj, nullptr);
	output(&sys, fp_erg, fp_traj);
	fclose(fp_erg);
	fclose(fp_traj);

	fp_erg = fopen(erg_filename, "r");
	ASSERT_NE(fp_erg, nullptr);
	double temp_sum;
    fscanf(fp_erg,"%d %lf %lf %lf %lf\n", &sys.nfi, &sys.temp, &sys.ekin, &sys.epot, &temp_sum);
	ASSERT_EQ(sys.nfi, 3);
	ASSERT_DOUBLE_EQ(sys.temp, 1.5);
	ASSERT_DOUBLE_EQ(sys.ekin, 2.6);
	ASSERT_DOUBLE_EQ(sys.epot, 3.7);
	ASSERT_DOUBLE_EQ(temp_sum, 6.3);
	fclose(fp_erg);

	fp_traj = fopen(traj_filename, "r");
	ASSERT_NE(fp_traj, nullptr);
	fscanf(fp_traj,"%d nfi=%d etot=%lf\n", &sys.natoms, &sys.nfi, &temp_sum);
	ASSERT_EQ(sys.natoms, 2);
	ASSERT_EQ(sys.nfi, 3);
	ASSERT_DOUBLE_EQ(temp_sum, 6.3);
	for (int i=0; i<sys.natoms; ++i) {
		fscanf(fp_traj,"Ar  %lf %lf %lf", &sys.rx[i], &sys.ry[i], &sys.rz[i]);
	}
	ASSERT_DOUBLE_EQ(sys.rx[0], 10.1);
	ASSERT_DOUBLE_EQ(sys.ry[0], 11.2);
	ASSERT_DOUBLE_EQ(sys.rz[0], 12.3);
	ASSERT_DOUBLE_EQ(sys.rx[1], 20.2);
	ASSERT_DOUBLE_EQ(sys.ry[1], 21.3);
	ASSERT_DOUBLE_EQ(sys.rz[1], 22.4);
	fclose(fp_traj);

	remove(erg_filename);
	remove(traj_filename);
}
