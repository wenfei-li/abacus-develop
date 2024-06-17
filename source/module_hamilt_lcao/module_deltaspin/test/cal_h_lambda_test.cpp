#include "../spin_constrain.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "module_cell/klist.h"
K_Vectors::K_Vectors(){}
K_Vectors::~K_Vectors(){}

/************************************************
 *  unit test of the function of cal_h_lambda
 ***********************************************/

/**
 * Tested function:
 *  - SpinConstrain::cal_h_lambda
 *    - this function calculates the h_lambda operator
 */

class SpinConstrainTest : public testing::Test
{
  protected:
    SpinConstrain<std::complex<double>, base_device::DEVICE_CPU>& sc
        = SpinConstrain<std::complex<double>, base_device::DEVICE_CPU>::getScInstance();
};

TEST_F(SpinConstrainTest, CalHLambda)
{
    Parallel_Orbitals paraV;
    std::ofstream ofs("test.log");
    int nrow = 2;
    int ncol = 2;
    sc.set_nspin(4);
    paraV.set_serial(nrow, ncol);
    sc.set_ParaV(&paraV);
    EXPECT_EQ(sc.ParaV->nloc, 4);
    std::map<int, int> atomCounts = {
        {0, 1}
    };
    std::map<int, int> orbitalCounts = {
        {0, 1}
    };
    sc.set_atomCounts(atomCounts);
    sc.set_orbitalCounts(orbitalCounts);
    sc.set_npol(2);
    std::vector<std::complex<double>> h_lambda(sc.ParaV->nloc);
    std::fill(h_lambda.begin(), h_lambda.end(), std::complex<double>(0, 0));
    std::vector<std::complex<double>> Sloc2 = {
        std::complex<double>{1.0, 0.0},
        std::complex<double>{0.0, 0.0},
        std::complex<double>{0.0, 0.0},
        std::complex<double>{1.0, 0.0}
    };
    ModuleBase::Vector3<double>* sc_lambda = new ModuleBase::Vector3<double>[1];
    sc_lambda[0][0] = 1.0;
    sc_lambda[0][1] = 1.0;
    sc_lambda[0][2] = 1.0;
    sc.set_sc_lambda(sc_lambda, 1);
    // column_major = true
    sc.cal_h_lambda(&h_lambda[0], Sloc2, true, 0);
    // h_lambda = - [lambda_x * sigma_x + lambda_y * sigma_y + lambda_z * sigma_z] * Sloc2
    std::vector<std::complex<double>> columnMajor_h_lambda = {
        std::complex<double>{-1.0, 0.0 },
        std::complex<double>{-1.0, 1.0 },
        std::complex<double>{-1.0, -1.0},
        std::complex<double>{1.0,  0.0 }
    };
    EXPECT_DOUBLE_EQ(h_lambda[0].real(), columnMajor_h_lambda[0].real());
    EXPECT_DOUBLE_EQ(h_lambda[0].imag(), columnMajor_h_lambda[0].imag());
    EXPECT_DOUBLE_EQ(h_lambda[1].real(), columnMajor_h_lambda[1].real());
    EXPECT_DOUBLE_EQ(h_lambda[1].imag(), columnMajor_h_lambda[1].imag());
    EXPECT_DOUBLE_EQ(h_lambda[2].real(), columnMajor_h_lambda[2].real());
    EXPECT_DOUBLE_EQ(h_lambda[2].imag(), columnMajor_h_lambda[2].imag());
    EXPECT_DOUBLE_EQ(h_lambda[3].real(), columnMajor_h_lambda[3].real());
    EXPECT_DOUBLE_EQ(h_lambda[3].imag(), columnMajor_h_lambda[3].imag());
    // column_major = false
    delete[] sc_lambda;
    sc.cal_h_lambda(&h_lambda[0], Sloc2, false, 0);
    std::vector<std::complex<double>> rowMajor_h_lambda = {
        std::complex<double>{-1.0, 0.0 },
        std::complex<double>{-1.0, -1.0},
        std::complex<double>{-1.0, 1.0 },
        std::complex<double>{1.0,  0.0 }
    };
    EXPECT_DOUBLE_EQ(h_lambda[0].real(), rowMajor_h_lambda[0].real());
    EXPECT_DOUBLE_EQ(h_lambda[0].imag(), rowMajor_h_lambda[0].imag());
    EXPECT_DOUBLE_EQ(h_lambda[1].real(), rowMajor_h_lambda[1].real());
    EXPECT_DOUBLE_EQ(h_lambda[1].imag(), rowMajor_h_lambda[1].imag());
    EXPECT_DOUBLE_EQ(h_lambda[2].real(), rowMajor_h_lambda[2].real());
    EXPECT_DOUBLE_EQ(h_lambda[2].imag(), rowMajor_h_lambda[2].imag());
    EXPECT_DOUBLE_EQ(h_lambda[3].real(), rowMajor_h_lambda[3].real());
    EXPECT_DOUBLE_EQ(h_lambda[3].imag(), rowMajor_h_lambda[3].imag());
    remove("test.log");
}

TEST_F(SpinConstrainTest, CalHLambdaS2)
{
    Parallel_Orbitals paraV;
    std::ofstream ofs("test.log");
    int nrow = 1;
    int ncol = 1;
    sc.set_nspin(2);
    paraV.set_serial(nrow, ncol);
    sc.set_ParaV(&paraV);
    EXPECT_EQ(sc.ParaV->nloc, 1);
    std::map<int, int> atomCounts = {
        {0, 1}
    };
    std::map<int, int> orbitalCounts = {
        {0, 1}
    };
    sc.set_atomCounts(atomCounts);
    sc.set_orbitalCounts(orbitalCounts);
    sc.set_npol(1);
    std::vector<std::complex<double>> h_lambda(sc.ParaV->nloc);
    std::fill(h_lambda.begin(), h_lambda.end(), std::complex<double>(0, 0));
    std::vector<std::complex<double>> Sloc2 = {
        std::complex<double>{1.0, 0.0}
    };
    ModuleBase::Vector3<double>* sc_lambda = new ModuleBase::Vector3<double>[1];
    sc_lambda[0][0] = 0.0;
    sc_lambda[0][1] = 0.0;
    sc_lambda[0][2] = 2.0;
    sc.set_sc_lambda(sc_lambda, 1);
    // column_major = true
    sc.cal_h_lambda(&h_lambda[0], Sloc2, true, 0);
    // h_lambda = -Sloc2[icc] * this->lambda_[iat2][2]
    std::vector<std::complex<double>> columnMajor_h_lambda = {
        std::complex<double>{-2.0, 0.0 },
    };
    EXPECT_DOUBLE_EQ(h_lambda[0].real(), columnMajor_h_lambda[0].real());
    EXPECT_DOUBLE_EQ(h_lambda[0].imag(), columnMajor_h_lambda[0].imag());
    sc.cal_h_lambda(&h_lambda[0], Sloc2, true, 1);
    // h_lambda = -Sloc2[icc] * (-this->lambda_[iat2][2])
    std::vector<std::complex<double>> columnMajor_h_lambda1 = {
        std::complex<double>{2.0, 0.0 },
    };
    EXPECT_DOUBLE_EQ(h_lambda[0].real(), columnMajor_h_lambda1[0].real());
    EXPECT_DOUBLE_EQ(h_lambda[0].imag(), columnMajor_h_lambda1[0].imag());
    // column_major = false
    delete[] sc_lambda;
    // h_lambda = -Sloc2[icc] * this->lambda_[iat2][2]
    sc.cal_h_lambda(&h_lambda[0], Sloc2, false, 0);
    std::vector<std::complex<double>> rowMajor_h_lambda = {
        std::complex<double>{-2.0, 0.0 },
    };
    EXPECT_DOUBLE_EQ(h_lambda[0].real(), rowMajor_h_lambda[0].real());
    EXPECT_DOUBLE_EQ(h_lambda[0].imag(), rowMajor_h_lambda[0].imag());
    // h_lambda = -Sloc2[icc] * (-this->lambda_[iat2][2])
    sc.cal_h_lambda(&h_lambda[0], Sloc2, false, 1);
    std::vector<std::complex<double>> rowMajor_h_lambda1 = {
        std::complex<double>{2.0, 0.0 },
    };
    EXPECT_DOUBLE_EQ(h_lambda[0].real(), rowMajor_h_lambda1[0].real());
    EXPECT_DOUBLE_EQ(h_lambda[0].imag(), rowMajor_h_lambda1[0].imag());
    remove("test.log");
}
