//
// Created by ubrdog on 12/25/18.
//

#include "CppNNet2.h"
#include <gtest/gtest.h>

using namespace CppNNet2;

TEST(Matrix, FloatMultiplication) {
  Matrix<float> smat(2, 2);
  smat(0, 0) = 1.0;
  smat(0, 1) = 2.0;
  smat(1, 0) = 3.0;
  smat(1, 1) = 4.0;

  Matrix<float> result(2, 2);
  result(0, 0) = 7.0;
  result(0, 1) = 10.0;
  result(1, 0) = 15.0;
  result(1, 1) = 22.0;

  Matrix<float> to_check = smat * smat;
  Matrix<bool> exptrue = to_check == result;

  for (size_t i = 0; i < exptrue.size(); i++) {
    EXPECT_TRUE(exptrue[i]);
  }
}

TEST(Matrix, DoubleMultiplication) {
  Matrix<double> smat(2, 2);
  smat(0, 0) = 1.0;
  smat(0, 1) = 2.0;
  smat(1, 0) = 3.0;
  smat(1, 1) = 4.0;

  Matrix<double> result(2, 2);
  result(0, 0) = 7.0;
  result(0, 1) = 10.0;
  result(1, 0) = 15.0;
  result(1, 1) = 22.0;

  Matrix<double> to_check = smat * smat;
  Matrix<bool> exptrue = to_check == result;

  for (size_t i = 0; i < exptrue.size(); i++) {
    EXPECT_TRUE(exptrue[i]);
  }
}

TEST(Matrix, IntegerMultiplication) {
  Matrix<int> smat(2, 2);
  smat(0, 0) = 1;
  smat(0, 1) = 2;
  smat(1, 0) = 3;
  smat(1, 1) = 4;

  Matrix<int> result(2, 2);
  result(0, 0) = 7;
  result(0, 1) = 10;
  result(1, 0) = 15;
  result(1, 1) = 22;

  Matrix<int> to_check = smat * smat;
  Matrix<bool> exptrue = to_check == result;

  for (size_t i = 0; i < exptrue.size(); i++) {
    EXPECT_TRUE(exptrue[i]);
  }
}

TEST(Matrix, HorizontalConcatenation) {
  Matrix<float> smat1(2, 2);
  Matrix<float> smat2(2, 2);
  Matrix<float> result(2, 4);
  smat1(0, 0) = 5.0;
  smat1(0, 1) = 2.5;
  smat1(1, 0) = 8.6;
  smat1(1, 1) = 8.9;
  smat2(0, 0) = 5.5;
  smat2(0, 1) = 2.8;
  smat2(1, 0) = 8.2;
  smat2(1, 1) = 8.6;
  //
  result(0, 0) = 5.0;
  result(0, 1) = 2.5;
  result(1, 0) = 8.6;
  result(1, 1) = 8.9;
  result(0, 2) = 5.5;
  result(0, 3) = 2.8;
  result(1, 2) = 8.2;
  result(1, 3) = 8.6;

  Matrix<float> to_check = concat(HORIZONTAL, smat1, smat2);
  EXPECT_EQ(result.cols(), to_check.cols());
  EXPECT_EQ(result.rows(), to_check.rows());

  Matrix<bool> boolMatrix = result == to_check;
  for (size_t i = 0; i < boolMatrix.size(); i++) {
    EXPECT_TRUE(boolMatrix[i]);
  }
}

TEST(Matrix, VerticalConcatenation) {
  Matrix<float> smat1(2, 2);
  Matrix<float> smat2(2, 2);
  Matrix<float> result(4, 2);
  smat1(0, 0) = 5.0;
  smat1(0, 1) = 2.5;
  smat1(1, 0) = 8.6;
  smat1(1, 1) = 8.9;
  smat2(0, 0) = 5.5;
  smat2(0, 1) = 2.8;
  smat2(1, 0) = 8.2;
  smat2(1, 1) = 8.6;
  //
  result(0, 0) = 5.0;
  result(0, 1) = 2.5;
  result(1, 0) = 8.6;
  result(1, 1) = 8.9;
  result(2, 0) = 5.5;
  result(2, 1) = 2.8;
  result(3, 0) = 8.2;
  result(3, 1) = 8.6;

  Matrix<float> to_check = concat(VERTICAL, smat1, smat2);
  EXPECT_EQ(result.cols(), to_check.cols());
  EXPECT_EQ(result.rows(), to_check.rows());

  Matrix<bool> boolMatrix = result == to_check;
  for (size_t i = 0; i < boolMatrix.size(); i++) {
    EXPECT_TRUE(boolMatrix[i]);
  }
}

TEST(Matrix, SimpleMathOperations) {
  Matrix<double> mat1(1, 1);
  Matrix<double> mat2(1, 1);
  double abs_error = 0.001;
  double pw = 1.5;

  mat1[0] = -1.0;
  mat2[0] = 5.0;

  EXPECT_EQ(abs(mat1)[0], 1);
  EXPECT_NEAR(acos(mat1)[0], 3.1415, abs_error);
  EXPECT_NEAR(asin(mat1)[0], -1.5707, abs_error);
  EXPECT_NEAR(cos(mat1)[0], 0.5403, abs_error);
  EXPECT_NEAR(sin(mat1)[0], -0.8415, abs_error);
  EXPECT_NEAR(cosh(mat1)[0], 1.5431, abs_error);
  EXPECT_NEAR(sinh(mat1)[0], -1.1752, abs_error);
  EXPECT_NEAR(tan(mat1)[0], -1.5574, abs_error);
  EXPECT_NEAR(tanh(mat1)[0], -0.7616, abs_error);
  EXPECT_NEAR(exp(mat1)[0], 0.3679, abs_error);
  EXPECT_NEAR(log(mat2)[0], 1.6094, abs_error);
  EXPECT_NEAR(log10(mat2)[0], 0.6989, abs_error);
  EXPECT_NEAR(sqrt(mat2)[0], 2.2361, abs_error);
  EXPECT_NEAR(pow(mat2, mat2)[0], 3125.0, abs_error);
  EXPECT_NEAR(pow(mat2, pw)[0], 11.1803, abs_error);
  EXPECT_NEAR(pow(pw, mat2)[0], 7.5937, abs_error);
}

TEST(Matrix, Transpose) {
  Matrix<int> mat1(3, 2);
  Matrix<int> mat2(2, 3);

  mat1(0, 0) = 0;
  mat1(0, 1) = 1;
  mat1(1, 0) = 2;
  mat1(1, 1) = 3;
  mat1(2, 0) = 4;
  mat1(2, 1) = 5;

  mat2(0, 0) = 0;
  mat2(0, 1) = 2;
  mat2(0, 2) = 4;
  mat2(1, 0) = 1;
  mat2(1, 1) = 3;
  mat2(1, 2) = 5;

  Matrix<int> mat3 = transpose(mat1);
  Matrix<bool> boolMatrix = mat3 == mat2;
  for (size_t i = 0; i < boolMatrix.size(); i++) {
    EXPECT_TRUE(boolMatrix[i]);
  }
}

TEST(Matrix, FloatSum) {
  Matrix<float> smat(2, 2);
  smat(0, 0) = 1.0;
  smat(0, 1) = 2.0;
  smat(1, 0) = 3.0;
  smat(1, 1) = 4.0;

  Matrix<float> result(2, 2);
  result(0, 0) = 2.0;
  result(0, 1) = 4.0;
  result(1, 0) = 6.0;
  result(1, 1) = 8.0;

  Matrix<float> to_check = smat + smat;
  Matrix<bool> exptrue = to_check == result;

  for (size_t i = 0; i < exptrue.size(); i++) {
    EXPECT_TRUE(exptrue[i]);
  }
}

TEST(Matrix, DoubleSum) {
  Matrix<double> smat(2, 2);
  smat(0, 0) = 1.0;
  smat(0, 1) = 2.0;
  smat(1, 0) = 3.0;
  smat(1, 1) = 4.0;

  Matrix<double> result(2, 2);
  result(0, 0) = 2.0;
  result(0, 1) = 4.0;
  result(1, 0) = 6.0;
  result(1, 1) = 8.0;

  Matrix<double> to_check = smat + smat;
  Matrix<bool> exptrue = to_check == result;

  for (size_t i = 0; i < exptrue.size(); i++) {
    EXPECT_TRUE(exptrue[i]);
  }
}

TEST(Matrix, IntegerSum) {
  Matrix<int> smat(2, 2);
  smat(0, 0) = 1;
  smat(0, 1) = 2;
  smat(1, 0) = 3;
  smat(1, 1) = 4;

  Matrix<int> result(2, 2);
  result(0, 0) = 2;
  result(0, 1) = 4;
  result(1, 0) = 6;
  result(1, 1) = 8;

  Matrix<int> to_check = smat + smat;
  Matrix<bool> exptrue = to_check == result;

  for (size_t i = 0; i < exptrue.size(); i++) {
    EXPECT_TRUE(exptrue[i]);
  }
}

TEST(Matrix, FloatSubtraction) {
  Matrix<float> smat(2, 2);
  smat(0, 0) = 1.0;
  smat(0, 1) = 2.0;
  smat(1, 0) = 3.0;
  smat(1, 1) = 4.0;

  Matrix<float> result(2, 2);
  result(0, 0) = 0.0;
  result(0, 1) = 0.0;
  result(1, 0) = 0.0;
  result(1, 1) = 0.0;

  Matrix<float> to_check = smat - smat;
  Matrix<bool> exptrue = to_check == result;

  for (size_t i = 0; i < exptrue.size(); i++) {
    EXPECT_TRUE(exptrue[i]);
  }
}

TEST(Matrix, DoubleSubtraction) {
  Matrix<double> smat(2, 2);
  smat(0, 0) = 1.0;
  smat(0, 1) = 2.0;
  smat(1, 0) = 3.0;
  smat(1, 1) = 4.0;

  Matrix<double> result(2, 2);
  result(0, 0) = 0.0;
  result(0, 1) = 0.0;
  result(1, 0) = 0.0;
  result(1, 1) = 0.0;

  Matrix<double> to_check = smat - smat;
  Matrix<bool> exptrue = to_check == result;

  for (size_t i = 0; i < exptrue.size(); i++) {
    EXPECT_TRUE(exptrue[i]);
  }
}

TEST(Matrix, IntegerSubtraction) {
  Matrix<int> smat(2, 2);
  smat(0, 0) = 1;
  smat(0, 1) = 2;
  smat(1, 0) = 3;
  smat(1, 1) = 4;

  Matrix<int> result(2, 2);
  result(0, 0) = 0;
  result(0, 1) = 0;
  result(1, 0) = 0;
  result(1, 1) = 0;

  Matrix<int> to_check = smat - smat;
  Matrix<bool> exptrue = to_check == result;

  for (size_t i = 0; i < exptrue.size(); i++) {
    EXPECT_TRUE(exptrue[i]);
  }
}

TEST(Matrix, KroneckerProduct) {
  Matrix<float> smat(2, 2);
  smat(0, 0) = 1.0;
  smat(0, 1) = 2.0;
  smat(1, 0) = 3.0;
  smat(1, 1) = 4.0;

  Matrix<float> result(4, 4);
  result(0, 0) = 1.0;
  result(0, 1) = 2.0;
  result(0, 2) = 2.0;
  result(0, 3) = 4.0;
  result(1, 0) = 3.0;
  result(1, 1) = 4.0;
  result(1, 2) = 6.0;
  result(1, 3) = 8.0;
  result(2, 0) = 3.0;
  result(2, 1) = 6.0;
  result(2, 2) = 4.0;
  result(2, 3) = 8.0;
  result(3, 0) = 9.0;
  result(3, 1) = 12.0;
  result(3, 2) = 12.0;
  result(3, 3) = 16.0;

  Matrix<float> to_check = kronecker_product(smat, smat);
  Matrix<bool> exptrue = to_check == result;

  for (size_t i = 0; i < exptrue.size(); i++) {
    EXPECT_TRUE(exptrue[i]);
  }
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}