//
// Created by ubrdog on 12/25/18.
//

#include "Matrix.h"
#include <gtest/gtest.h>

using namespace CppNNet2;

TEST(MatrixMultiplication, Float) {
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

TEST(MatrixMultiplication, Double) {
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

TEST(MatrixMultiplication, Integer) {
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

TEST(MatrixConcatenation, Horizontal) {
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

TEST(MatrixConcatenation, Vertical) {
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

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}