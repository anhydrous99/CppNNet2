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

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}