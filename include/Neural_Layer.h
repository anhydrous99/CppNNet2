//
// Created by ubrdog on 12/25/18.
//

#ifndef CPPNNET2_NEURAL_LAYER_H
#define CPPNNET2_NEURAL_LAYER_H

#include "Layer.h"

namespace CppNNet2 {

  class Neural_Layer : public Layer {
    Matrix <netfl> Weights;
    Matrix <netfl> n;
    bool n_evaluated = false;
  public:

    Matrix <netfl> feedforward(Matrix <netfl> &input) override;
  };


  Matrix <netfl> Neural_Layer::feedforward(Matrix <netfl> &input) {
    if (!n_evaluated) {
      n = Weights * ((Previous_layer) ? feedforward(input) : input) + Biases;
      n_evaluated = true;
    }
    return n;
  }
}

#endif //CPPNNET2_NEURAL_LAYER_H
