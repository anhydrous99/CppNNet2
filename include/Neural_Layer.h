//
// Created by ubrdog on 12/25/18.
//

#ifndef CPPNNET2_NEURAL_LAYER_H
#define CPPNNET2_NEURAL_LAYER_H

#include "Layer.h"

namespace CppNNet2 {

  class Neural_Layer : public Layer {
    Matrix <netfl> Weights;
    Matrix <netfl> Biases;
  public:
    Neural_Layer(Matrix <netfl> &weights, Matrix <netfl> &biases);
    Neural_Layer(Matrix <netfl> &weights, Matrix <netfl> &biases, layer_placement placement,
                 std::shared_ptr<Layer> layer);
    Neural_Layer(Matrix <netfl> &weights, Matrix <netfl> &biases, std::shared_ptr<Layer> previous_layer,
                 std::shared_ptr<Layer> next_layer);

    Matrix <netfl> feedforward(Matrix <netfl> &input);
    Matrix <netfl> feedbackward(Matrix <netfl> &input);
  };

  Neural_Layer::Neural_Layer(CppNNet2::Matrix<netfl> &weights, CppNNet2::Matrix<netfl> &biases) {
    Weights = weights;
    Biases = biases;
  }

  Neural_Layer::Neural_Layer(CppNNet2::Matrix<netfl> &weights, CppNNet2::Matrix<netfl> &biases,
                             CppNNet2::layer_placement placement, std::shared_ptr<CppNNet2::Layer> layer) :
      Layer(placement, layer) {
    Weights = weights;
    Biases = biases;
  }

  Neural_Layer::Neural_Layer(Matrix <netfl> &weights, Matrix <netfl> &biases,
                             std::shared_ptr<Layer> previous_layer,
                             std::shared_ptr<Layer> next_layer) : Layer(previous_layer, next_layer) {
    Weights = weights;
    Biases = biases;
  }

  Matrix <netfl> Neural_Layer::feedforward(Matrix <netfl> &input) {
    return Weights * input + Biases;
  }

  Matrix <netfl> Neural_Layer::feedbackward(Matrix <netfl> &input) {
    // TODO
  }
}

#endif //CPPNNET2_NEURAL_LAYER_H
