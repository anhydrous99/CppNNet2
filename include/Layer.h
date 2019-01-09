//
// Created by ubrdog on 12/25/18.
//

#ifndef CPPNNET2_LAYER_H
#define CPPNNET2_LAYER_H

#include "typedefs.h"
#include "Matrix.h"
#include <memory>

namespace CppNNet2 {
  enum layer_placement { NEXT, PREVIOUS };
    class Layer : private std::enable_shared_from_this<Layer> {
        // data:
        std::shared_ptr<Layer> Previous_layer = nullptr;
        std::shared_ptr<Layer> Next_layer = nullptr;

    public:
        // constructor & destructor:
        Layer() = default;
        Layer(std::shared_ptr<Layer> &previous_layer, std::shared_ptr<Layer> &next_layer);
        Layer(layer_placement placement, std::shared_ptr<Layer> &layer);
        ~Layer() = default;

        // feeding functions:
        virtual Matrix<netfl> feedforward(Matrix<netfl> &input);
        virtual Matrix<netfl> feedbackward(Matrix<netfl> &input);

        // msc:
        std::shared_ptr<Layer> getptr();
    };

    Layer::Layer(std::shared_ptr<CppNNet2::Layer> &previous_layer, std::shared_ptr<CppNNet2::Layer> &next_layer) {
        Previous_layer = previous_layer;
        Next_layer = next_layer;
    }

    Layer::Layer(CppNNet2::layer_placement placement, std::shared_ptr<CppNNet2::Layer> &layer) {
      if (placement == NEXT)
        Next_layer = layer;
      else
        Previous_layer = layer;
    }

    Matrix<netfl> Layer::feedforward(Matrix<netfl> &input) {
        return (Previous_layer) ? Previous_layer->feedforward(input) : input;
    }

    Matrix<netfl> Layer::feedbackward(Matrix<netfl> &input) {
        return (Next_layer) ? Next_layer->feedbackward(input) : input;
    }

    std::shared_ptr<Layer> Layer::getptr() {
        return shared_from_this();
    }
}

#endif //CPPNNET2_LAYER_H
