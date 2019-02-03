//
// Created by ubrdog on 12/25/18.
//

#ifndef CPPNNET2_LAYER_H
#define CPPNNET2_LAYER_H

#include "typedefs.h"
#include "Matrix.h"
#include <memory>

namespace CppNNet2 {
  enum layer_placement {
    NEXT, PREVIOUS
  };

  class Layer : private std::enable_shared_from_this<Layer> {
  protected:
    // data:
    std::shared_ptr<Layer> Previous_layer = nullptr;

  public:
    // constructor & destructor:
    Layer() = default;

    explicit Layer(std::shared_ptr<Layer> &previous_layer);

    ~Layer() = default;

    // feeding functions:
    virtual Matrix <netfl> feedforward(Matrix <netfl> &input);

    // msc:
    std::shared_ptr<Layer> Get_Previous_Layer();

    std::shared_ptr<Layer> getptr();
  };

  Layer::Layer(std::shared_ptr<CppNNet2::Layer> &previous_layer) {
    Previous_layer = previous_layer;
  }

  Matrix <netfl> Layer::feedforward(Matrix <netfl> &input) {
    return (Previous_layer) ? Previous_layer->feedforward(input) : input;
  }

  std::shared_ptr<Layer> Layer::Get_Previous_Layer() {
    return Previous_layer;
  }

  std::shared_ptr<Layer> Layer::getptr() {
    return shared_from_this();
  }
}

#endif //CPPNNET2_LAYER_H
