//
// Created by ubrdog on 12/25/18.
//

#ifndef CPPNNET2_LAYER_H
#define CPPNNET2_LAYER_H

#include <memory>

class Layer : private std::enable_shared_from_this<Layer> {
    std::shared_ptr<Layer> previous_layer = nullptr;
    std::shared_ptr<Layer> next_layer = nullptr;
};

#endif //CPPNNET2_LAYER_H
