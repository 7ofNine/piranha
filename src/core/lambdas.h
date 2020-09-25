#ifndef PIRANHA_LAMBDAS_H
#define PIRANHA_LAMBDAS_H

#include <cstddef>


namespace piranha {

    auto increment = [](const std::size_t x) -> std::size_t { return x + 1; }; // increment functor

}

#endif PIRANHA_LAMBDAS_H