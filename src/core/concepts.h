#pragma once

//
// first attempt to play with simple concepts in piramha
//

#include<type_traits>
#include<concepts>

#include "base_classes/base_series.h"


namespace piranha {

// test on valid EchelonLevel before we ultimately attempt to remove it
//
// the max echelonLevel is set to 1
// we intend to actually remove it completely
//     doesn't work for non type parameters
//template <unsigned int T>
//concept EchelonLevel = T <= 1; 
    //not really a concept. non type parameters don't work in concepts but in requires statements
    constexpr unsigned int maxEchelonLevel = 1;

    //template <typename T, typename S>
    //concept SameLevel = T::echelonLevel == S::echelonLevel;
    
}
