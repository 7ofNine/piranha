/***************************************************************************
 *   Copyright (C) 2007, 2008 by Francesco Biscani   *
 *   bluescarni@gmail.com   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <boost/numeric/conversion/cast.hpp>
#include <climits>
#include <cstddef>
#include <thread>

#include "core/config.h"
#include "core/exceptions.h"
#include "core/memory.h"
#include "core/settings.h"

/////////////////////////////////////////////////////////////////////////////////////
//
// This file is generated from settings.cpp.template
//
// Any changes directly to this file will get lost
/////////////////////////////////////////////////////////////////////////////////////

namespace piranha
{
        // Settings' static members.
        std::size_t settings::m_memory_limit = 15500000000u; // ~ 15.5GByte
        double settings::m_numerical_zero = 1E-80;
        bool settings::m_debug = false;
        const std::string settings::m_version = "@PIRANHA_VERSION_STRING@";
        const std::size_t settings::cache_size = PIRANHA_CACHE_SIZE;
        bool settings::blocker = false;
        unsigned settings::m_nthread = std::thread::hardware_concurrency() ? std::thread::hardware_concurrency() : 1;
        settings::startup_class settings::startup;
        std::size_t settings::m_max_pretty_print_size = 500;
        settings::MultiplicationAlgorithm settings::multiplicationAlgorithm = settings::MultiplicationAlgorithm::AUTOMATIC;

        settings::startup_class::startup_class()
        {
                // Some static checks.
                // This check is here because in integer_typedefs we select integers based on the number of bits.
                static_assert(CHAR_BIT == 8, "Invalid size of char.");
                // This is because somewhere we use chars as bools.
                static_assert(sizeof(char) == sizeof(bool), "Wrong char-bool size ratio.");
                // TODO: where is this used? Find out and wipe it.
                static_assert(sizeof(std::size_t) == sizeof(void *), "std::size_t and void * are not the same size.");
                static_assert(PIRANHA_MAX_ECHELON_LEVEL >= 0, "Max echelon level must be nonnegative.");
                static_assert(settings::cache_size > 0 && lg<settings::cache_size>::value > 1, "Invalid value for cache size.");
                // Startup report.
                std::cout << "Piranha version         : " << "@PIRANHA_VERSION_STRING@" << std::endl;
                std::cout << "Piranha GIT revision    : " << "@PIRANHA_GIT_REVISION@" << std::endl;
                std::cout << "Number of cores detected: " <<
                        (std::thread::hardware_concurrency() ? std::thread::hardware_concurrency() : 1) << std::endl;
                std::cout << std::endl << std::endl << 
                                                         "--------------------------" << std::endl <<
                                                         "! Piranha core is ready. ! " << std::endl;
                                            std::cout << "__________________________" << std::endl << std::flush;
        }

        /// Get version.
        const char *settings::get_version()
        {
                return m_version.c_str();
        }

        std::size_t settings::get_max_pretty_print_size()
        {
                return m_max_pretty_print_size;
        }

        void settings::set_max_pretty_print_size(std::size_t n)
        {
                if (n < 10) {
                        PIRANHA_THROW(value_error,"Invalid max size for pretty printing, "
                                "please use an integer greater than 10");
                }
                m_max_pretty_print_size = boost::numeric_cast<std::size_t>(n);
        }

        void settings::set_nthread(const unsigned int &n)
        {
                m_nthread = n;
        }

        std::size_t settings::get_nthread()
        {
                return m_nthread;
        }
}
