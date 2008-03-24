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

#include <boost/integer_traits.hpp>

#include "core/memory.h"
#include "core/piranha_gsl.h"
#include "core/settings_manager.h"
#include "core/stream_manager.h"
#include "core/version.h"

namespace piranha
{
  // Settings manager's static members.
  double settings_manager::hash_max_load_factor = 0.5;
  double settings_manager::numerical_zero = 1E-80;
  const max_fast_uint settings_manager::min_u = boost::integer_traits<max_fast_uint>::min();
  const max_fast_uint settings_manager::max_u = boost::integer_traits<max_fast_uint>::max();
  const max_fast_int settings_manager::min_i = boost::integer_traits<max_fast_int>::min();
  const max_fast_int settings_manager::max_i = boost::integer_traits<max_fast_int>::max();
  const unsigned int settings_manager::jacang_limit = 20;
  std::string settings_manager::path = _PIRANHA_DEFAULT_PATH;
  const std::string settings_manager::default_path = _PIRANHA_DEFAULT_PATH;
  const std::string settings_manager::version = __PIRANHA_VERSION;
  bool settings_manager::enable_progress_display = true;
  settings_manager::startup_class settings_manager::startup;
#ifdef _PIRANHA_TBB
  const tbb::task_scheduler_init settings_manager::tbb_init;
#endif

  settings_manager::startup_class::startup_class()
  {
    // Startup report.
    std::cout << "This is Piranha version " << get_version() << std::endl;
    std::cout << "Default parameters initialized:" << std::endl;
    std::cout << "Print precision\t\t\t=\t" << stream_manager::digits() << std::endl;
    std::cout << "Numerical zero\t\t\t=\t" << get_numerical_zero() << std::endl;
    std::cout << "Fast unsigned int range\t\t=\t" << "[0," << max_u << ']' << std::endl;
    std::cout << "Fast signed int range\t\t=\t" << '[' << min_i << ',' << max_i << ']' << std::endl;
    std::cout << "Default path\t=\t" << get_default_path() << std::endl;
    std::cout << "Piranha is ready." << std::endl;
    std::cout << "_______________________________" << std::endl << std::endl;
    // Setup cout.
    stream_manager::setup_print(std::cout);
    // Custom allocators for gmp.
    mp_set_memory_functions(mp_alloc,mp_realloc,mp_free);
#ifdef _PIRANHA_GSL
    // Turn off gsl error handler.
    gsl_set_error_handler_off();
#endif
#ifdef _PIRANHA_TBB
    std::cout << "TBB version strings:" << std::endl;
    std::cout << __TBB_VERSION_STRINGS;
    std::cout << "TBB:  BUILD_DATE\t\t" << __TBB_DATETIME << std::endl;
#endif
  }

  /// Set path to theories of motion.
  void settings_manager::set_path(const std::string &str)
  {
    path=str;
  }

  /// Get version.
  const std::string &settings_manager::get_version()
  {
    return version;
  }
}
