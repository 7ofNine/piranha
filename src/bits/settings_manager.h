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

#ifndef PIRANHA_SETTINGS_MANAGER_H
#define PIRANHA_SETTINGS_MANAGER_H

#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
#include <string>

#include "common_typedefs.h"
#include "compile_switches.h"
#include "piranha_tbb.h" // For task scheduler init.

namespace piranha
{
/// Manager class for piranha-specific settings.
/**
 * This class manages parameters specific to piranha classes.
 */
  class settings_manager
  {
    public:
/// Return maximum load factor for hashed containers.
/**
 * Relevant only to those containers with support such setting.
 */
      static const double &load_factor()
      {
        return hash_max_load_factor_;
      }
// Getters.
/// Get Jacobi-Anger expansion limit.
      static unsigned int jacang_lim()
      {
        return jacang_limit_;
      }
/// Get numerical zero.
      static const double &numerical_zero()
      {
        return numerical_zero_;
      }
/// Get path to theories of motion files.
      static const std::string &theories_path()
      {
        return theories_path_;
      }
/// Get Piranha version.
      static const std::string &version();
/// Get default path to theories of motion files.
      static const std::string &default_theories_path()
      {
        return default_theories_path_;
      }
/// Set maximum load factor for hashed containers.
/**
 * @see settings_manager::load_factor().
 */
      static void set_load_factor(const double &value)
      {
        if (value <= 0 or value >= 1)
        {
          std::cout << "Please insert a real number in the ]0,1[ interval." << std::endl;
          return;
        }
        hash_max_load_factor_ = value;
      }
      static void set_theories_path(const std::string &);
      static bool display_progress()
      {
        return enable_progress_display_;
      }
      static void set_display_progress(bool flag)
      {
        if (!(compile_switches::display_progress))
        {
          std::cout << "Progress bar has not been built." << std::endl;
        }
        enable_progress_display_=flag;
      }
      static unsigned long int mp_default_prec()
      {
        return mpf_get_default_prec();
      }
      static void set_mp_default_prec(int n)
      {
        if (n <= 0)
        {
          std::cout << "Please insert a strictly positive value." << std::endl;
          return;
        }
        mpf_set_default_prec(n);
      }
    private:
/// Private ctor.
      settings_manager() {}
/// Startup class.
/**
 * Startup class is constructed at piranha invocation and sets default parameters.
 */
      class startup
      {
        public:
          startup();
      };
/// Load factor for hashed containers.
      static double                         hash_max_load_factor_;
/// Numerical zero.
      static double                         numerical_zero_;
/// Minimum fast unsigned integer.
      static const max_fast_uint            min_u_;
/// Maximum fast unsigned integer.
      static const max_fast_uint            max_u_;
/// Minimum fast integer.
      static const max_fast_int             min_i_;
/// Maximum fast integer.
      static const max_fast_int             max_i_;
/// Jacobi Anger expansion limit.
      static const unsigned int             jacang_limit_;
/// Path to theories of motion.
      static std::string                    theories_path_;
      static const std::string              default_theories_path_;
      static const std::string              version_;
      static bool                           enable_progress_display_;
      static startup                        startup_;
#ifdef _PIRANHA_TBB
      static const tbb::task_scheduler_init tbb_init;
#endif
  };
}

#endif
