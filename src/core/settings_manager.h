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
      static const double &get_load_factor()
      {
        return hash_max_load_factor;
      }
      // Getters.
      /// Get Jacobi-Anger expansion limit.
      static unsigned int jacang_lim()
      {
        return jacang_limit;
      }
      /// Get numerical zero.
      static const double &get_numerical_zero()
      {
        return numerical_zero;
      }
      /// Get path to theories of motion files.
      static const std::string &get_path()
      {
        return path;
      }
      /// Get Piranha version.
      static const std::string &get_version();
      /// Get default path to theories of motion files.
      static const std::string &get_default_path()
      {
        return default_path;
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
        hash_max_load_factor = value;
      }
      static void set_path(const std::string &);
      static bool display_progress()
      {
        return enable_progress_display;
      }
      static void set_display_progress(bool flag)
      {
        if (!(compile_switches::display_progress))
        {
          std::cout << "Warning: progress bar has not been built." << std::endl;
        }
        enable_progress_display=flag;
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
      class startup_class
      {
        public:
          startup_class();
      };
      /// Load factor for hashed containers.
      static double                         hash_max_load_factor;
      /// Numerical zero.
      static double                         numerical_zero;
      /// Minimum fast unsigned integer.
      static const max_fast_uint            min_u;
      /// Maximum fast unsigned integer.
      static const max_fast_uint            max_u;
      /// Minimum fast integer.
      static const max_fast_int             min_i;
      /// Maximum fast integer.
      static const max_fast_int             max_i;
      /// Jacobi Anger expansion limit.
      static const unsigned int             jacang_limit;
      /// Path to theories of motion.
      static std::string                    path;
      static const std::string              default_path;
      static const std::string              version;
      static bool                           enable_progress_display;
      static startup_class                  startup;
#ifdef _PIRANHA_TBB
      static const tbb::task_scheduler_init tbb_init;
#endif
  };
}
#endif
