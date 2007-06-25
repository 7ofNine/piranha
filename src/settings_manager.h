/***************************************************************************
 *   Copyright (C) 2007 by Francesco Biscani   *
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

#include <boost/thread/mutex.hpp>
#include <string>

namespace piranha
{
/// Manager class for piranha-specific settings.
/**
 * This class manages parameters specific to piranha classes.
 */
  class settings_manager
  {
    public:
/// Greeter class.
/**
 * Greeter class is constructed at piranha invocation and sets default parameters.
 */
      class greeter
      {
        public:
          greeter();
      };
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
/// Get default path to theories of motion files.
      static const std::string &default_theories_path()
      {
        return default_theories_path_;
      }
/// Get calculations precision.
      static const double &prec()
      {
        return prec_;
      }
      static void set_prec(const double &);
      static void set_theories_path(const std::string &);
    private:
/// Numerical zero.
      static double             numerical_zero_;
/// Jacobi Anger expansion limit.
      static const unsigned int jacang_limit_;
/// Path to theories of motion.
      static std::string        theories_path_;
      static const std::string  default_theories_path_;
/// Relative precision of series multiplications.
      static double             prec_;
      static greeter            grt_;
      static boost::mutex       mutex_;
  };
}
#endif                                            // PIRANHA_SETTINGS_H
