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

#include "settings_manager.h"
#include "stream_manager.h"
#include "version.h"

namespace piranha
  {
  // Settings manager static members.
  double settings_manager::numerical_zero_ = 1E-80;
  const unsigned int settings_manager::jacang_limit_ = 20;
  double settings_manager::prec_ = 1E-6;
  std::string settings_manager::theories_path_ = _PIRANHA_THEORIES_DIR;
  const std::string settings_manager::default_theories_path_ = _PIRANHA_THEORIES_DIR;
  settings_manager::greeter settings_manager::grt_;
  boost::mutex settings_manager::mutex_;

  settings_manager::greeter::greeter()
  {
    // Startup report.
    std::cout << "This is Piranha version " << __piranha_version << std::endl;
    std::cout << "Default parameters initialized:" << std::endl;
    std::cout << "Print precision\t\t\t=\t" << stream_manager::digits() << std::endl;
    std::cout << "Numerical zero\t\t\t=\t" << numerical_zero() << std::endl;
    std::cout << "Multiplication precision\t=\t" << prec() << std::endl;
    std::cout << "Theories of motion directory\t=\t" << theories_path() << std::endl;
    std::cout << "Piranha is ready." << std::endl;
    std::cout << "_______________________________" << std::endl << std::endl;
  }


  /// Set precision in multiplication of series.
  void settings_manager::set_prec(const double &value)
  {
    boost::mutex::scoped_lock lock(mutex_)
      ;
    prec_=value;
  }


  /// Set path to theories of motion.
  void settings_manager::set_theories_path(const std::string &str)
  {
    boost::mutex::scoped_lock lock(mutex_)
      ;
    theories_path_=str;
  }
}
