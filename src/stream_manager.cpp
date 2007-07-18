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

#include "stream_manager.h"

namespace piranha
{
  unsigned int stream_manager::digits_ = 15;
  const std::string stream_manager::data_separator_ = ";";
  stream_manager::out_format stream_manager::format_ = stream_manager::plain;
  const unsigned int stream_manager::min_digits_ = 0;
  const unsigned int stream_manager::max_digits_ = 50;
  stream_manager::fp_representation stream_manager::fp_rep_ = stream_manager::scientific;

  unsigned int stream_manager::digits()
  {
    return digits_;
  }

  unsigned int stream_manager::min_digits()
  {
    return min_digits_;
  }

  unsigned int stream_manager::max_digits()
  {
    return max_digits_;
  }

  const std::string &stream_manager::data_separator()
  {
    return data_separator_;
  }

  void stream_manager::set_digits(int n)
  {
    if (n<(int)min_digits_ || n>(int)max_digits_)
    {
      std::cout << "Invalid number of digits." << std::endl;
    }
    else
    {
      digits_=(unsigned int)n;
    }
  }

  void stream_manager::set_fp_rep(fp_representation fpr)
  {
    fp_rep_=fpr;
  }

  stream_manager::fp_representation stream_manager::fp_rep()
  {
    return fp_rep_;
  }

  void stream_manager::setup_print(std::ostream &out_stream)
  {
    out_stream << std::setprecision(digits_);
    switch (fp_rep_)
    {
      case scientific:
        out_stream << std::scientific;
        break;
      case decimal:
        out_stream << std::fixed;
        break;
    }
  }

  stream_manager::out_format stream_manager::format()
  {
    return format_;
  }

  void stream_manager::set_format(out_format fmt)
  {
    format_=fmt;
  }
}
