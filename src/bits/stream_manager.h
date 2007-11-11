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

#ifndef PIRANHA_STREAM_MANAGER_H
#define PIRANHA_STREAM_MANAGER_H

#include <iomanip>
#include <iostream>
#include <string>

namespace piranha
{
  class stream_manager
  {
    public:
      enum out_format
      {
        plain,
        latex
      };
      enum fp_representation
      {
        scientific,
        decimal
      };
// Getters.
      static unsigned int digits();
      static const std::string &data_separator();
      static unsigned int min_digits();
      static unsigned int max_digits();
// Setters
      static void set_digits(int n);
      static void setup_print(std::ostream &out_stream=std::cout);
      static out_format format();
      static void set_format(out_format);
      static fp_representation fp_rep();
      static void set_fp_rep(fp_representation);
    private:
/// Minimum number of digits for output streams.
      static const unsigned int       min_digits_;
/// Maximum number of digits for output streams.
      static const unsigned int       max_digits_;
/// Number of digits to display in output stream.
      static unsigned int             digits_;
/// Data separator for formatted I/O.
      static const std::string        data_separator_;
/// Format for output.
      static out_format               format_;
/// Floating point representation.
      static fp_representation        fp_rep_;
  };
}
#endif
