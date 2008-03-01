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
      static unsigned int min_digits();
      static unsigned int max_digits();
      // Setters
      static void set_digits(int n);
      static void setup_print(std::ostream &out_stream=std::cout)
      {
        out_stream << std::setprecision(m_digits);
        switch (m_fp_rep)
        {
          case scientific:
            out_stream << std::scientific;
            break;
          case decimal:
            out_stream << std::fixed;
            break;
        }
      }
      static out_format format();
      static void set_format(out_format);
      static fp_representation fp_rep();
      static void set_fp_rep(fp_representation);
    private:
      /// Minimum number of digits for output streams.
      static const unsigned int       m_min_digits;
      /// Maximum number of digits for output streams.
      static const unsigned int       m_max_digits;
      /// Number of digits to display in output stream.
      static unsigned int             m_digits;
      /// Format for output.
      static out_format               m_format;
      /// Floating point representation.
      static fp_representation        m_fp_rep;
  };

  // Static initializations.
  unsigned int stream_manager::m_digits = 15;
  stream_manager::out_format stream_manager::m_format = stream_manager::plain;
  stream_manager::fp_representation stream_manager::m_fp_rep = stream_manager::scientific;
}

#endif
