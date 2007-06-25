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

#ifndef PIRANHA_PS_TERM_IO_H
#define PIRANHA_PS_TERM_IO_H

namespace piranha
{
/// Print in plain format.
  template <class Cf,class Trig>
    inline void ps_term<Cf,Trig>::print_plain(std::ostream &out_stream, const vector_psym_p &cv,
    const vector_psym_p &tv) const
  {
// Setup formatting.
    stream_manager::setup_print(out_stream);
    c_.print_plain(out_stream,cv);
    out_stream << stream_manager::data_separator();
    trig_args_.print_plain(out_stream,tv);
    switch (g_flavour())
    {
      case true:
        out_stream << "c";
        break;
      case false:
        out_stream << "s";
    }
  }

/// Print in latex format.
  template <class Cf,class Trig>
    inline void ps_term<Cf,Trig>::print_latex(std::ostream &out_stream, const vector_psym_p &cv,
    const vector_psym_p &tv) const
  {
// Setup formatting
    stream_manager::setup_print(out_stream);
    c_.print_latex(out_stream,cv);
    out_stream << "&";
    out_stream << "$" << phase(tv) << "$" << "&" << "$" << freq(tv) << "$" << "&";
    switch (g_flavour())
    {
      case true:
        out_stream << "c&";
        break;
      case false:
        out_stream << "s&";
    }
    trig_args_.print_latex(out_stream,tv);
  }
}
#endif
