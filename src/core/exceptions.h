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

#ifndef PIRANHA_EXCEPTIONS_H
#define PIRANHA_EXCEPTIONS_H

#include <string>

namespace piranha
{
  class base_exception
  {
    public:
      base_exception(const std::string &s):m_what(s) {}
      const std::string &what() const
      {
        return m_what;
      }
    private:
      std::string m_what;
  };

  struct bad_input:public base_exception
  {
    bad_input(const std::string &s):base_exception(s) {}
  };

  struct term_not_insertable:public base_exception
  {
    term_not_insertable(const std::string &s):base_exception(s) {}
  };

  struct unsuitable:public base_exception
  {
    unsuitable(const std::string &s):base_exception(s) {}
  };

  struct not_existing:public base_exception
  {
    not_existing(const std::string &s):base_exception(s) {}
  };
}

#endif
