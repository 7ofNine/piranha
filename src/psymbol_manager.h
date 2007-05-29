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

#ifndef PIRANHA_PSYMBOL_MANAGER_H
#define PIRANHA_PSYMBOL_MANAGER_H

#include <boost/thread/mutex.hpp>
#include <set>

#include "psymbol.h"

namespace piranha
  {
  /// Manage a global set of psymbols.
  /**
   * piranha::psymbol must be registered here to keep a global list of symbol for use in Poisson series.
   */
  class psymbol_manager
    {
    public:
      /// Functor used in psymbol comparison in set.
      struct ltpsymbol
        {
          size_t operator()(const psymbol &psym1, const psymbol &psym2) const
            {
              return (psym1.name().compare(psym2.name())<0);
            }
        };
      /// Alias for symbol set.
      typedef std::set
        <psymbol,ltpsymbol> set_type;
      /// Alias for standard iterator.
      typedef set_type::iterator iterator;
      static void print(std::ostream &stream=std::cout);
      /// Print to screen.
      static void put()
      {
        print();
      }
      static iterator reg(const psymbol &);
    private:
      static set_type             p_set_;
      static boost::mutex         mutex_;
    };


  /// Print to stream.
  inline void psymbol_manager::print(std::ostream &stream)
  {
    stream_manager::setup_print(stream);
    const iterator it_f=p_set_.end();
    for (iterator it1=p_set_.begin();it1!=it_f;++it1)
      {
        stream << "Symbol:" << std::endl;
        it1->print(stream);
      }
  }


  /// Register a symbol in the manager.
  inline psymbol_manager::iterator psymbol_manager::reg(const psymbol &psym)
  {
    // Guard with mutex, we could be registering symbols from more than one thread.
    boost::mutex::scoped_lock lock(mutex_)
      ;
    const iterator it=p_set_.find(psym);
    if (it==p_set_.end())
      {
        // Symbol is not already present, add it.
        std::pair<iterator,bool> result=p_set_.insert(psym);
        if (!result.second)
          {
            std::cout << "Damn!" << std::endl;
            std::exit(1);
          }
        return result.first;
      }
    else
      {
        // Symbol name is already registered, check that it is really the same symbol.
        if (psym!=*it)
          {
            std::cout << "Warning: you tried to add a psymbol with the same name of another one "
            << "but with different properties. The original one will be used." << std::endl;
            psym.print();
            std::cout << std::endl;
            it->print();
            std::cout << std::endl;
            std::exit(1);
          }
        return it;
      }
  }

  /// Typedefs used in series, terms, coefficients and trigonometric parts.
  typedef psymbol_manager::iterator psym_p;
  typedef std::vector<psym_p> vector_psym_p;
}

#endif
