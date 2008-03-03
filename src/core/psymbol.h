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

#ifndef PIRANHA_PSYMBOL_H
#define PIRANHA_PSYMBOL_H

#include <boost/algorithm/string.hpp>
#include <boost/array.hpp>
#include <cmath>
#include <set>
#include <string>
#include <vector>

#include "common_typedefs.h"
#include "p_assert.h"
#include "stream_manager.h"
#include "utils.h"

namespace piranha
{
  /// Manage a global set of psymbols.
  /**
   * piranha::psymbol must be registered here to keep a global list of symbol for use in Poisson series.
   */
  class psymbol_manager
  {
    public:
      /// Simple symbol class.
      /**
       * This class is used in piranha to represent symbolic entitites (arguments and, in a later stage,
       * symbolic coefficients). It features a string representing the symbol's name and a numerical vector
       * which is used to evaluate the symbol in a polynomial fashion. For instance, if the numerical vector
       * has a size of three and its elements are named \f$ \alpha \f$, \f$ \beta \f$ and \f$ \gamma \f$,
       * it means that the symbol is evaluated as
       * \f[
       * \alpha + \beta t + \gamma t^2,
       * \f]
       * where \f$ t \f$ is time.
       * @see piranha::base_pseries::m_arguments tuple of arguments of a Poisson series.
       */
      class psymbol
      {
        public:
          // Ctors
          /// Constructor from std::string and vector of doubles.
          /**
           * Assigns both psymbol_manager::psymbol::m_name and psymbol_manager::psymbol::m_time_eval.
           * @param[in] str symbol's name.
           * @param[in] pol symbol's polynomial evaluation vector.
           */
          psymbol(const std::string &name, const std::string &te):m_name(name),
            m_time_eval(utils::str_to_vector<double>(te,separator))
            {psymbol_manager::reg(*this);}
          psymbol(const std::string &);
          psymbol(const std::string &s, const double &x1):m_name(s),m_time_eval((size_t)1)
          {
            boost::array<double,1> tmp =
            {
              {
                x1
              }
            };
            build_from_array(tmp);
          }
          psymbol(const std::string &s, const double &x1, const double &x2):m_name(s),m_time_eval((size_t)2)
          {
            boost::array<double,2> tmp =
            {
              {
                x1, x2
              }
            };
            build_from_array(tmp);
          }
          psymbol(const std::string &s, const double &x1, const double &x2, const double &x3):
          m_name(s),m_time_eval((size_t)3)
          {
            boost::array<double,3> tmp =
            {
              {
                x1, x2, x3
              }
            };
            build_from_array(tmp);
          }
          psymbol(const std::string &s, const double &x1, const double &x2, const double &x3, const double &x4):
          m_name(s),m_time_eval((size_t)4)
          {
            boost::array<double,4> tmp =
            {
              {
                x1, x2, x3, x4
              }
            };
            build_from_array(tmp);
          }
          psymbol(const std::string &s, const double &x1, const double &x2, const double &x3,
            const double &x4, const double &x5):m_name(s),m_time_eval((size_t)5)
          {
            boost::array<double,5> tmp =
            {
              {
                x1, x2, x3, x4, x5
              }
            };
            build_from_array(tmp);
          }
          /// Copy function used in python bindings.
          psymbol copy() const
          {
            return psymbol(*this);
          }
          bool operator==(const psymbol &) const;
          void print(std::ostream& out_stream=std::cout) const;
          /// Print to screen.
          void put() const {print(std::cout);}
          double t_eval(const double &) const;
          // Getters
          /// Get symbol's name.
          const std::string &name() const {return m_name;}
          /// Get polynomial evaluation vector.
          const std::vector<double> &time_eval() const {return m_time_eval;}
          /// Get symbol's phase.
          double phase() const {return get_time_eval<0>();}
          /// Get symbol's frequency.
          double freq() const {return get_time_eval<1>();}
          std::string powers_string() const {return utils::vector_to_str(m_time_eval,psymbol::separator);}
        private:
          // Default constructor. Privatized because it mustn't be called.
          psymbol() {p_assert(false);}
          // Helper for getting time evals.
          template <int N>
          double get_time_eval() const
          {
            if (N >= m_time_eval.size())
            {
              std::cout << "WARNING: Time eval element out of bounds, returning 0." << std::endl;
              return 0;
            }
            return m_time_eval[N];
          }
          // Helper for ctor from boost::array.
          template <class T>
            void build_from_array(const T &);
          // Data members.
        private:
          std::string                 m_name;
          std::vector<double>         m_time_eval;
          static const std::string    separator;
      };
      /// Functor used in psymbol comparison in set.
      struct ltpsymbol
      {
        size_t operator()(const psymbol &psym1, const psymbol &psym2) const
        {
          return (psym1.name().compare(psym2.name())<0);
        }
      };
      /// Alias for symbol set.
      // TODO: turn into multiindex container, so that we can search by name and drop that ugly O(n)
      // operation somewhere below.
      typedef std::set
        <psymbol,ltpsymbol> set_type;
      /// Alias for standard iterator, to be used in pyranha.
      typedef set_type::const_iterator iterator;
      /// Alias for iterator.
      typedef iterator psym_p;
      /// Print to stream.
      static void print(std::ostream &stream=std::cout)
      {
        stream_manager::setup_print(stream);
        const iterator it_f=p_set_.end();
        for (iterator it=p_set_.begin();it!=it_f;++it)
        {
          stream << "Symbol:" << std::endl;
          it->print(stream);
        }
      }
      /// Print to screen.
      static void put() {print();}
      static std::pair<bool,iterator> get_pointer(const psymbol &psym)
      {
        std::pair<bool,iterator> retval=get_pointer(psym.name());
        p_assert(retval.first);
        return retval;
      }
      // NOTICE: this is an O(n) operation, maybe it can be sped up.
      static std::pair<bool,iterator> get_pointer(const std::string &name)
      {
        std::pair<bool,iterator> retval(false,end());
        for (iterator it=begin();it!=retval.second;++it)
        {
          if (it->name() == name)
          {
            retval.first=true;
            retval.second=it;
            break;
          }
        }
        return retval;
      }
      // Function to iterate over the manager. To be used in pyranha.
      static iterator begin() {return p_set_.begin();}
      static iterator end() {return p_set_.end();}
      static size_t length() {return p_set_.size();}
    private:
      /// Register a symbol.
      /**
       * Searches for symbol: if it is not found, inserts new, otherwise will override its properties.
       */
      static void reg(const psymbol &psym)
      {
        const iterator it=p_set_.find(psym);
        if (it == p_set_.end())
        {
          // Symbol is not already present, add it.
          std::pair<iterator,bool> result=p_set_.insert(psym);
          p_assert(result.second);
        }
        else
        {
          // Symbol name is already registered, check that it is really the same symbol.
          if (!(psym == *it))
          {
            // TODO: overwrite here instead of using old.
            std::cout << "Warning: you tried to add a psymbol with the same name of another one "
              << "but with different properties. The original one will be used." << std::endl;
            psym.print();
            std::cout << std::endl;
            it->print();
            std::cout << std::endl;
            std::exit(1);
          }
        }
      }
      // Data members.
    private:
      static set_type             p_set_;
  };

  /// Constructor from std::string.
  /**
   * Searches for psymbol in piranha::psymbol_manager. If found, it builds a copy of symbol, otherwise
   * psymbol with zero m_time_eval is built.
   */
  inline psymbol_manager::psymbol::psymbol(const std::string &str):m_name(str),m_time_eval()
  {
    std::pair<bool,psym_p> res=psymbol_manager::get_pointer(str);
    if (!(res.first))
    {
      psymbol_manager::reg(*this);
    }
    else
    {
      m_name=res.second->name();
      m_time_eval=res.second->time_eval();
    }
  }

  // Helper for ctor from boost::array.
  template <class T>
    inline void psymbol_manager::psymbol::build_from_array(const T &a)
  {
    p_assert(a.size()==m_time_eval.size());
    for (size_t i=0;i<a.size();++i)
    {
      m_time_eval[i]=a[i];
    }
    psymbol_manager::reg(*this);
  }

  /// Print to stream.
  inline void psymbol_manager::psymbol::print(std::ostream &out_stream) const
  {
    stream_manager::setup_print(out_stream);
    out_stream << "name=" << m_name << std::endl;
    out_stream << "time_eval=";
    for (unsigned int j=0;j<m_time_eval.size();++j)
    {
      out_stream << m_time_eval[j];
      if (j != m_time_eval.size()-1)
      {
        out_stream << separator;
      }
    }
    out_stream << std::endl;
  }

  /// Test for equality.
  inline bool psymbol_manager::psymbol::operator==(const psymbol &psym2) const
  {
    if (m_name != psym2.m_name)
    {
      std::cout << "Different name!" << std::endl;
      return false;
    }
    else if (m_time_eval != psym2.m_time_eval)
    {
      std::cout << "Different time_eval!" << std::endl;
      return false;
    }
    return true;
  }

  /// Time evaluation.
  inline double psymbol_manager::psymbol::t_eval(const double &t) const
  {
    double retval=0.;
    const size_t w=m_time_eval.size();
    for (size_t i=0;i<w;++i)
    {
      // FIXME: use natural_pow or the like here, to speed up?
      retval+=std::pow(t,(int)i)*m_time_eval[i];
    }
    return retval;
  }

  /// Typedefs used in series, terms, coefficients and trigonometric parts.
  typedef psymbol_manager::psym_p psym_p;
  typedef std::vector<psym_p> vector_psym_p;
  typedef psymbol_manager::psymbol psymbol;
}
#endif
