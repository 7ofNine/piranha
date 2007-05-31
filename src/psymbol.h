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

#ifndef PIRANHA_PSYMBOL_H
#define PIRANHA_PSYMBOL_H

#include <boost/algorithm/string.hpp>
#include <boost/thread/mutex.hpp>
#include <cmath>
#include <set>
#include <string>
#include <vector>

#include "stream_manager.h"

namespace piranha
  {
  /// Simple symbol class
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
  * @see piranha::base_pseries::cf_s_vec_, the vector of symbols embedded in a Poisson series representing
  * its trigonometric arguments.
  * @see piranha::base_pseries::trig_s_vec_, the vector of symbols embedded in a Poisson series representing
  * its trigonometric arguments.
  */
  class psymbol
    {
    public:
      // Ctors
      /// Default constructor.
      /**
      * Assigns 'null' as name of the symbol.
      */
      //psymbol():name_(null_name_)
      //{}
      /// Constructor from std::string.
      /**
      * It just assigns psymbol::name_. psymbol::poly_eval_ will be empty.
      * @param[in] str symbol's name.
      */
      //psymbol(const std::string &str):name_(valid_name(str))
      //{}
      /// Constructor from std::string and std::vector<double>.
      /**
      * Assigns both psymbol::name_ and psymbol::poly_eval_.
      * @param[in] str symbol's name.
      * @param[in] pol symbol's polynomial evaluation vector.
      */
      //psymbol(const std::string &str, const std::vector<double> &pol):name_(valid_name(str)),poly_eval_(pol)
      //{}
      /// Copy constructor.
      psymbol(const psymbol &psym):name_(psym.name_),poly_eval_(psym.poly_eval_)
      {}
      /// Copy function used in python bindings.
      psymbol copy() const
        {
          return psymbol(*this);
        }
      bool operator==(const psymbol &) const;
      bool operator!=(const psymbol &psym) const
        {
          return !(*this==psym);
        }
      void print(std::ostream& out_stream=std::cout) const;
      /// Print to screen.
      void put() const
        {
          print(std::cout);
        }
      double eval (double) const;
      // Setters
      /// Set symbol's name.
      void set_name(const std::string &str)
      {
        name_=valid_name(str);
      }
      /// Set polynomial evaluation vector.
      void set_poly_eval(const std::vector<double> &pol)
      {
        poly_eval_=pol;
      }
      // Getters
      /// Get symbol's name.
      const std::string &name() const
        {
          return name_;
        }
      /// Get polynomial evaluation vector.
      const std::vector<double> &poly_eval() const
        {
          return poly_eval_;
        }
      /// Get symbol's phase.
      double phase() const
        {
          return get_poly_eval_elem(0);
        }
      /// Get symbol's frequency.
      double freq() const
        {
          return get_poly_eval_elem(1);
        }
      /// Get null name for symbols.
      static const std::string &null_name()
      {
        return null_name_;
      }
    private:
      // Helper for getting poly evals.
      double get_poly_eval_elem(const size_t &n) const
        {
          if (n>=poly_eval_.size())
            {
              std::cout << "WARNING: Poly eval element out of bounds, returning 0." << std::endl;
              return 0;
            }
          return poly_eval_[n];
        }
      // Helper for ctor.
      std::string valid_name(const std::string &str)
      {
        std::string tmp_str=str;
        boost::trim(tmp_str);
        if (tmp_str=="" || tmp_str.empty())
          {
            return std::string(null_name_);
          }
        else
          {
            return tmp_str;
          }
      }
      // Data members.
    private:
      std::string                 name_;
      std::vector<double>         poly_eval_;
      static const std::string    null_name_;
    };


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


  /// Print to stream.
  inline void psymbol::print(std::ostream &out_stream) const
    {
      stream_manager::setup_print(out_stream);
      out_stream << "name=" << name_ << std::endl;
      out_stream << "poly_eval=";
      for (unsigned int j=0;j<poly_eval_.size();++j)
        {
          out_stream << poly_eval_[j];
          if (j!=poly_eval_.size()-1)
            {
              out_stream << ';';
            }
        }
      out_stream << std::endl;
    }


  /// Test for equality.
  inline bool psymbol::operator==(const psymbol &psym2) const
    {
      if (name_!=psym2.name_)
        {
          std::cout << "Different name!" << std::endl;
          return false;
        }
      else if (poly_eval_!=psym2.poly_eval_)
        {
          std::cout << "Different poly_eval!" << std::endl;
          return false;
        }
      return true;
    }


  /// Time evaluation.
  //FIXME: rename to t_eval?
  inline double psymbol::eval(double t) const
    {
      double retval=0.;
      const size_t w=poly_eval_.size();
      for (size_t i=0;i<w;++i)
        {
          // FIXME: use natural_pow or the like here, to speed up?
          retval+=std::pow(t,(int)i)*poly_eval_[i];
        }
      return retval;
    }
}

#endif
