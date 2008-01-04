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

#ifndef PIRANHA_PHASE_LIST_H
#define PIRANHA_PHASE_LIST_H

#include <boost/algorithm/string.hpp>             // to_lower FIXME: why is it used here??
                                                  // Check those file open routines, uniform them.

#include "../common_typedefs.h"                      // vector_double.
#include "../utils.h"                                // open_file.

namespace piranha
{
/// List of phases
/**
 * A list of real phases, to be read by a file (each line is a phase). This is used in
 * series to insert phases into series' terms with piranha::base_pseries::insert_phases.
 * Phases are stored inside the class in a STL vector.
 * For each phase \f$ \phi_i \f$ inserted a term is modified and a new one is created. For instance:
 * \f[
 * C\cos\left(x+\phi_i\right)=C\cos x \cos \phi_i - C \sin x \sin \phi_i.
 * \f]
 * Phases can either be added or assigned. If they are added \f$ \phi_i \f$ is simply an element of the
 * vector. If phases are assigned, the phase of each term is calculated (using the properties of each
 * argument) and subtracted from the corresponding phase in the vector. The difference between the two
 * becomes \f$ \phi_i \f$, so that the final effect is the same as if the arguments' combined phases
 * were originally the ones provided in the phase list.
 * Phase assignment was created with TASS in mind, since TASS provides the phases of each term which are
 * not necessarily the phases deriving from the combination of arguments. In other words, arbitrary
 * phases are present and not mentioned explicitly, so we need to infer them using this method.
 *
 * Whether phases are added or assigned is establised by the piranha::phase_list::operation_ parameter.
 * @see piranha::base_pseries::insert_phases.
 */
  class phase_list
  {
    public:
/// Types of operations that can be performed with phases.
      enum op
      {
        add
        ,
        assign
      };
      typedef vector_double::const_iterator const_iterator;
      typedef vector_double::const_reverse_iterator const_reverse_iterator;
      phase_list():operation_(add) {}
      phase_list(const std::string &);
/// Assignment operator.
      phase_list &operator=(const phase_list &pl)
      {
        if (&pl!=this)
        {
          phases_=pl.phases_;
          operation_=pl.operation_;
        }
        return *this;
      }
/// Return begin iterator.
      const_iterator begin() const
      {
        return phases_.begin();
      }
/// Return end iterator.
      const_iterator end() const
      {
        return phases_.end();
      }
/// Return reverse begin iterator.
      const_reverse_iterator rbegin() const
      {
        return phases_.rbegin();
      }
/// Return reverse end iterator.
      const_reverse_iterator rend() const
      {
        return phases_.rend();
      }
/// Get the type of operation to be performed with phases.
      op operation() const
      {
        return operation_;
      }
/// Print phase_list to screen.
      void put() const
      {
        std::cout << operation_ << std::endl;
        for (const_iterator it=begin();it!=end();++it)
        {
          std::cout << *it << std::endl;
        }
      }
    private:
/// List of phases.
      vector_double    phases_;
/// Operation to be performed with phases.
      op          operation_;
  };

/// Constructor from file.
/**
 * Read a file line by line, converting each line in a phase. The first optional line of the file
 * is a string which specifies whether the phases have to be added ('add') or assigned ('assign').
 * If no such string is found the default is to add phases.
 * @param[in] fn std::string representing filename.
 */
  inline phase_list::phase_list(const std::string &fn):operation_(add)
  {
    std::ifstream inf;
    utils::open_file(fn,inf);
// Read from file
    if (inf.is_open())
    {
      std::string temp;
      double tmp_phase;
      bool op_read=false;
      while (utils::get_valid_string(inf,temp))
      {
        try
        {
          tmp_phase=boost::lexical_cast<double>(temp);
        }
        catch(boost::bad_lexical_cast &)
        {
          if (op_read)
          {
            std::cout << "Error reading phase from file, returning 0." << std::endl;
            tmp_phase=0.;
          }
          else
          {
            boost::algorithm::to_lower(temp);
            if (temp=="add")
            {
              operation_=add
                ;
            }
            else if (temp=="assign")
            {
              operation_=assign;
            }
            else
            {
              std::cout << "Unrecoginzed operation type, defaulting to 'add'." << std::endl;
            }
            op_read=true;
            continue;
          }
        }
        phases_.push_back(tmp_phase);
      }
// Close file
      inf.close();
    }
  }
}
#endif
