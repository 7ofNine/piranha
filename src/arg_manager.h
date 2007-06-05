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

#ifndef PIRANHA_ARG_MANAGER_H
#define PIRANHA_ARG_MANAGER_H

#include "psymbol.h"

namespace piranha
  {
  /// Manager for arguments.
  /**
   * This class is used to manage information about arguments in those context where such information
   * is not available. For instance when managing the elements of the multiindex container in a series
   * class the functors used to order the terms do not know anything about arguments, since they do not
   * appear inside each term. In those cases, before manipulating a multiindex container, an
   * arg_manager::arg_assigner instance should be created, so that proper arguments are made available
   * through the arg_manager::cf_args and arg_manager::trig_args methods.
   */
  class arg_manager
    {
    public:
      /// Argument assigner.
      /**
       * Create an instance of this class whenever there is need to refer to arguments but the context
       * does not provide a method to pass such arguments.
       */
      class arg_assigner
        {
        public:
          /// Constructor from pointers to piranha::vector_psym_p.
          arg_assigner(vector_psym_p const *cf_a, vector_psym_p const *trig_a)
              :was_assigned_(assigned_)
          {
            p_assert(!was_assigned_);
            cf_args_=cf_a;
            trig_args_=trig_a;
            assigned_=true;
          }
          /// Destructor: undoes assignment and unlocks resources.
          ~arg_assigner()
          {
            assigned_=false;
          }
        private:
          bool was_assigned_;
        }
      ;
      /// Check whether arguments were assigned or not.
      static bool assigned()
      {
        return assigned_;
      }
      /// Retrieve pointer to coefficient arguments vector.
      static vector_psym_p const *cf_args()
      {
        return cf_args_;
      }
      /// Retrieve pointer to trigonometric arguments vector.
      static vector_psym_p const *trig_args()
      {
        return trig_args_;
      }
    private:
      static bool                       assigned_;
      static boost::mutex               mutex_;
      static boost::mutex::scoped_lock  lock_;
      static vector_psym_p const        *cf_args_;
      static vector_psym_p const        *trig_args_;
    }
  ;
}

#endif
