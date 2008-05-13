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

#ifndef PIRANHA_ARG_MANAGER_H
#define PIRANHA_ARG_MANAGER_H

#include <boost/static_assert.hpp>
#include <boost/tuple/tuple.hpp>

#include "arg_manager_mp.h"
#include "ntuple.h"
#include "psym.h"

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
	template <class Term>
	class arg_manager
	{
			typedef typename ntuple<vector_psym_p &, 10>::type temp_tuple_type;
		public:
			/// Argument assigner.
			/**
			 * Create an instance of this class whenever there is need to refer to arguments but the context
			 * does not provide a method to pass such arguments.
			 */
			class arg_assigner
			{
				public:
					/// Constructor from tuple of arguments vectors.
					template <class ArgsTuple>
					arg_assigner(const ArgsTuple &args_tuple)
							: m_was_assigned(m_assigned) {
						p_assert(!m_was_assigned);
						arg_manager_assign_tuple<ArgsTuple, temp_tuple_type>::run(args_tuple, m_args_tuple);
						m_assigned = true;
					}
					/// Destructor: undoes assignment.
					~arg_assigner() {
						m_assigned = false;
					}
				private:
					bool m_was_assigned;
			};
			/// Check whether arguments were assigned or not.
			static bool assigned() {
				return m_assigned;
			}
			/// Retrieve const reference to arguments tuple.
			static const temp_tuple_type &get() {
				return m_args_tuple;
			}
		private:
			arg_manager() {}
			~arg_manager() {}
		private:
			static bool                 m_assigned;
			static vector_psym_p        m_arg_vector;
			static temp_tuple_type      m_args_tuple;
	};

	template <class Term>
	bool arg_manager<Term>::m_assigned = false;
	template <class Term>
	vector_psym_p arg_manager<Term>::m_arg_vector;
	template <class Term>
	typename arg_manager<Term>::temp_tuple_type arg_manager<Term>::m_args_tuple =
		typename arg_manager<Term>::temp_tuple_type(
			arg_manager<Term>::m_arg_vector,
			arg_manager<Term>::m_arg_vector,
			arg_manager<Term>::m_arg_vector,
			arg_manager<Term>::m_arg_vector,
			arg_manager<Term>::m_arg_vector,
			arg_manager<Term>::m_arg_vector,
			arg_manager<Term>::m_arg_vector,
			arg_manager<Term>::m_arg_vector,
			arg_manager<Term>::m_arg_vector,
			arg_manager<Term>::m_arg_vector
		);
}
#endif
