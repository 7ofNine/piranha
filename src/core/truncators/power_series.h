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

#ifndef PIRANHA_POWER_SERIES_TRUNCATOR
#define PIRANHA_POWER_SERIES_TRUNCATOR

#include <string>
#include <vector>

#include "../exceptions.h"
#include "degree.h"
#include "norm.h"

namespace piranha
{
namespace truncators
{
	/// Truncator for power series.
	class power_series
	{
		public:
			template <class Multiplier>
			class get_type:
				public degree::template get_type<Multiplier>,
				public norm::template get_type<Multiplier>
			{
					typedef typename degree::template get_type<Multiplier> degree_ancestor;
					typedef typename norm::template get_type<Multiplier> norm_ancestor;
					enum selected_truncator {deg_t, norm_t, null_t};
				public:
					typedef get_type type;
					get_type(Multiplier &m):degree_ancestor(m,false),
						norm_ancestor(m,false),m_active_truncator(deg_t) {
						degree_ancestor::init();
						if (!degree_ancestor::is_effective()) {
							norm_ancestor::init();
							if (!norm_ancestor::is_effective()) {
								m_active_truncator = null_t;
							} else {
								m_active_truncator = norm_t;
							}
						} else {
							m_active_truncator = deg_t;
						}
					}
					template <class T, class ArgsTuple>
					static size_t power_series_iterations(const T &x, const int &start, const int &step_size,
						const ArgsTuple &args_tuple) {
						std::string msg("No useful truncation limit for a power series expansion could be "
							"established by the power series truncator. The reported errors were:\n");
						try {
							return degree_ancestor::power_series_iterations(x,start,step_size,args_tuple);
						}
						catch (const value_error &ve) {
							msg += std::string(ve.what()) + "\n";
						}
						try {
							return norm_ancestor::power_series_iterations(x,start,step_size,args_tuple);
						}
						catch (const value_error &ve) {
							msg += std::string(ve.what()) + "\n";
						}
						piranha_throw(value_error,msg);
					}
					template <class Series, class ArgsTuple>
					static std::vector<typename Series::term_type const *> get_sorted_pointer_vector(const Series &s, const ArgsTuple &args_tuple)
					{
						std::string msg("The power series truncator was not able to establish a series ordering. The reported errors were:\n");
						try {
							return degree_ancestor::get_sorted_pointer_vector(s,args_tuple);
						} catch (const value_error &ve) {
									msg += std::string(ve.what()) + "\n";
						}
						try {
							return norm_ancestor::get_sorted_pointer_vector(s,args_tuple);
						} catch (const value_error &ve) {
									msg += std::string(ve.what()) + "\n";
						}
						piranha_throw(value_error,msg);
					}
					bool is_effective() const {
						return m_active_truncator != null_t;
					}
					template <class T>
					bool accept(const T &x) const {
						switch (m_active_truncator) {
							case deg_t:
								return degree_ancestor::accept(x);
							case norm_t:
								return norm_ancestor::accept(x);
							case null_t:
								piranha_assert(false);
						}
						piranha_assert(false);
						return true;
					}
					template <class T, class U>
					bool skip(const T &x1, const U &x2) const {
						switch (m_active_truncator) {
							case deg_t:
								return degree_ancestor::skip(x1,x2);
							case norm_t:
								return norm_ancestor::skip(x1,x2);
							case null_t:
								piranha_assert(false);
						}
						piranha_assert(false);
						return false;
					}
				private:
					selected_truncator m_active_truncator;
			};
	};
}
}

#endif