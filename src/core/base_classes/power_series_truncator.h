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

#include "../exceptions.h"
#include "../p_assert.h"
#include "degree_truncator.h"
#include "expo_truncator.h"
#include "norm_truncator.h"

namespace piranha
{
	/// Truncator for power series.
	class power_series_truncator
	{
		public:
			template <class Multiplier>
			class get_type:
				public expo_truncator::template get_type<Multiplier>,
				public degree_truncator::template get_type<Multiplier>,
				public norm_truncator::template get_type<Multiplier>
			{
					typedef typename expo_truncator::template get_type<Multiplier> expo_ancestor;
					typedef typename degree_truncator::template get_type<Multiplier> degree_ancestor;
					typedef typename norm_truncator::template get_type<Multiplier> norm_ancestor;
					enum selected_truncator {exp_t, deg_t, norm_t, null_t};
				public:
					typedef get_type type;
					get_type(Multiplier &m):expo_ancestor(m),degree_ancestor(m,false),
						norm_ancestor(m,false),m_active_truncator(exp_t) {
						if (!expo_ancestor::is_effective()) {
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
					}
					template <class T, class ArgsTuple>
					static size_t power_series_limit(const T &x, const ArgsTuple &args_tuple,
						const int &start = 0, const int &step_size = 1) {
						std::string msg("No useful truncation limit for a power series expansion could be "
							"established by the power series truncator. The reported errors where:\n");
						try {
							return expo_ancestor::power_series_limit(x,args_tuple,start,step_size);
						}
						catch (const base_exception &b) {
							msg += b.what() + "\n";
						}
						try {
							return degree_ancestor::power_series_limit(x,args_tuple,start,step_size);
						}
						catch (const base_exception &b) {
							msg += b.what() + "\n";
						}
						try {
							return norm_ancestor::power_series_limit(x,args_tuple,start,step_size);
						}
						catch (const base_exception &b) {
							msg += b.what() + "\n";
						}
						throw unsuitable(msg);
					}
					bool is_effective() const {
						return m_active_truncator != null_t;
					}
					template <class T>
					bool accept(const T &x) const {
						switch (m_active_truncator) {
							case exp_t:
								return expo_ancestor::accept(x);
							case deg_t:
								return degree_ancestor::accept(x);
							case norm_t:
								return norm_ancestor::accept(x);
							case null_t:
								p_assert(false);
						}
						p_assert(false);
						return true;
					}
					template <class T, class U>
					bool skip(const T &x1, const U &x2) const {
						switch (m_active_truncator) {
							case exp_t:
								return expo_ancestor::skip(x1,x2);
							case deg_t:
								return degree_ancestor::skip(x1,x2);
							case norm_t:
								return norm_ancestor::skip(x1,x2);
							case null_t:
								p_assert(false);
						}
						p_assert(false);
						return false;
					}
				private:
					selected_truncator m_active_truncator;
			};
	};
}

#endif
