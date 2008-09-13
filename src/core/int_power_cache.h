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

#ifndef PIRANHA_INT_POWER_CACHE_H
#define PIRANHA_INT_POWER_CACHE_H

#include <boost/unordered_map.hpp>
#include <vector>

#include "integer_typedefs.h"
#include "p_assert.h"
#include "settings.h"

namespace piranha
{
	// Uses exponentiation by squaring (EBS) internally. Computed values are stored in a hash
	// map and re-used to calculate new values. To use it, simply construct and request the
	// value with operator[].
	template <class T, class ArithmeticFunctor>
	class int_power_cache
	{
		protected:
			typedef boost::unordered_map<max_fast_int,T> container_type;
			typedef typename container_type::iterator iterator;
			typedef ArithmeticFunctor arith_functor_type;
			int_power_cache():m_container(),m_arith_functor()
			{}
		public:
			int_power_cache(const T &x_1):
				m_container(),
				m_arith_functor() {
				init(x_1);
			}
			int_power_cache(const T &x_1, const T &inv_x):
				m_container(),
				m_arith_functor() {
				init(x_1);
				m_container[-1] = inv_x;
			}
			const T &operator[](const max_fast_int &n) {
				iterator it = m_container.find(n);
				if (it == m_container.end()) {
					p_assert(n != 0 && n != 1);
					if (n > 0) {
						return insert_new<true>(n);
					} else {
						return insert_new<false>(n);
					}
				} else {
					return it->second;
				}
			}
		private:
			void init(const T &x_1) {
				m_container[0] = T(static_cast<max_fast_int>(1));
				m_container[1] = x_1;
			}
			template <bool Sign>
			T &insert_new(max_fast_int n) {
				p_assert((Sign && n > 0) || (!Sign && n < 0));
				const max_fast_int orig_n = n;
				std::vector<max_fast_int> ebs_sequence;
				// First we find the first available value in the EBS sequence.
				const iterator it_f = m_container.end();
				__PDEBUG(std::cout << "Int power cache pushing " << orig_n << '\n');
				ebs_sequence.push_back(orig_n);
				do {
					// n == -1 may not have been provided in the ctor. If we reached this point, it means
					// that our negative EBS sequence has come to the end.
					if (n == -1) {
						m_container[-1] = m_arith_functor.pow(m_container[1],static_cast<max_fast_int>(-1));
						break;
					}
					if (n & 1) {
						if (Sign) {
							--n;
						} else {
							++n;
						}
					} else {
						n >>= 1;
					}
					__PDEBUG(std::cout << "Int power cache pushing " << n << '\n');
					ebs_sequence.push_back(n);
				} while (m_container.find(n) == it_f);
				// At this point n is the lowest known integer of our EBS sequence. It is a kind of seed value.
				// We are now going to pop back the sequence and multiply by the last popped element or x_1
				// at each iteration.
				max_fast_int index = ebs_sequence.back();
				while (index != orig_n) {
					// Check that we are not going to pop too much.
					p_assert(ebs_sequence.size() >= 2);
					max_fast_int tmp = ebs_sequence.back();
					__PDEBUG(std::cout << "Int power cache popping " << ebs_sequence.back() << '\n');
					ebs_sequence.pop_back();
					index = ebs_sequence.back();
					p_assert((Sign && index > tmp) || (!Sign && index < tmp));
					m_container[index] = m_container[tmp];
					if (std::abs(index - tmp) == 1) {
						if (Sign) {
							m_arith_functor.multiply(m_container[index], m_container[1]);
						} else {
							m_arith_functor.multiply(m_container[index], m_container[-1]);
						}
					} else {
						p_assert(index % tmp == 0 && std::abs(index / tmp) == 2)
						m_arith_functor.multiply(m_container[index], m_container[tmp]);
					}
				}
				return m_container[index];
			}
		protected:
			container_type				m_container;
			const arith_functor_type	m_arith_functor;
	};
}

#endif
