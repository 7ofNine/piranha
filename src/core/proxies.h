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

#ifndef PIRANHA_PROXIES_H
#define PIRANHA_PROXIES_H

/*! \file proxies.h
    \brief Proxies.

    Classes implementing the proxy pattern.
*/

namespace piranha
{
	template <class T>
	struct copy_proxy {
		copy_proxy(): m_element() {}
		const T &get() const {
			return m_element;
		}
		void assignment(const T &x) {
			m_element = x;
		}
protected:
		T m_element;
	};

	template <class T>
	struct reference_proxy {
		reference_proxy(): m_element(0) {}
		const T &get() const {
			p_assert(m_element != 0);
			return *m_element;
		}
		void assignment(const T &x) {
			m_element = &x;
		}
protected:
		T const *m_element;
	};

	/// Proxy for coefficients during series multiplication.
	/**
	 * Defaults to copying. Use template specialization to change the behaviour.
	 */
	template <class Cf>
	class cf_mult_proxy: public copy_proxy<Cf>
	{
			typedef copy_proxy<Cf> ancestor;
		public:
			cf_mult_proxy(): ancestor() {}
			void operator=(const Cf &cf) {
				ancestor::assignment(cf);
			}
	};

	/// Proxy for keys during series multiplication.
	/**
	 * Defaults to referencing. Use template specialization to change the behaviour.
	 */
	template <class Key>
	class key_mult_proxy: public reference_proxy<Key>
	{
			typedef reference_proxy<Key> ancestor;
		public:
			key_mult_proxy(): ancestor() {}
			void operator=(const Key &key) {
				ancestor::assignment(key);
			}
	};
}

#endif
