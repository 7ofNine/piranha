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

#ifndef PIRANHA_SERIES_UNORDERED_SET_BACKEND_H
#define PIRANHA_SERIES_UNORDERED_SET_BACKEND_H

#include <boost/static_assert.hpp>
#include <memory>
#include <utility> // For std::pair.
#include <tr1/unordered_set>

#include "../p_assert.h"
#include "../settings.h"

namespace piranha
{
	template < class Term, template <class> class, class Allocator = std::allocator<char> >
	class series_unordered_set_backend
	{
			typedef Term term_type;
			typedef Allocator allocator_type;
			typedef std::tr1::unordered_set < term_type, typename term_type::hasher,
			std::equal_to<term_type>, allocator_type > container_type;
		public:
			typedef typename container_type::iterator iterator;
			typedef typename container_type::const_iterator const_iterator;
			template <int N>
			class nth_index
			{
					BOOST_STATIC_ASSERT(N == 0 or N == 1);
				public:
					typedef series_unordered_set_backend type;
			};
			series_unordered_set_backend() {
				m_container.max_load_factor(settings::load_factor());
			}
			template <int N>
			series_unordered_set_backend &get() {
				BOOST_STATIC_ASSERT(N == 0 or N == 1);
				return *this;
			}
			template <int N>
			const series_unordered_set_backend &get() const {
				BOOST_STATIC_ASSERT(N == 0 or N == 1);
				return *this;
			}
			size_t size() const {
				return m_container.size();
			}
			const_iterator begin() const {
				return m_container.begin();
			}
			const_iterator end() const {
				return m_container.end();
			}
			iterator begin() {
				return m_container.begin();
			}
			iterator end() {
				return m_container.end();
			}
			bool empty() const {
				return m_container.empty();
			}
			void swap(series_unordered_set_backend &other) {
				m_container.swap(other.m_container);
			}
			iterator find(const term_type &t) {
				return m_container.find(t);
			}
			iterator insert(const const_iterator &, const term_type &t) {
				std::pair<iterator, bool> res(m_container.insert(t));
				p_assert(res.second);
				// TODO: understand here if it is better to return res.first or not.
				return res.first;
				//return m_container.end();
			}
			template <class Modifier>
			bool modify(iterator &it, Modifier &m) {
				m(*it);
				return true;
			}
			void erase(const const_iterator &it) {
				m_container.erase(*it);
			}
			static const int n_indices = 1;
		private:
			container_type	m_container;
	};
}

#endif
