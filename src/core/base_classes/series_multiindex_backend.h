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

#ifndef PIRANHA_SERIES_MULTIINDEX_BACKEND_H
#define PIRANHA_SERIES_MULTIINDEX_BACKEND_H

#include <boost/multi_index_container.hpp>
#include <boost/tuple/tuple.hpp>
#include <memory>

#include "../settings.h"

namespace piranha
{
	// This assumes that second index is hashed.
	template < class Term, template <class> class I, class Allocator = std::allocator<Term> >
	class series_multiindex_backend
	{
			typedef Term term_type;
			typedef Allocator allocator_type;
			typedef typename I<term_type>::type index_type;
			typedef boost::multi_index_container<term_type, index_type, allocator_type> container_type;
		public:
			template <int N>
			class nth_index
			{
				public:
					typedef typename container_type::template nth_index<N>::type type;
			};
			series_multiindex_backend() {
				m_container.template get<1>().max_load_factor(settings::load_factor());
			}
			template <int N>
			typename nth_index<N>::type &get() {
				return m_container.template get<N>();
			}
			template <int N>
			const typename nth_index<N>::type &get() const {
				return m_container.template get<N>();
			}
			void swap(series_multiindex_backend &other) {
				m_container.swap(other.m_container);
			}
			static const int n_indices = boost::mpl::size<index_type>::value;
		private:
			container_type	m_container;
	};
}

#endif
