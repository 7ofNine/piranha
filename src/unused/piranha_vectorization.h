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

#ifndef PIRANHA_VECTORIZATION_H
#define PIRANHA_VECTORIZATION_H

#include "config.h"
#include "integer_typedefs.h"

#ifdef _PIRANHA_SSE2
#ifdef __GNUC__

#include <emmintrin.h>

namespace piranha
{
	template <int Size>
	class mfiav_selector
	{
			p_static_check(Size == 4, "");
		public:
			static const size_t value = 4;
			static __m128i add(__m128i a, __m128i b) {
				return _mm_add_epi32(a,b);
			}
			static __m128i sub(__m128i a, __m128i b) {
				return _mm_sub_epi32(a,b);
			}
	};

	template <>
	class mfiav_selector<8>
	{
		public:
			static const size_t value = 2;
			static __m128i add(__m128i a, __m128i b) {
				return _mm_add_epi64(a,b);
			}
			static __m128i sub(__m128i a, __m128i b) {
				return _mm_sub_epi64(a,b);
			}
	};

	class max_fast_int_atom_vector
	{
			typedef mfiav_selector<sizeof(max_fast_int)> selector;
		public:
			static const size_t size = selector::value;
			void add(const max_fast_int_atom_vector &v1, const max_fast_int_atom_vector &v2) {
				m_container.m = selector::add(v1.m_container.m,v2.m_container.m);
			}
			void sub(const max_fast_int_atom_vector &v1, const max_fast_int_atom_vector &v2) {
				m_container.m = selector::sub(v1.m_container.m,v2.m_container.m);
			}
			const max_fast_int &operator[](const size_t &n) const {
				return m_container.v[n];
			}
			max_fast_int &operator[](const size_t &n) {
				return m_container.v[n];
			}
		private:
			union container_type
			{
				max_fast_int	v[size];
				__m128i			m;
			};
			__attribute__ ((aligned (16))) container_type m_container;
	};

	class double_atom_vector
	{
		public:
			static const size_t size = 2;
			void mul(const double_atom_vector &v1, const double_atom_vector &v2) {
				m_container.m = _mm_mul_pd(v1.m_container.m,v2.m_container.m);
			}
			void div(const double_atom_vector &v1, const double_atom_vector &v2) {
				m_container.m = _mm_div_pd(v1.m_container.m,v2.m_container.m);
			}
			const double &operator[](const size_t &n) const {
				return m_container.v[n];
			}
			double &operator[](const size_t &n) {
				return m_container.v[n];
			}
		private:
			union container_type
			{
				double	v[size];
				__m128d	m;
			};
			__attribute__ ((aligned (16))) container_type 				m_container;
			__attribute__ ((aligned (16))) static const container_type	m_v2;
	};

	extern __attribute__ ((aligned (16))) const double_atom_vector::container_type double_atom_vector::m_v2 = { { 2, 2 } };
}

#else
#error No SSE2 support available for this compiler
#endif // __GNUC__
#endif // _PIRANHA_SSE2

#endif // PIRANHA_VECTORIZATION_H
