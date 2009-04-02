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

#ifndef PIRANHA_COMMON_COMPARISONS_H
#define PIRANHA_COMMON_COMPARISONS_H

namespace piranha
{
	template <class ArgsTuple>
	class cf_norm_comparison
	{
		public:
			cf_norm_comparison(const ArgsTuple &args_tuple):m_args_tuple(args_tuple) {}
			template <class Term>
			bool operator()(const Term &t1, const Term &t2) const {
				return t1.m_cf.norm_(m_args_tuple) > t2.m_cf.norm_(m_args_tuple);
			}
		private:
			const ArgsTuple &m_args_tuple;
	};

	template <class ArgsTuple>
	class cf_norm_comparison_reverse
	{
		public:
			cf_norm_comparison_reverse(const ArgsTuple &args_tuple):m_args_tuple(args_tuple) {}
			template <class Term>
			bool operator()(const Term *t1, const Term *t2) const {
				return t1->m_cf.norm_(m_args_tuple) < t2->m_cf.norm_(m_args_tuple);
			}
		private:
			const ArgsTuple &m_args_tuple;
	};

	template <class ArgsTuple>
	class term_cf_min_degree_comparison
	{
		public:
			term_cf_min_degree_comparison(const ArgsTuple &args_tuple):m_args_tuple(args_tuple) {}
			template <class Term>
			bool operator()(const Term &t1, const Term &t2) const {
				const int d1 = t1.m_cf.min_degree(), d2 = t2.m_cf.min_degree();
				if (d1 == d2) {
					return t1.m_cf.norm_(m_args_tuple) > t2.m_cf.norm_(m_args_tuple);
				} else {
					return d1 < d2;
				}
			}
		private:
			const ArgsTuple &m_args_tuple;
	};

	template <class ArgsTuple>
	class term_key_degree_comparison
	{
		public:
			term_key_degree_comparison(const ArgsTuple &) {}
			template <class Term>
			bool operator()(const Term *t1, const Term *t2) const {
				return t1->m_key.degree() < t2->m_key.degree();
			}
	};
}

#endif
