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

#ifndef PIRANHA_TYPE_TRAITS_H
#define PIRANHA_TYPE_TRAITS_H

#include <boost/type_traits/integral_constant.hpp>
#include <complex>

namespace piranha
{
	// Fwd declaration.
	template <class Cf, class Key> class tetd_helper;

	template <class Term>
	class term_eval_type_determiner
	{
		public:
			typedef typename tetd_helper<typename Term::cf_type::eval_type, typename Term::key_type::eval_type>::type type;
	};

	template <class CfEval, class KeyEval>
	class tetd_helper
	{
		public:
			typedef double type;
	};

	template <>
	class tetd_helper<double,std::complex<double> >
	{
		public:
			typedef std::complex<double> type;
	};

	template <>
	class tetd_helper<std::complex<double>,double>
	{
		public:
			typedef std::complex<double> type;
	};

	template <>
	class tetd_helper<std::complex<double>,std::complex<double> >
	{
		public:
			typedef std::complex<double> type;
	};

	template <class>
	struct is_lightweight: boost::false_type {};
}
#endif
