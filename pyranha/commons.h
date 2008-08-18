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

#ifndef PYRANHA_COMMONS_H
#define PYRANHA_COMMONS_H

#include <boost/python/class.hpp>
#include <boost/tuple/tuple.hpp>
#include <sstream>
#include <string>

#include "../src/core/exceptions.h"
#include "../src/core/shared_args.h"

namespace pyranha
{
	template <class T>
	inline T py_copy(const T &x) {
		return T(x);
	}

	template <class Vector>
	inline typename Vector::value_type py_vector_getitem(const Vector &v, const piranha::max_fast_int &n_)
	{
		int n = n_;
		const size_t size = v.size();
		if (n_ < 0) {
			n = n_ + size;
		}
		if (n < 0 || static_cast<size_t>(n) >= size) {
			std::ostringstream stream;
			stream << "Index " << n << " is out of range.";
			throw piranha::unsuitable(stream.str());
		}
		return v[n];
	}

	template <class T>
	std::string py_print_to_string(const T &origin)
	{
		std::ostringstream stream;
		origin.print(stream);
		return stream.str();
	}

	template <class Series>
	inline typename Series::const_iterator py_series_begin(const Series &s)
	{
		return s.begin();
	}

	template <class Series>
	inline typename Series::const_iterator py_series_end(const Series &s)
	{
		return s.end();
	}

	template <class Series, class Term>
	inline void py_series_append(Series &s, const Term &t) {
		s.insert(t,piranha::shared_args::get());
	}

	template <class Term>
	inline boost::python::class_<Term> py_series_term(const std::string &name, const std::string &description)
	{
		boost::python::class_<Term> term_inst((name + "_term").c_str(),
				(std::string("Term for: ") + description).c_str());
		term_inst.def_readwrite("cf", &Term::m_cf);
		term_inst.def_readwrite("key", &Term::m_key);
		return term_inst;
	}

	template <class Series>
	inline Series py_series_get_real(const std::complex<Series> &c)
	{
		return c.real();
	}

	template <class Series>
	inline void py_series_set_real(std::complex<Series> &c, const Series &r)
	{
		c.real(r);
	}

	template <class Series>
	inline Series py_series_get_imag(const std::complex<Series> &c)
	{
		return c.imag();
	}

	template <class Series>
	inline void py_series_set_imag(std::complex<Series> &c, const Series &i)
	{
		c.imag(i);
	}

	template <class ArgsDescr>
	class arguments_type_report_helper
	{
		public:
			template <class ArgsTuple>
			static void run(const ArgsTuple &args_tuple, std::string &report) {
				report += ArgsDescr::head_type::name;
				report += "\n";
				arguments_type_report_helper<typename ArgsDescr::tail_type>::run(args_tuple.get_tail(), report);
			}
	};

	template <>
	class arguments_type_report_helper<boost::tuples::null_type>
	{
		public:
			template <class ArgsTuple>
			static void run(const ArgsTuple &, const std::string &) {}
	};

	template <class Series>
	inline std::string py_series_arguments_description(const Series &s)
	{
		std::string retval;
		arguments_type_report_helper<typename Series::arguments_description>::run(s.arguments(), retval);
		return retval;
	}

	template <class Series>
	inline typename Series::args_tuple_type py_series_arguments(const Series &s)
	{
		return s.arguments();
	}

	template <class Series>
	inline void py_series_set_shared_arguments(const Series &s)
	{
		piranha::shared_args::set(s.arguments());
	}
}

#endif
