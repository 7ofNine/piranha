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

#include <boost/lexical_cast.hpp>
#include <boost/python/class.hpp>
#include <boost/tuple/tuple.hpp>
#include <sstream>
#include <string>

namespace pyranha
{
	template <class T>
	inline T py_copy(const T &x) {
		return T(x);
	}

	template <class T>
	std::string py_print_to_string(const T &origin)
	{
		std::ostringstream stream;
		origin.print(stream);
		return stream.str();
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

	template <class ArgsDescr>
	class arguments_type_report_helper
	{
		public:
			template <class ArgsTuple>
			static void run(const ArgsTuple &args_tuple, std::string &report) {
				report += ArgsDescr::head_type::name;
				report += ":";
				report += boost::lexical_cast<std::string>(args_tuple.get_head().size());
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

	template <class NamedSeries>
	inline size_t psi0(const NamedSeries &s) {
		return s.psi(0, 1);
	}

	template <class NamedSeries>
	inline size_t psi1(const NamedSeries &s, const int &n) {
		return s.psi(n, 1);
	}

	template <class NamedSeries>
	inline size_t psi2(const NamedSeries &s, const int &n, const int &step) {
		return s.psi(n, step);
	}
}

#endif
