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

#ifndef PIRANHA_NAMED_SERIES_IO_H
#define PIRANHA_NAMED_SERIES_IO_H

#include <boost/algorithm/string.hpp>
#include <boost/tuple/tuple.hpp>
#include <complex>
#include <iostream>
#include <sstream>
#include <string>

#include "../shared_args.h"
#include "../config.h"
#include "../exceptions.h"
#include "../psym.h"

namespace piranha
{
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline std::complex<Derived> named_series<__PIRANHA_NAMED_SERIES_TP>::complex() const
	{
		return std::complex<Derived>(*derived_const_cast);
	}

	// TMP for series printing.
	template <class ArgsDescr>
	inline void named_series_print_plain(std::ostream &stream,
										 const typename ntuple<vector_psym_p, boost::tuples::length<ArgsDescr>::value>::type &args_tuple)
	{
		for (size_t i = 0; i < args_tuple.get_head().size(); ++i) {
			stream << "[" << ArgsDescr::head_type::name << "_arg]" << '\n';
			args_tuple.get_head()[i]->print(stream);
		}
		named_series_print_plain<typename ArgsDescr::tail_type>(stream, args_tuple.get_tail());
	}

	template <>
	inline void named_series_print_plain<boost::tuples::null_type>(std::ostream &,
			const ntuple<vector_psym_p, boost::tuples::length<boost::tuples::null_type>::value>::type &)
	{}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline void named_series<__PIRANHA_NAMED_SERIES_TP>::print_plain(std::ostream &stream, int limit) const
	{
		named_series_print_plain<arguments_description>(stream, m_arguments);
		stream << "[terms]" << std::endl;
		derived_const_cast->print_terms_plain(stream, m_arguments, limit);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline void named_series<__PIRANHA_NAMED_SERIES_TP>::print_latex(std::ostream &, int) const
	{
		// TODO: implement.
	}

	/// Print series to stream.
	/**
	 * Print first "limit" terms. If limit is negative, print all terms. The output format is read
	 * from the piranha::stream_manager class.
	 */
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline void named_series<__PIRANHA_NAMED_SERIES_TP>::print(std::ostream &out_stream, int limit) const
	{
		switch (stream_manager::format()) {
		case stream_manager::plain:
			print_plain(out_stream, limit);
			break;
		case stream_manager::latex:
			print_latex(out_stream, limit);
		}
	}

	/// Construct from file.
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline void named_series<__PIRANHA_NAMED_SERIES_TP>::construct_from_file(const std::string &fn)
	{
		std::ifstream inf;
		std::string filename = utils::open_file(fn, inf);
		// Read from file
		if (inf.is_open()) {
			// Clear the stack of unknown data.
			unknown_data.clear();
			read_sections(inf);
			std::cout << "EOF" << std::endl;
			// Close file
			inf.close();
		}
		trim();
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline void named_series<__PIRANHA_NAMED_SERIES_TP>::read_sections(std::ifstream &inf)
	{
		std::string temp;
		while (utils::get_valid_string(inf, temp)) {
			if (temp.size() > 2 && temp[0] == '[' && temp[temp.size()-1] == ']') {
				std::string sec_name = temp;
				boost::trim_if(sec_name, boost::is_any_of("[]"));
				std::cout << "New section found: " << sec_name << std::endl;
				std::vector<std::string> split_v;
				boost::split(split_v, sec_name, boost::is_any_of("_"));
				if (split_v.size() == 2 && split_v[1] == "arg") {
					read_arg(inf, split_v[0]);
				} else if (sec_name == "terms") {
					read_terms(inf);
					// If we found the data, then we don't want any more sections.
					return;
				} else {
					std::cout << "Found unknown section '" << sec_name << "', ignoring." << std::endl;
					unknown_data.push_back(temp);
				}
			} else {
				std::cout << "Found string not belonging to any (known) section: " << temp << std::endl;
				unknown_data.push_back(temp);
			}
		}
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline void named_series<__PIRANHA_NAMED_SERIES_TP>::read_arg(std::ifstream &inf, const std::string &name)
	{
		// Temporary attributes for the argument.
		std::string temp_name, temp_time_eval;
		// Record stream position, so we can rewind once we finish parsing the argument.
		std::streampos cur_pos = inf.tellg();
		// Temporary storage.
		std::string temp;
		while (utils::get_valid_string(inf, temp)) {
			// If we found a new section, step back the cursor before exiting.
			if (temp.size() > 2 && temp[0] == '[' && temp[temp.size()-1] == ']') {
				std::cout << "Finished parsing " << name << " argument." << std::endl;
				inf.seekg(cur_pos);
				append_arg(name, psyms::get_pointer(psym(temp_name, temp_time_eval)));
				return;
			}
			std::vector<std::string> split_v;
			boost::split(split_v, temp, boost::is_any_of("="));
			if (split_v.size() != 2) {
				std::cout << "Invalid line in " << name << " argument section: \"" << temp << "\"" << std::endl;
			} else if (split_v[0] == "name") {
				std::cout << "name=" << split_v[1] << std::endl;
				temp_name = split_v[1];
			} else if (split_v[0] == "time_eval") {
				std::cout << "time_eval=" << split_v[1] << std::endl;
				temp_time_eval = split_v[1];
			} else {
				std::cout << "Unknown field in " << name << " argument section: \"" << split_v[0] << "\"" << std::endl;
			}
			cur_pos = inf.tellg();
		}
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline void named_series<__PIRANHA_NAMED_SERIES_TP>::read_terms(std::ifstream &inf)
	{
		typedef typename Derived::term_type term_type;
		std::string temp;
		while (!inf.eof()) {
			getline(inf, temp, derived_const_cast->separator);
			boost::trim(temp);
			// Ignore empty lines.
			if (temp.empty()) {
				continue;
			}
			try {
				term_type term(temp, m_arguments);
				if (!term.is_insertable(m_arguments)) {
					throw bad_input("Term not insertable in named series.");
				}
				derived_cast->insert(term, derived_const_cast->m_arguments);
			} catch (bad_input &b) {
				std::cout << b.what() << std::endl;
			}
		}
	}

	/// Save series to file.
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline void named_series<__PIRANHA_NAMED_SERIES_TP>::save_to(const std::string &filename) const
	{
		std::ofstream outf(filename.c_str(), std::ios_base::trunc);
		if (outf.fail()) {
			std::cout << "Error saving to file " << filename << "." << std::endl;
			outf.close();
			return;
		}
		print(outf);
		outf.close();
	}

	/// Print the series, including header, to screen.
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline std::string named_series<__PIRANHA_NAMED_SERIES_TP>::py_repr() const
	{
		std::ostringstream stream;
		print(stream, -1);
		std::string retval(stream.str());
		return retval;
	}

	/// Constructor from psym and from position in the arguments set.
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	template <int N>
	inline void named_series<__PIRANHA_NAMED_SERIES_TP>::construct_from_psym(const psym &psym)
	{
		p_assert(derived_const_cast->empty());
		psym_p p(psyms::get_pointer(psym));
		append_arg<N>(p);
		derived_cast->construct_from_psym_p(p, N, m_arguments);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline const typename named_series<__PIRANHA_NAMED_SERIES_TP>::args_tuple_type &
	named_series<__PIRANHA_NAMED_SERIES_TP>::arguments() const
	{
		return m_arguments;
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline typename named_series<__PIRANHA_NAMED_SERIES_TP>::args_tuple_type
	named_series<__PIRANHA_NAMED_SERIES_TP>::py_arguments() const
	{
		return m_arguments;
	}

	template <class ArgsDescr>
	class arguments_type_report_helper
	{
		public:
			template <class ArgsTuple>
			static void run(const ArgsTuple &args_tuple, std::string &report) {
				report += ArgsDescr::head_type::name;
				report += "\n";
				arguments_type_report_helper<typename ArgsDescr::tail_type>::run(args_tuple.get_tail(),report);
			}
	};

	template <>
	class arguments_type_report_helper<boost::tuples::null_type>
	{
		public:
			template <class ArgsTuple>
			static void run(const ArgsTuple &, const std::string &) {}
	};

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline std::string named_series<__PIRANHA_NAMED_SERIES_TP>::py_arguments_description() const
	{
		std::string retval;
		arguments_type_report_helper<arguments_description>::run(m_arguments,retval);
		return retval;
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline void named_series<__PIRANHA_NAMED_SERIES_TP>::py_shared_arguments_set() const
	{
		shared_args::set(m_arguments);
	}
}

#endif
