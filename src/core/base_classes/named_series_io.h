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
#include <cstddef>
#include <iostream>
#include <stdexcept> // For runtime error exception.
#include <string>

#include "../config.h"
#include "../exceptions.h"
#include "../Psym.h"
#include "named_series_def.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast       static_cast<Derived *>(this)

namespace piranha
{
	// TODO: move out in static method? Where is this used?
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline std::complex<Derived> NamedSeries<__PIRANHA_NAMED_SERIES_TP>::complex() const
	{
		return std::complex<Derived>(*derived_const_cast);
	}


	// TMP for series printing.
	template <class ArgsDescr>
	inline void named_series_print_plain(std::ostream &stream,
		 const typename Ntuple<VectorPsym, boost::tuples::length<ArgsDescr>::value>::type &argsTuple)
	{
		for (std::size_t i = 0; i < argsTuple.get_head().size(); ++i) 
        {
			stream << "[" << ArgsDescr::head_type::name << "_arg]" << '\n';
			argsTuple.get_head()[i].print(stream);
		}

		named_series_print_plain<typename ArgsDescr::tail_type>(stream, argsTuple.get_tail());
	}


	template <>
	inline void named_series_print_plain<boost::tuples::null_type>(std::ostream &,
			const Ntuple<VectorPsym, boost::tuples::length<boost::tuples::null_type>::value>::type &)
	{}


	/// Print series to stream in plain format.
	/**
	 * This is the same text format used when saving series to file.
	 */
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<__PIRANHA_NAMED_SERIES_TP>::print_plain(std::ostream &stream) const
	{
		named_series_print_plain<arguments_description>(stream, argumentsTuple);
		stream << "[terms]" << std::endl;
		derived_const_cast->print_terms_plain(stream, argumentsTuple);
	}


	/// Print series to stream in pretty format.
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<__PIRANHA_NAMED_SERIES_TP>::print_pretty(std::ostream &stream) const
	{
		derived_const_cast->print_terms_pretty(stream, argumentsTuple);
	}


	/// Print in tex format.
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<__PIRANHA_NAMED_SERIES_TP>::print_tex(std::ostream &stream) const
	{
		derived_const_cast->print_terms_tex(stream, argumentsTuple);
	}


	/// Print series to stream.
	/**
	 * Equivalent to NamedSeries::print_pretty.
	 */
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<__PIRANHA_NAMED_SERIES_TP>::print(std::ostream &outStream) const
	{
		print_pretty(outStream);
	}


	/// Construct from file.
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<__PIRANHA_NAMED_SERIES_TP>::construct_from_file(const std::string &fn)
	{
		std::ifstream inf;
		inf.open(fn.c_str(), std::ios::in | std::ios::binary);
		if (inf.fail()) 
        {
			piranha_throw(std::runtime_error,"unable to open file");
		}

		// Clear the stack of unknown data.
		unknown_data.clear();
		read_sections(inf);
		trim();
	}


	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<__PIRANHA_NAMED_SERIES_TP>::read_sections(std::ifstream &inf)
	{
		std::string temp;
		while (utils::get_valid_string(inf, temp)) 
        {
			if (temp.size() > 2 && temp[0] == '[' && temp[temp.size() - 1] == ']') 
            {
				std::string sectionName = temp;
				boost::trim_if(sectionName, boost::is_any_of("[]"));
				std::cout << "New section found: " << sectionName << std::endl;
				std::vector<std::string> split_v;
				boost::split(split_v, sectionName, boost::is_any_of("_"));

				if (split_v.size() == 2 && split_v[1] == "arg") 
                {
					read_arg(inf, split_v[0]);

				} else if (sectionName == "terms") 
                {
					read_terms(inf);
					// If we found the data, then we don't want any more sections.
					return;

				} else 
                {
					std::cout << "Found unknown section '" << sectionName << "', ignoring." << std::endl;
					unknown_data.push_back(temp);
				}
			} else 
            {
				std::cout << "Found string not belonging to any (known) section: " << temp << std::endl;
				unknown_data.push_back(temp);
			}
		}
	}


	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<__PIRANHA_NAMED_SERIES_TP>::read_arg(std::ifstream &inf, const std::string &name)
	{
		// Temporary attributes for the argument.
		std::string temp_name, temp_time_eval;
		// Record stream position, so we can rewind once we finish parsing the argument.
		std::streampos cur_pos = inf.tellg();
		// Temporary storage.
		std::string temp;
		while (utils::get_valid_string(inf, temp)) 
        {
			// If we found a new section, step back the cursor before exiting.
			if (temp.size() > 2 && temp[0] == '[' && temp[temp.size() - 1] == ']') 
            {
				std::cout << "Finished parsing " << name << " argument." << std::endl;
				inf.seekg(cur_pos);
				append_arg(name, Psym(temp_name, temp_time_eval));
				return;
			}

			std::vector<std::string> split_v;
			boost::split(split_v, temp, boost::is_any_of("="));
			if (split_v.size() != 2) 
            {
				std::cout << "Invalid line in " << name << " argument section: \"" << temp << "\"" << std::endl;

			} else if (split_v[0] == "name") 
            {
				std::cout << "name=" << split_v[1] << std::endl;
				temp_name = split_v[1];

			} else if (split_v[0] == "time_eval") 
            {
				std::cout << "time_eval=" << split_v[1] << std::endl;
				temp_time_eval = split_v[1];

			} else 
            {
				std::cout << "Unknown field in " << name << " argument section: \"" << split_v[0] << "\"" << std::endl;
			}

			cur_pos = inf.tellg();
		}
	}


	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<__PIRANHA_NAMED_SERIES_TP>::read_terms(std::ifstream &inf)
	{
		typedef typename Derived::term_type term_type;
		std::string temp;
		while (!inf.eof())
        {
			getline(inf, temp, derived_const_cast->separator);
			boost::trim(temp);
			// Ignore empty lines.
			if (temp.empty()) 
            {
				continue;
			}

			try {
				term_type term(temp, argumentsTuple);
				if (!term.m_cf.is_insertable(argumentsTuple) || !term.m_key.is_insertable(argumentsTuple)) 
                {
					piranha_throw(value_error,"term not insertable in series");
				}

				derived_cast->insert(term, derived_const_cast->argumentsTuple);

			} catch (value_error &ve) 
            {
				std::cout << ve.what() << std::endl;
			}
		}
	}


	/// Save series to file.
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<__PIRANHA_NAMED_SERIES_TP>::save_to(const std::string &filename) const
	{
		std::ofstream outf(filename.c_str(), std::ios_base::trunc);
		if (outf.fail()) 
        {
			std::cout << "Error saving to file " << filename << "." << std::endl;
			outf.close();
			return;
		}

		print_plain(outf);
		outf.close();
	}


	/// Constructor from Psym and from position in the arguments set.
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	template <int N>
	inline void NamedSeries<__PIRANHA_NAMED_SERIES_TP>::construct_from_psym(const Psym &p)
	{
		piranha_assert(derived_const_cast->empty());
		append_arg<N>(p);
		derived_cast->base_construct_from_psym(p, N, argumentsTuple);
	}


	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline const typename NamedSeries<__PIRANHA_NAMED_SERIES_TP>::ArgsTupleType &
	NamedSeries<__PIRANHA_NAMED_SERIES_TP>::arguments() const
	{
		return argumentsTuple;
	}


	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<__PIRANHA_NAMED_SERIES_TP>::set_arguments(const ArgsTupleType &argsTuple)
	{
		typedef typename Derived::const_iterator const_iterator;

		const const_iterator it_f = derived_const_cast->end();
		for (const_iterator it = derived_const_cast->begin(); it != it_f; ++it) 
        {
			if (!it->m_cf.is_insertable(argsTuple) || !it->m_key.is_insertable(argsTuple) ||
				it->m_cf.needs_padding(argsTuple) || it->m_key.needs_padding(argsTuple))
			{
				piranha_throw(value_error,"incompatible arguments tuple in set_arguments()");
			}
		}

		argumentsTuple = argsTuple;
	}


	/// Construct series from a key.
	/**
	 * @see piranha::BaseSeries::base_series_from_key.
	 */
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	template <class Key>
	inline Derived NamedSeries<__PIRANHA_NAMED_SERIES_TP>::series_from_key(const Key &key) const
	{
		Derived retval(derived_const_cast->base_series_from_key(key, m_arguments));
		retval.m_arguments = m_arguments;
		retval.trim();
		return retval;
	}


	/// Construct series from a cf.
	/**
	 * @see piranha::BaseSeries::base_series_from_cf.
	 */
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	template <class Cf>
	inline Derived NamedSeries<__PIRANHA_NAMED_SERIES_TP>::series_from_cf(const Cf &cf) const
	{
		Derived retval(derived_const_cast->base_series_from_cf(cf, m_arguments));
		retval.m_arguments = m_arguments;
		retval.trim();
		return retval;
	}


	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline typename NamedSeries<__PIRANHA_NAMED_SERIES_TP>::s_iterator NamedSeries<__PIRANHA_NAMED_SERIES_TP>::s_begin() const
	{
		return s_iterator(derived_const_cast->begin(), s_iterator_generator(*derived_const_cast));
	}


	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline typename NamedSeries<__PIRANHA_NAMED_SERIES_TP>::s_iterator NamedSeries<__PIRANHA_NAMED_SERIES_TP>::s_end() const
	{
		return s_iterator(derived_const_cast->end(), s_iterator_generator(*derived_const_cast));
	}


	// Trivial destructor. It's here only to enforce a static check that cannot stay in the class definition.
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline NamedSeries<__PIRANHA_NAMED_SERIES_TP>::~NamedSeries()
	{
		p_static_check(boost::tuples::length<arguments_description>::value == Derived::echelon_level + 1, "");
	}


	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline std::ostream &operator<<(std::ostream &os, const NamedSeries<__PIRANHA_NAMED_SERIES_TP> &series)
	{
		series.print(os);
		return os;
	}
}

#undef derived_const_cast
#undef derived_cast

#endif
