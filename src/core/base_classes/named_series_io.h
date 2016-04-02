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

#include <complex>
#include <cstddef>
#include <iostream>
#include <stdexcept> // For runtime error exception.
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/tuple/tuple.hpp>

#include "../config.h"
#include "../exceptions.h"
#include "../Psym.h"
#include "named_series_def.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast       static_cast<Derived *>(this)

namespace piranha
{
	// TODO: move out in static method? Where is this used?
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline std::complex<Derived> NamedSeries<PIRANHA_NAMED_SERIES_TP>::complex() const
	{
		return std::complex<Derived>(*derived_const_cast);
	}


	// TMP for series printing.
	template <class ArgsDescr>
	inline void namedSeriesPrintPlain(std::ostream &stream, typename NTuple<VectorPsym, boost::tuples::length<ArgsDescr>::value>::Type const &argsTuple)
	{
		for (std::size_t i = 0; i < argsTuple.get_head().size(); ++i) 
        {
			stream << "[" << ArgsDescr::head_type::name << "_arg]" << '\n';
			argsTuple.get_head()[i].print(stream);
		}

		namedSeriesPrintPlain<typename ArgsDescr::tail_type>(stream, argsTuple.get_tail());
	}


	template <>
	inline void namedSeriesPrintPlain<boost::tuples::null_type>(std::ostream &, NTuple<VectorPsym, boost::tuples::length<boost::tuples::null_type>::value>::Type const &)
	{}


	/// Print series to stream in plain format.
	/**
	 * This is the same text format used when saving series to file.
	 */
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<PIRANHA_NAMED_SERIES_TP>::printPlain(std::ostream &stream) const
	{
		namedSeriesPrintPlain<ArgumentsDescription>(stream, argumentsTuple);
		stream << "[terms]" << std::endl;
		derived_const_cast->printTermsPlain(stream, argumentsTuple);
	}


	/// Print series to stream in pretty format.
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<PIRANHA_NAMED_SERIES_TP>::printPretty(std::ostream &stream) const
	{
		derived_const_cast->printTermsPretty(stream, argumentsTuple);
	}


	/// Print in tex format.
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<PIRANHA_NAMED_SERIES_TP>::printTex(std::ostream &stream) const
	{
		derived_const_cast->printTermsTEX(stream, argumentsTuple);
	}


	/// Print series to stream.
	/**
	 * Equivalent to NamedSeries::print_pretty.
	 */
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<PIRANHA_NAMED_SERIES_TP>::print(std::ostream &outStream) const
	{
		printPretty(outStream);
	}


	/// Construct from file.
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<PIRANHA_NAMED_SERIES_TP>::constructFromFile(std::string const &fileName)
	{
		std::ifstream infile;
		infile.open(fileName.c_str(), std::ios::in | std::ios::binary);
		if (infile.fail()) 
        {
			PIRANHA_THROW(std::runtime_error, "Unable to open file " + fileName);
		}

		// Clear the stack of unknown data.
		unknownData.clear();
		readSections(infile);
		trim();
	}


	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<PIRANHA_NAMED_SERIES_TP>::readSections(std::ifstream &infile)
	{
		std::string input;
		while (utils::get_valid_string(infile, input)) 
        {
			if (input.size() > 2 && input[0] == '[' && input[input.size() - 1] == ']') 
            {
				std::string sectionName = input;
				boost::trim_if(sectionName, boost::is_any_of("[]"));

				std::cout << "New section found: " << sectionName << std::endl;
				
                std::vector<std::string> splitInput;
				boost::split(splitInput, sectionName, boost::is_any_of("_"));

				if (splitInput.size() == 2 && splitInput[1] == "arg") 
                {
					readArg(infile, splitInput[0]);

				} else if (sectionName == "terms") 
                {
					readTerms(infile);
					// If we found the data, then we don't want any more sections.
					return;

				} else 
                {
					std::cout << "Found unknown section '" << sectionName << "', ignoring." << std::endl;
					unknownData.push_back(input);
				}
			} else 
            {
				std::cout << "Found string not belonging to any (known) section: " << input << std::endl;
				unknownData.push_back(input);
			}
		}
	}


	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<PIRANHA_NAMED_SERIES_TP>::readArg(std::ifstream &infile, std::string const &argumentTypeName)
	{
		// Temporary attributes for the argument.
		std::string name;
        std::string timeEval;
		// Record stream position, so we can rewind once we finish parsing the argument.
		std::streampos currentPosition = infile.tellg();
		// Temporary storage.
		std::string line;
		while (utils::get_valid_string(infile, line)) 
        {
			// If we found a new section, step back the cursor before exiting.
			if (line.size() > 2 && line[0] == '[' && line[line.size() - 1] == ']') 
            {
				std::cout << "Finished parsing " << argumentTypeName << " argument." << std::endl;
				infile.seekg(currentPosition);
				appendArg(line, Psym(name, timeEval));
				return;
			}

			std::vector<std::string> splitArgument;
			boost::split(splitArgument, line, boost::is_any_of("="));
			if (splitArgument.size() != 2) 
            {
				std::cout << "Invalid line in " << argumentTypeName << " argument section: \"" << line << "\"" << std::endl;

			} else if (splitArgument[0] == "name") 
            {
				std::cout << "name = " << splitArgument[1] << std::endl;
				name = splitArgument[1];

			} else if (splitArgument[0] == "time_eval") 
            {
				std::cout << "time_eval = " << splitArgument[1] << std::endl;
				timeEval = splitArgument[1];

			} else 
            {
				std::cout << "Unknown field in " << line << " argument section: \"" << splitArgument[0] << "\"" << std::endl;
			}

			currentPosition = infile.tellg();
		}
	}


	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<PIRANHA_NAMED_SERIES_TP>::readTerms(std::ifstream &infile)
	{
		typedef typename Derived::TermType TermType;
		std::string line;

		while (!infile.eof())
        {
			getline(infile, line, derived_const_cast->separator);
			boost::trim(line);
			// Ignore empty lines.
			if (line.empty()) 
            {
				continue;
			}

			try {
				TermType term(line, argumentsTuple);
				if (!term.cf.isInsertable(argumentsTuple) || !term.key.isInsertable(argumentsTuple)) 
                {
					PIRANHA_THROW(value_error, "Term not insertable in series");
				}

				derived_cast->insert(term, derived_const_cast->argumentsTuple);

			} catch (value_error &ve) 
            {
				std::cout << ve.what() << std::endl;
			}
		}
	}


	/// Save series to file.
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<PIRANHA_NAMED_SERIES_TP>::saveTo(std::string const &filename) const
	{
		std::ofstream outfile(filename.c_str(), std::ios_base::trunc);
		if (outfile.fail()) 
        {
			std::cout << "Error saving to file " << filename << "." << std::endl;
			outfile.close();
			return;
		}

		printPlain(outfile);
		outfile.close();
	}


	/// Constructor from Psym and from position in the arguments set.
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <int N>
	inline void NamedSeries<PIRANHA_NAMED_SERIES_TP>::constructFromPsym(Psym const &p)
	{
		PIRANHA_ASSERT(derived_const_cast->empty());
		appendArg<N>(p);
		derived_cast->baseConstructFromPsym(p, N, argumentsTuple);
	}


	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline const typename NamedSeries<PIRANHA_NAMED_SERIES_TP>::ArgsTupleType &
	NamedSeries<PIRANHA_NAMED_SERIES_TP>::arguments() const
	{
		return argumentsTuple;
	}


	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline void NamedSeries<PIRANHA_NAMED_SERIES_TP>::setArguments(const ArgsTupleType &argsTuple)
	{
		typedef typename Derived::const_iterator Iterator;

		const Iterator itEnd= derived_const_cast->end();
		for (Iterator it = derived_const_cast->begin(); it != itEnd; ++it) 
        {
			if (!it->cf.isInsertable(argsTuple) || !it->key.isInsertable(argsTuple) ||
			     it->cf.needsPadding(argsTuple) ||  it->key.needsPadding(argsTuple))
			{
				PIRANHA_THROW(value_error, "Incompatible arguments tuple in set_arguments()");
			}
		}

		argumentsTuple = argsTuple;
	}


	/// Construct series from a key.
	/**
	 * @see piranha::BaseSeries::baseSeriesFromKey.
	 */
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <class Key>
	inline Derived NamedSeries<PIRANHA_NAMED_SERIES_TP>::seriesFromKey(Key const &key) const
	{
		Derived retval(derived_const_cast->baseSeriesFromKey(key, argumentsTuple));
		retval.argumentsTuple = argumentsTuple;
		retval.trim();

		return retval;
	}


	/// Construct series from a cf.
	/**
	 * @see piranha::BaseSeries::baseSeriesFromCf.
	 */
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	template <class Cf>
	inline Derived NamedSeries<PIRANHA_NAMED_SERIES_TP>::seriesFromCf(Cf const &cf) const
	{
		Derived retval(derived_const_cast->baseSeriesFromCf(cf, argumentsTuple));
		retval.argumentsTuple = argumentsTuple;
		retval.trim();

		return retval;
	}


	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline typename NamedSeries<PIRANHA_NAMED_SERIES_TP>::SeriesIterator NamedSeries<PIRANHA_NAMED_SERIES_TP>::itBegin() const
	{
		return SeriesIterator(derived_const_cast->begin(), SeriesIteratorGenerator(*derived_const_cast));
	}


	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline typename NamedSeries<PIRANHA_NAMED_SERIES_TP>::SeriesIterator NamedSeries<PIRANHA_NAMED_SERIES_TP>::itEnd() const
	{
		return SeriesIterator(derived_const_cast->end(), SeriesIteratorGenerator(*derived_const_cast));
	}


	// Trivial destructor. It's here only to enforce a static check that cannot stay in the class definition.
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline NamedSeries<PIRANHA_NAMED_SERIES_TP>::~NamedSeries()
	{
		PIRANHA_STATIC_CHECK(boost::tuples::length<ArgumentsDescription>::value == Derived::echelonLevel + 1, "");
	}


	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline std::ostream &operator<<(std::ostream &os, NamedSeries<PIRANHA_NAMED_SERIES_TP> const &series)
	{
		series.print(os);
		return os;
	}
}

#undef derived_const_cast
#undef derived_cast

#endif
