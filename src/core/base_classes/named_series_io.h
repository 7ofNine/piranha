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

#include <algorithm>
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

namespace {
    bool createFile(std::string const & fileName, std::ofstream & outfile)
    {
        outfile = std::ofstream(fileName.c_str(), std::ios_base::trunc);
        if (outfile.fail())
        {
            std::cout << "Error printing series to file " << fileName << "." << std::endl;
            outfile.close();
            return false;
        }else
        {
            return true;
        }

    }
}

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

	// print series in a sortet and structured format
	//
	// expSymbols: a vector of Psyms that are exponential variables. The index gives the sorting for printing it
	// trigSymbols: a vector of Psyms that are trigonometric variables. The index gives the sorting for printing it
	// in general the sorting is done according to trigonometric first and then according to exponential
	// this is how they are typically printed in publications but there is no canonical form one can use. It is 
	// all personal preference
	// we have to do multiple sorting. first on the trigonometric arguments and then on the exponential arguments
	// 
	// the printed file is not intended to be read back into piranha. For that purpose use saveTo.
	//

	template<PIRANHA_NAMED_SERIES_TP_DECL>
	inline void  NamedSeries<PIRANHA_NAMED_SERIES_TP>::printToSorted(std::string const & fileName, VectorPsym const & expSymbols, VectorPsym const & trigSymbols) const
	{
        std::ofstream outfile;
        if (!createFile(fileName, outfile))
        {
            return;
        }


		auto positionTuple = psyms2pos(trigSymbols, argumentsTuple); // get a tuple for where the symbols are in the terms
		// get trigonometric positions first
		auto const positions = positionTuple.get<1>();               // TODO: how to get to the trig and exponential ones separately (ie. where is the position hidden in the class)
		                                                             // we should be able to get that from the args descriptor 
		                                                             // this is a vector of pairs fist: true/false = present; second: index into vector 
		auto terms = getTrigSortedTerms(positions);                  // a vector of pointers to the terms sorted by trigonometric argument
		
		// now we are sorted by trigonometric argument
		// next we have to sort the coefficient by exponents
		// and then finaly print it
		printToSorted(outfile, terms, expSymbols, trigSymbols, positions);
		//printPlain(outfile);
		outfile.close(); 
	}

	template<PIRANHA_NAMED_SERIES_TP_DECL>
	inline void  NamedSeries<PIRANHA_NAMED_SERIES_TP>::printToSorted(std::ofstream & outfile, std::vector<Term const*> trigSortedTerms, VectorPsym const & expSymbols, 
                                                                     VectorPsym const & trigSymbols, std::vector<std::pair<bool, std::size_t> > const & trigPositions) const
	{
		typename Term::CfType check;
		typedef std::vector<Term const*>::size_type Index;
		//get positions of exponent arguments
		auto const positionTuple = psyms2pos(expSymbols, argumentsTuple);
		auto const expPositions = positionTuple.get<0>();

        // determine the max print length of the coefficient, to make it look pretty.
        // should be done more efficiently, and not by scanning the series twice
        std::size_t maxCoeffLength = 0;
        for (Index i = 0; i < trigSortedTerms.size(); ++i)
        {
            typename Term::CfType coeff = trigSortedTerms[i]->get<0>();
            auto terms = getExpoSortedCoefficient(coeff, expPositions);
            for (std::size_t j = 0; j < terms.size(); ++j)
            {
                maxCoeffLength = std::max(maxCoeffLength, (terms[j]->get<0>()).printLength());
            }
        }

        static std::string const coefficientText("Coefficient");
        static std::string const termIdText("TermId");
        std::size_t const printCoeffWidth = std::max(coefficientText.length() + 2 + 2, maxCoeffLength); // length of word 'Coefficient' and 2 leading and trailing blanks each.

        // write header of file. 
        outfile << std::setw(termIdText.length()) << termIdText << " " << std::setw(printCoeffWidth) << std::right << coefficientText << " ";
        for (std::size_t i = 0; i < expSymbols.size(); ++i)
        {
            outfile << std::setw(6) << std::right << expSymbols[i].getName() << " ";
        }

        outfile << "  |  ";

        for (std::size_t i = 0; i < trigSymbols.size(); ++i)
        {
            outfile << std::setw(6) << std::right << trigSymbols[i].getName() << " ";
        }
        outfile << std::endl;


		int termId = 0;
		for (Index i = 0; i < trigSortedTerms.size(); ++i)
		{
			typename Term::CfType coeff = trigSortedTerms[i]->get<0>();

            outfile << i + 1 <<":" << endl; // group index, makes it easier to find in the listing
			auto terms = getExpoSortedCoefficient(coeff, expPositions); // these are paointers to the single terms as they are in the coefficient split out and sorted according to exponent and position
			for (decltype(terms.size()) j = 0, e = terms.size(); j < e; ++j)
			{

				outfile << std::setw(6) << (++termId) << " ";
				outfile << std::setw(printCoeffWidth) << std::right << terms[j]->get<0>(); // the coefficient
				outfile << " ";
				//outfile << setw(5) << right;
				terms[j]->get<1>().printPlainSorted(outfile, expPositions, argumentsTuple);
                outfile << "  |  ";
				trigSortedTerms[i]->get<1>().printPlainSorted(outfile, trigPositions, argumentsTuple);
				outfile << std::endl;
			}
			outfile << std::endl;
		}
	}

	template<PIRANHA_NAMED_SERIES_TP_DECL>
	inline std::vector<typename Term::CfType::TermType const*>  NamedSeries<PIRANHA_NAMED_SERIES_TP>::getExpoSortedCoefficient(typename Term::CfType const & coeff, std::vector<std::pair<bool, size_t> > const & positions) const
	{

		class CompareExpo
		{
		public:
			explicit CompareExpo(std::vector<std::pair<bool, std::size_t> > const & positions):positions(positions) {}

			bool operator()(typename Term::CfType::TermType const * const t1, typename Term::CfType::TermType const * const t2) const
			{
				auto const expoKey1 = t1->get<1>();
				auto const expoKey2 = t2->get<1>();
				PIRANHA_ASSERT(positions.size() == expoKey1.size());
				PIRANHA_ASSERT(positions.size() == expoKey2.size());
				for (decltype(positions.size()) i = 0, e = positions.size(); i < e; ++i)
				{
					if (positions[i].first)
					{
						auto currentIndex = positions[i].second;
						
						if (expoKey1[currentIndex] < expoKey2[currentIndex])
						{
							return true;
						} else if (expoKey1[currentIndex] > expoKey2[currentIndex])
						{
							return false;
						}
					}
				}
				return false;
			}

		private:
			std::vector<std::pair<bool, std::size_t> > const & positions;
		};
		

		// first create the vector of pointers
		typedef std::vector<typename Term::CfType::TermType const*> RetValType;
		RetValType retval;

		// create vector of pointers to the series terms in the retval vector
		// the _1 means the first parameter of the functional operator (boost::lambda). The expresion determines the address (&) of the series term
		// iterated through (series.begin()..series.end() and this address is inserted (insert_iterator) into the retval vector starting at retval.begin() 
		std::transform(coeff.begin(), coeff.end(), std::insert_iterator< RetValType >(retval, retval.begin()), &boost::lambda::_1);
		std::sort(retval.begin(), retval.end(), CompareExpo(positions));

		// the ultimate result
		return retval;

	}

	//

	template <PIRANHA_NAMED_SERIES_TP_DECL>
	inline std::vector<Term const *>  NamedSeries<PIRANHA_NAMED_SERIES_TP>::getTrigSortedTerms(std::vector<std::pair<bool, std::size_t> > const &positions) const 
	{
    	// first create the vector of pointers
		typedef std::vector<Term const *> RetValType;
		RetValType retval;

		// create vector of pointers to the series terms in the retval vector
		// the _1 means the first parameter of the functional operator (boost::lambda). The expresion determines the address (&) of the series term
		// iterated through (series.begin()..series.end() and this address is inserted (insert_iterator) into the retval vector starting at retval.begin() 
		std::transform(derived_const_cast->begin(), derived_const_cast->end(), std::insert_iterator< RetValType >(retval, retval.begin()), &boost::lambda::_1);
		std::sort(retval.begin(), retval.end(), CompareTrig(positions));

		// the ultimate result
		return retval;
	}

    template <PIRANHA_NAMED_SERIES_TP_DECL>
    void NamedSeries<PIRANHA_NAMED_SERIES_TP>::printToSequenced(std::string const & fileName, VectorPsym const & expSymbols, VectorPsym const & trigSymbols, PrintSequenceType const & sequence) const
    {
        std::ofstream outfile;
        if (!createFile(fileName, outfile))
        {
            return;
        }

        auto positionTuple = psyms2pos(trigSymbols, argumentsTuple); // get a tuple for where the symbols are in the terms
                                                                     // get trigonometric positions first
        auto const positions = positionTuple.get<1>();               // TODO: how to get to the trig and exponential ones separately (ie. where is the position hidden in the class)
                                                                     // we should be able to get that from the args descriptor 
                                                                     // this is a vector of pairs fist: true/false = present; second: index into vector 
        auto terms = getTrigSequencedTerms(positions, sequence);     // a vector of pointers to the terms sorted by trigonometric argument

                                                                     // now we are sorted by trigonometric argument
                                                                     // next we have to sort the coefficient by exponents
                                                                     // and then finaly print it
        printToSorted(outfile, terms, expSymbols, trigSymbols, positions);
        //printPlain(outfile);
        outfile.close(); // bad programing outfile is hidden in createFile
    }


    // sort terms according to a pre-given sequence (mainly for comparison with printed results to avoid long manual searches)
    // the idea is to give a sequence of trigonmetric keys and return the series with the terms in that sequence and the remaining terms sorted as usual
    // This should make it easier to comapre historic printed results with a new derivation of the result. This is slow and should be rarely used
    template <PIRANHA_NAMED_SERIES_TP_DECL>
    inline std::vector<Term const *>   NamedSeries<PIRANHA_NAMED_SERIES_TP>::getTrigSequencedTerms(std::vector<std::pair<bool, std::size_t> > const &positions, PrintSequenceType const & sequence) const
    {
    
        // first: morph all the input sequence terms into the proper sequence according to position
        PrintSequenceType normalSequence;
        std::transform(sequence.begin(), sequence.end(), std::insert_iterator< PrintSequenceType >(normalSequence, normalSequence.begin()),
                       [&positions](typename Term::KeyType const & trigKey) -> typename Term::KeyType 
                       {    
                            Term::KeyType newKey;
                            PIRANHA_ASSERT(positions.size() == trigKey.size())
                            newKey.resize(trigKey.size());
                                for (decltype(positions.size()) i = 0, e = positions.size(); i < e; ++i)
                                {
                                    if (positions[i].first)
                                    {
                                        newKey[positions[i].second] = trigKey[i];
                                    }
                                }
                            newKey.setFlavour(trigKey.getFlavour());
                            return newKey;
                       });

        // first create the vector of pointers
        typedef std::vector<Term const *> RetValType;

        // we don't want to destroy the original series and we work with pointers to the terms anyways
        RetValType tempResult;
        std::transform(derived_const_cast->begin(), derived_const_cast->end(), std::insert_iterator< RetValType >(tempResult, tempResult.begin()), &boost::lambda::_1);
        
        RetValType retval;
        // now find all the elements given by the sequence. By definition a trigKey uniquely identifies an element of the series (echelon level 1)
        // if nor found put it in a temporaryresult that will get standard sorted
        for (Term::KeyType t : normalSequence)
        {
            decltype(tempResult.end()) it;
            it = std::find_if(tempResult.begin(), tempResult.end(),
                              [&t](RetValType::value_type const &st)->bool
                              {
                                  auto const s = st->get<1>();
                                  return (s == t) || (-s == t);
                              } 
                              );
            // found the sought for element. It doesn't necessarily exist 
            if (it != tempResult.end())
            {
                retval.push_back(*it);
                tempResult.erase(it);
            }
        }

        std::sort(tempResult.begin(), tempResult.end(), CompareTrig(positions));// sort the remainder i.e. those elemenst that are not in the sequence
        retval.insert(retval.end(), tempResult.begin(), tempResult.end());

        return retval;
            }


	//
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
				appendArg(argumentTypeName, Psym(name, timeEval)); // TODO: No protection against empty name ???
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
		static_assert(boost::tuples::length<ArgumentsDescription>::value == Derived::echelonLevel + 1, "");
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
