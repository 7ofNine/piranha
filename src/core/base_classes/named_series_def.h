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

#ifndef PIRANHA_NAMED_SERIES_DEF_H
#define PIRANHA_NAMED_SERIES_DEF_H

#include <complex>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <boost/functional/hash.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/unordered_map.hpp>

#include "../config.h"
#include "../exceptions.h"
#include "../mp.h"
#include "../ntuple.h"
#include "../Psym.h"
#include "../settings.h"
#include "../type_traits.h"
#include "base_series_def.h"
#include "named_series_mp.h"

#define PIRANHA_NAMED_SERIES_TP_DECL class ArgsDescr, class Term, class Derived
#define PIRANHA_NAMED_SERIES_TP            ArgsDescr,       Term,       Derived

namespace piranha
{
	/// Dictionary for evaluation.
	typedef boost::unordered_map<std::string, double> EvalDict;

	/// Named series toolbox.
	/**
	 * Toolbox for generating series with arguments.
	 * ArgsDescr must be a boost::tuple of structures each one containing a static const string
	 * called "name" naming the arguments of the series.
	 */

     //types for managing layouts of named series arguments. For details of layouts see named_series_manip.h
     typedef std::pair<bool, std::size_t> LayoutElement;
     typedef std::vector<LayoutElement> Layout;

     template<class ArgsTuple>
     class LayoutTuple
     {
        public:
        typedef typename NTuple< Layout, boost::tuples::length<ArgsTuple>::value >::Type Type;   
     };
     //typedef typename NTuple< Layout, boost::tuples::length<ArgsTuple>::value >::Type LayoutTuple;



	template <PIRANHA_NAMED_SERIES_TP_DECL>
	class NamedSeries
	{
			template <class T, class Enable>
			friend class NamedSeriesAddSelector;

			template <class T, class Enable>
			friend class NamedSeriesSubtractSelector;

			template <class T, class Enable>
			friend class NamedSeriesMultiplySelector;

			template <class T, class Enable>
			friend class NamedSeriesEqualitySelector;

		public:
			typedef ArgsDescr ArgumentsDescription; //ArgsDescr is expected to be a boost::tuples (NTuple) of echelonLevel+1 length. 
			typedef typename NTuple<VectorPsym, boost::tuples::length<ArgumentsDescription>::value>::Type ArgsTupleType;

		private:
			class SeriesIteratorGenerator
			{
                public:

				typedef Derived result_type; // don't change. required by boost::transform_iterator

				Derived operator()(Term const &t) const
				{
					Derived retval;
					retval.argumentsTuple = series.argumentsTuple;
					retval.insert(t, retval.argumentsTuple);
					return retval;
				}

				SeriesIteratorGenerator(Derived const &series):series(series) {}

                private:

				const Derived &series;
			};

            class CompareTrig
            {
            public:
                explicit CompareTrig(std::vector<std::pair<bool, std::size_t> > const & positions) :positions(positions) {}

                bool operator()(Term const * const t1, Term const * const t2) const
                {
                    //TODO: flavour?? cosine before sine which makes constants appear at the top
                    auto const trigkey1 = t1->get<1>();  //TODO: make this actually safe // a vector of exponents How do we make sure it is the trig arguments we want, Template parameter?
                    auto const trigkey2 = t2->get<1>();  // 1: should be the trig arguments
                    PIRANHA_ASSERT(positions.size() == trigkey1.size());
                    PIRANHA_ASSERT(positions.size() == trigkey2.size());
                    for (decltype(positions.size()) i = 0, e = positions.size(); i < e; ++i)
                    {
                        if (positions[i].first)
                        {
                            auto const currentIndex = positions[i].second;
                            if (trigkey1[currentIndex] < trigkey2[currentIndex])
                            {
                                return true;
                            }
                            else if (trigkey1[currentIndex] > trigkey2[currentIndex])
                            {
                                return false;
                            }
                            // they are equal, check the next one
                        }
                    }
                    return false; // all where the same
                }


            private:
                std::vector<std::pair<bool, std::size_t> > const & positions;
            };


		public:

			typedef boost::transform_iterator<SeriesIteratorGenerator, typename SeriesContainer<Term>::Type::const_iterator> SeriesIterator;
			SeriesIterator itBegin() const;
			SeriesIterator itEnd() const;
			std::complex<Derived> complex() const;

			void print(std::ostream &stream = std::cout) const;
			void printPlain(std::ostream &) const;
			void printTex(std::ostream &) const;

			void saveTo(const std::string &) const;
			void printToSorted(std::string const & fileName, VectorPsym const & expSymbols, VectorPsym const & trigSymbols) const;


            typedef std::vector<typename Term::KeyType> PrintSequenceType;
            void printToSequenced(std::string const & fileName, VectorPsym const & expSymbols, VectorPsym const & trigSymbols, PrintSequenceType const & sequence) const;

//		    Rework this.
// 			template <class Filter>
// 			Derived filter(const Filter &) const;
			
            void swap(Derived &);
			double norm() const;
			typename TermEvalTypeDeterminer<Term>::Type eval(double const &) const;
			typename TermEvalTypeDeterminer<Term>::Type eval(EvalDict const &) const;
			
            int psi(int const start = 0, int const step = 1) const;
			
            const ArgsTupleType &arguments() const;
			
            void setArguments(ArgsTupleType const &);

			template <class T>
			bool operator==(T const &) const;

			template <class T>
			bool operator!=(T const &) const;

			template <class T>
			Derived &operator+=(T const &);

			template <class T>
			Derived &operator-=(T const &);

			Derived operator-() const;

			template <class T>
			Derived &operator*=(T const &);

			template <class T>
			Derived &operator/=(T const &);

			Derived pow(double const) const;
			Derived pow(mp_rational const &) const;
			Derived root(int const) const;
			Derived partial(std::string const &, int const n = 1) const;

			template <class SubstitutionSeries>
			Derived sub(std::string const &, SubstitutionSeries const &) const;

			template <class Key>
			Derived seriesFromKey(Key const &) const;

			template <class Cf>
			Derived seriesFromCf(Cf const &) const;

			std::vector<std::vector<Derived> > split(int const n = 0) const;

			std::vector<Derived> flatten() const;
			
            ~NamedSeries();
			
            //protected:
			void trim(); // remove all arguments that are no longer in the series from the argumentsTuple
			
            template <class Derived2>
			void mergeArgs(Derived2 const &);
			
            // TODO: check these protected methods, some of them can be moved into private
			// with proper friendship in manipulator classes.
			void constructFromFile(std::string const &);

			template <int N>
			void constructFromPsym(Psym const &);

		private:

			template <class T>
			bool isEqualTo(T const &) const;

			void appendArg(std::string const &, Psym const &);

			template <int N>
			void appendArg(Psym const &);

			template <class Derived2>
			Derived & multiplyBySeries(Derived2 const &);

			template <class Number>
			Derived & multiplyNumberHelper(Number const &);

			template <class Number>
			Derived & divideNumberHelper(Number const &);

			void printPretty(std::ostream &) const;
			void readFromFile(std::ifstream &, std::string const &);
			void readSections(std::ifstream &);
			void readArg(std::ifstream &, std::string const &);
			void readTerms(std::ifstream &);

			template <class Derived2>
			bool isArgsCompatible(Derived2 const &) const;

			template <class Derived2>
			void mergeIncompatibleArgs(Derived2 const &);

			template <bool, class Derived2>
			Derived & mergeWithSeries(Derived2 const &);

			std::vector<Term const*> getTrigSortedTerms(std::vector<std::pair<bool, std::size_t> > const & positions) const; // for printing
			void printToSorted(std::ofstream & outfile, std::vector<Term const*> trigSortedTerms, VectorPsym const & expSymbols, VectorPsym const & trigSymbols, std::vector<std::pair<bool, std::size_t> > const & trigPositions) const; //for printing
			std::vector<typename Term::CfType::TermType const*> getExpoSortedCoefficient(typename Term::CfType const & coeff, std::vector<std::pair<bool, size_t> > const & expSymbols) const;//for printing return type??
            std::vector<Term const *> getTrigSequencedTerms(std::vector<std::pair<bool, std::size_t> > const &positions, PrintSequenceType const & sequence) const;


		protected:

			// Data members.
			ArgsTupleType                   argumentsTuple;  // the arguments of the NamedSeries. 
			static std::vector<std::string> unknownData;
	};

	// Initialization of static member.
	template <PIRANHA_NAMED_SERIES_TP_DECL>
	std::vector<std::string> NamedSeries<PIRANHA_NAMED_SERIES_TP>::unknownData;

// Useful macros for named series.
#define E0_SERIES_NAMED_ANCESTOR(Args, TermName, SeriesName) piranha::NamedSeries<Args, TermName, E0_SERIES(SeriesName)>

#define E1_SERIES_NAMED_ANCESTOR(Args1, Args2, TermName, SeriesName) piranha::NamedSeries<boost::tuple<Args1, Args2>, TermName, SeriesName>

#define NAMED_SERIES_BOILERPLATE(SeriesName, N) \
public: \
	explicit SeriesName() {} \
	explicit SeriesName(const piranha::Psym &p) { \
		this->template constructFromPsym<N>(p); \
	} \
	explicit SeriesName(const std::string &filename) \
	{ \
		this->constructFromFile(filename); \
	} \
	explicit SeriesName(const double &x) \
	{ \
		*this = this->baseSeriesFromNumber(x, this->argumentsTuple); \
		this->trim(); \
	} \
	explicit SeriesName(const piranha::mp_rational &q) \
	{ \
		*this = this->baseSeriesFromNumber(q, this->argumentsTuple); \
		this->trim(); \
	} \
	explicit SeriesName(const piranha::mp_integer &z) \
	{ \
		*this = this->baseSeriesFromNumber(z, this->argumentsTuple); \
		this->trim(); \
	}
}

#endif
