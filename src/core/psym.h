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

#ifndef PIRANHA_PSYMBOL_H
#define PIRANHA_PSYMBOL_H

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <set>
#include <string>
#include <utility> // For std::pair.
#include <vector>
#include <cmath>

#include <boost/lexical_cast.hpp>

#include "config.h"
#include "exceptions.h"
#include "ntuple.h"
#include "settings.h"
#include "utils.h"

namespace piranha
{
	// the global "dictionary" of Psym(s)
	struct  PsymManager
    {

        //implementation of a Psym stored in the PsymManager
		class  PsymImpl
        {
            friend class Psym;

            public:

			PsymImpl(std::string const &name, std::vector<double> const &timeEval = std::vector<double>(), int order = 1)
                    : name(name), timeEval(timeEval), order(order) {}
			
			// construct internal Psymbol of order order
			PsymImpl(std::string const &name, const int order)
				    :name(name), timeEval(std::vector<double>()), order(order) {}

			// Print to stream.
			void print(std::ostream &outStream) const
			{
				outStream << "name=" << name << '\n' << "time_eval=";

				for (std::size_t j = 0; j < timeEval.size(); ++j) 
				{
					outStream << boost::lexical_cast<std::string>(timeEval[j]);  // IS that really a good idea??? It creates funny digits
					if (j != timeEval.size() - 1) 
					{
						outStream << separator;
					}
				}

				outStream << '\n';
				outStream << "order=" << order << "\n";
			}

			// This is needed for storing them in a set. The set is sorted by name
			bool operator<(PsymImpl const &other) const
			{
				return (name < other.name);
			}


			// Time evaluation.
			double eval(double const &t) const
			{
				double retval = 0.0;
				const std::size_t w = timeEval.size();
				for (std::size_t i = 0; i < w; ++i) 
				{
					// FIXME: use natural_pow or the like here, to speed up?
					retval += std::pow(t, (int)i) * timeEval[i];
				}

				return retval;
			}

			void setOrder(const int order)
			{
				this->order = order;
			}

            private:
           
            // must not be changeable. The sorting of the set relies on it 
			const std::string		    name;     // the symbol name
			                                      // Mutable because we want to be able to freely change it in the Psym manager.
			mutable std::vector<double>	timeEval; // the time evolution polynomial. The order is the index into the vector 

			mutable int order;                    // the order of the symbol. Used for truncations. to make truncations specific to the symbol
                                                  // do we have to save it? 
            static const std::string	separator; // separator between timeValue elements. Leave in this position for readability during debug
			
		};

		// PsymManager global storage variable
		typedef std::set<PsymImpl>  ContainerType; // the set is sorted by name in the PsymImpl
		static ContainerType        container;
	};

	// Forward declaration for use in the typedef below.
	class Psym;

	typedef std::vector<Psym> VectorPsym;

	/// Literal symbol class.
	/**
	 * This class is used represent symbolic arguments. It features a string representing the symbol's name and a
	 * numerical vector
	 * which is used to evaluate the symbol in time in a polynomial fashion. For instance, if the numerical vector
	 * has a size of three and its elements are named \f$ \alpha \f$, \f$ \beta \f$ and \f$ \gamma \f$,
	 * it means that the symbol is evaluated as
	 * \f[
	 * \alpha + \beta t + \gamma t^2,
	 * \f]
	 * where \f$ t \f$ is time.
	 */
	class PIRANHA_VISIBLE Psym
    {

			typedef PsymManager::ContainerType::iterator Iterator;
			typedef PsymManager::PsymImpl                      PsymImpl;

			struct PushBackTo 
			{
				PushBackTo(VectorPsym &v):vectorPsym(&v) {}


				template <class T>
				void operator()(const T &x) const
				{
					vectorPsym->push_back(Psym(x.name, x.timeEval));
				}


				mutable VectorPsym *vectorPsym;
			};

		public:

			/// Constructor from name and time evaluation in string form.
			Psym(std::string const &name, std::string const &timeEval, int order = 1)
			{
				constructFromImpl(PsymImpl(name, utils::str_to_vector<double>(timeEval, PsymImpl::separator), order));
			}


			/// Constructor from name.
			/**
			 * If the symbol is already present in the Psym manager then use it, otherwise initialise
			 * a new Psym with the given name and an empty time evaluation vector.
			 * possiby give an order of the name (used for truncation)
			 */
			explicit Psym(std::string const &name, const int order = 1)
			{
				constructFromImpl(PsymImpl(name, order));
    //            // this is not constructFromImpl!
    //            PsymImpl pImpl(name, order);
    //            const Iterator itFound = PsymManager::container.find(pImpl);
				//if (itFound == PsymManager::container.end()) 
				//{
				//	std::pair<Iterator, bool> const res = PsymManager::container.insert(pImpl);
				//	PIRANHA_ASSERT(res.second);
				//	it = res.first;
				//} else 
				//{
				//	itFound->timeEval.
				//	itFound->order    = order;
				//	it = itFound;
				//}

			}


			/// Constructor from name and time value.
			Psym(std::string const &name, std::vector<double> const &timeEval)
			{
				constructFromImpl(PsymImpl(name, timeEval));
			}

			/// Constructor from name and time values and order.
			Psym(std::string const &name, std::vector<double> const &timeEval, unsigned int const order )
			{
				constructFromImpl(PsymImpl(name, timeEval, order));
			}


			/// Constructor from name and constant value.
			Psym(std::string const &name, double const &value)
			{
				constructFromImpl(PsymImpl(name, std::vector<double>(std::size_t(1), value)));
			}


			bool operator<(Psym const &other) const
			{
				return ((*it) < *(other.it));
			}


			bool operator==(Psym const &other) const
			{
				return it == other.it;
			}


			bool operator!=(Psym const &other) const
			{
				return !(*this == other);
			}


			static VectorPsym list()
			{
				VectorPsym retval;
				retval.reserve(PsymManager::container.size());

				std::for_each(PsymManager::container.begin(), PsymManager::container.end(), PushBackTo(retval));
				return retval;
			}


			// TODO: move to operator<< for streams? Along with other classes...
			/// Print to stream.
			void print(std::ostream &s = std::cout) const
			{
				it->print(s);
			}


			/// Evaluate symbol at time t.
			double eval(double const &t) const
			{
				return it->eval(t);
			}


			/// Name getter.
			const std::string &getName() const
			{
				return it->name;
			}


			/// Time evaluation vector getter.
			const std::vector<double> &getTimeEval() const
			{
				return it->timeEval;
			}

			// get order of the symbol
			int order() const
			{
				return it->order;
			}

			// set order to a new value
			// The underlying element is mutable. The set is not affected because it is sorted by name
			void setOrder(int const order)
			{
				it->order = order;
			}

			/// Time evaluation vector setter.
			// How is this possible without messing up the set.
			// the underlying field in the set element is mutable. The set is not affected because the set is ordered by name
			void setTimeEval(std::vector<double> const &t) const
			{
				it->timeEval = t;
			}
// TODO:     this does not work with the Python interfaces. How to fix
//            void set_time_eval(const std::string &t) const
//            {
//                
//                m_it->m_time_eval = utils::str_to_vector<double>(t, PsymImpl::separator);
//            }

		private:

			void constructFromImpl(PsymImpl const &p)
			{
				const Iterator itFound = PsymManager::container.find(p);
				if (itFound == PsymManager::container.end()) 
				{
					std::pair<Iterator, bool> const res = PsymManager::container.insert(p);
					PIRANHA_ASSERT(res.second);
					it = res.first;
				} else 
				{
					itFound->timeEval = p.timeEval; //overwrite time values
					itFound->order = p.order;       // overwrite order   
					it = itFound;
				}
			}

		private:

			// a Psym is an iterator into the set of PsymImpl managed in a global set in PsymManager . A PIMPL pattern
			Iterator it;
	};

	/////////////////////////////////////////////////////////////
	// the following methods are global in the piranha namespace
	/////////////////////////////////////////////////////////////

    // 
    // from a given vector of psyms return a vector of indicators where and if the elemnt in the vectorPsym is present in the
    // corresponding argsTuple component. ArgsTuple is a tuple of VectorPsyms, one vector per descriptortype (typ.: poly, trig -> 2-tuple)
	// for each descriptor in argsTuple a vector of pairs is returned.
	// -The size of the returned vectors is the size of the vector 'VectoPsym' entered.
    // -The index for vectorPsym element corresponds to the index in the vector in positionTuple. Elements of the vector in position
    // tuple are pairs (bool, int). 
    // (false , d) : not present in the argsTuple 
    // (true,   j) : present in the argsTuple vector at index j of the corresponding descriptor;
    //

	template <class ArgsTuple>
	class Psyms2posImpl
    {
        public:
		
        template <class PositionTuple>
		static void run(VectorPsym const &vectorPsym, PositionTuple &positionTuple, ArgsTuple const &argsTuple)
		{
			std::size_t const argsSize   = argsTuple.get_head().size();
			std::size_t const vectorSize = vectorPsym.size();

            positionTuple.get_head().reserve(vectorSize); // reserve result position vector in positionTuple
			
            // For each psymbol, test presence.
			for (std::size_t i = 0; i < vectorSize; ++i) 
			{
				// Initially set the symbol to not found.
				positionTuple.get_head().push_back(std::make_pair(false, std::size_t(0)));

				for (std::size_t j = 0; j < argsSize; ++j) 
				{
					if (argsTuple.get_head()[j] == vectorPsym[i]) 
					{
						positionTuple.get_head().back().first  = true; // present: update the last element in the positionVector
						positionTuple.get_head().back().second = j;    // position in argsTuple component
						// No need to continue, the symbol is (supposed to be) unique.
						break;
					}
				}
			}

            // got to next component in tuples
			Psyms2posImpl<typename ArgsTuple::tail_type>::run(vectorPsym, positionTuple.get_tail(), argsTuple.get_tail());
		}
	};


    //terminate recursion on tuple
	template <>
	class Psyms2posImpl<boost::tuples::null_type> 
	{
        public: 
		static void run(VectorPsym const &,  boost::tuples::null_type const &, boost::tuples::null_type const &) {}
	};


	// Transform a vector of psyms into a tuple of vectors of (position flag, position) pairs, given a reference arguments tuple.
	// Return value will be a tuple of vectors, each of size v.size(), containing (presence, position) pairs for the corresponding symbols
	// in v.
	template <class ArgsTuple>
	inline typename NTuple<std::vector<std::pair<bool, std::size_t> >, boost::tuples::length<ArgsTuple>::value>::Type  psyms2pos(VectorPsym const &vectorPsym, ArgsTuple const &argsTuple)
	{
		typedef typename NTuple<std::vector<std::pair<bool, std::size_t> >, boost::tuples::length<ArgsTuple>::value>::Type PositionTupleType;

        // TODO: Does this really work properly??
        // It is only used with vectors of one element i.e. find a symbols position in the argsTuple
        // or with partialDegree truncation which is doubtfull if that really works and is properly implemented. Haven't found out, yet or don't understand how that might actually be used??
        //
        // Problem is the reference from the vectorPsym to the retval vector.
        // std::set may change the sorting and because of the uniquenes in set retval may have fewer elements than vectorPsym???
        // 
        // I think the solution is that there is no reference back again to the original symbols. Only the position is of interest for further use, i.e. this would work
        // but is rather dangerous without documentation
        //
		// the field 'second' is only valid if 'first' is set to true. If set to 'false', 'second' is set to 0 wich would be a valid value

		// First we want to make sure that the vector of symbols does not contain duplicate elements.
		std::set<Psym> const uniquesSet(vectorPsym.begin(), vectorPsym.end());
		//set is only used to check on uniqueness. we don't want to destroy the sorting of the incomming symbols. This is used e.g. during printing
		PIRANHA_ASSERT(uniquesSet.size() == vectorPsym.size());
		//VectorPsym     const uniquesVector(uniquesSet.begin(), uniquesSet.end());
		PositionTupleType retval;

		Psyms2posImpl<ArgsTuple>::run(vectorPsym, retval, argsTuple);

		return retval;
	}


	// Transform a vector of names into a vector of symbols.
    // If name already exists the existing symbol is returned. If not a new symbol is created and returned
	inline VectorPsym names2psyms(std::vector<std::string> const &vectorNames)
	{
		const std::size_t size = vectorNames.size();
		VectorPsym vectorPsym;
		vectorPsym.reserve(size);

		for (std::size_t i = 0; i < size; ++i) 
		{
			vectorPsym.push_back(Psym(vectorNames[i]));
		}

		return vectorPsym;
	}
}
#endif
