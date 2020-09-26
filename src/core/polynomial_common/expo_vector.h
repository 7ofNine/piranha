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

#ifndef PIRANHA_EXPO_VECTOR_H
#define PIRANHA_EXPO_VECTOR_H

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/type_traits/integral_constant.hpp>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "../base_classes/vector_key.h"
#include "../common_functors.h"
#include "../exceptions.h"
#include "../power_cache.h"
#include "../Psym.h"
#include "../type_traits.h"
#include "expo_vector_mp.h"

#define __PIRANHA_EXPO_VECTOR_TP_DECL class T, int Pos
#define __PIRANHA_EXPO_VECTOR_TP T, Pos

namespace piranha
{
	/// Exponents vector.
	template < __PIRANHA_EXPO_VECTOR_TP_DECL >
	class ExpoVector: public VectorKey<__PIRANHA_EXPO_VECTOR_TP, ExpoVector<__PIRANHA_EXPO_VECTOR_TP> >
	{
		private:

			typedef VectorKey<__PIRANHA_EXPO_VECTOR_TP, ExpoVector<__PIRANHA_EXPO_VECTOR_TP> > Ancestor;

			template <class SubSeries, class ArgsTuple>
			class SubCache: public PowerCache<SubSeries, T, BaseSeriesArithmetics<SubSeries, ArgsTuple> >
			{
					typedef PowerCache<SubSeries, T, BaseSeriesArithmetics<SubSeries, ArgsTuple> > Ancestor;

				public:

					SubCache():Ancestor() {}

					void setup(const SubSeries &s, const ArgsTuple *argsTuple)
					{
						this->arithmeticFunctor.argsTuple = argsTuple;
						// NOTE: move semantics here.
						SubSeries tmp;
						tmp.baseAdd(1, *argsTuple);
						this->container[T(0)] = tmp;
						this->container[T(1)] = s;
					}
			};

			// This is just a stub.
			template <class SubSeries, class ArgsTuple>
			class EiSubCache
			{
				public:

					void setup(const SubSeries &, const ArgsTuple *) {}
			};

			template <class, class>
			friend struct ExpoVectorPowDouble;

			template <class, class>
			friend struct ExpoVectorPowRational;

		public:

			typedef typename Ancestor::value_type value_type;
			typedef value_type                    DegreeType;
			typedef typename Ancestor::size_type  size_type;
			typedef double                        EvalType;

			template <class SubSeries, class SubCachesCons, class ArgsTuple>
			struct SubstitutionCacheSelector 
            {
				typedef boost::tuples::cons<SubCache<SubSeries, ArgsTuple>, SubCachesCons> Type;
			};

			template <class SubSeries, class SubCachesCons, class ArgsTuple>
			struct EiSubstitutionCacheSelector 
            {
				typedef boost::tuples::cons<EiSubCache<SubSeries, ArgsTuple>, SubCachesCons> Type;
			};


			// Ctors.
			/// Default ctor.
			ExpoVector(): Ancestor() {}


			/// Ctor from string of exponent values.
			// what is the ArgsTuple good for. The association between argstuple and acutal value is done by the using class
			// the created vector of exponents has the length of different detected values in the string 's' and is not related to ArgsTuple
			// construct an expo vector from a string of numeric values. The relation between symbol and position of the 
			// (index) of the numerical value is up to the user.
			// Is this actually used anywhere?
			template <class ArgsTuple>
			explicit ExpoVector(const std::string &s, const ArgsTuple &): Ancestor()
			{
				PIRANHA_ASSERT(!s.empty());
				std::vector<std::string> sd;
				boost::split(sd, s, boost::is_any_of(std::string(1, this->separator)));
				const size_type w = boost::numeric_cast<size_type>(sd.size());

				for (size_type i = 0; i < w; ++i)
				{
					this->container.push_back(boost::lexical_cast<value_type>(boost::algorithm::trim_copy(sd[i])));
				}
			}


			/// Ctor from Psym.
			// construct an expoVector from a Psym. There will only be one element in the vector.
			// n is the position of the key in the key hierarchy (echelon level)
			template <class ArgsTuple>
			explicit ExpoVector(const Psym &p, const int n, const ArgsTuple &a): Ancestor(p, n, a) {}


			// Math.
			// Multiplication.
            // multiply this ExpoVector with expoVector into result
            // result is allocated externally
            // the exponents have to represent the same symbol in order to make sense. This
            // has to asserted before one can use this method. See merge arguments
			// This is asymmetric. The size odf this  has to be bigger than the size of
			// the input expoVector.
			// Should that be made symmetric?
			template <class ResultType>
			void multiply(const ExpoVector &expoVector, ResultType &result) const
			{
				const size_type maxw = this->size();
				const size_type minw = expoVector.size();
				// Resize, if needed.
				result.resize(maxw);

				// Assert widths, *this should always come from a polynomial, and its width should hence be
				// already adjusted my mergeArgs in multiplication routines.
				PIRANHA_ASSERT(maxw >= minw);
				PIRANHA_ASSERT(result.size() == maxw);

                size_type i;
				for (i = 0; i < minw; ++i) 
                {
					result[i] = (*this)[i] + expoVector[i]; //multiplying means adding exponents
				}

                //i comes from previous for loop and starts with i == minw
                // Here it is where the merge of different length expo vectors is happening.
				for (; i < maxw; ++i) 
                {
					result[i] = (*this)[i]; 
				}
			}


			// I/O.
			template <class ArgsTuple>
			void printPlain(std::ostream &outStream, const ArgsTuple &argsTuple) const
			{
				PIRANHA_ASSERT(argsTuple.template get<Ancestor::position>().size() == this->size());
                (void)argsTuple;

				this->printElements(outStream);
			}

			template <class ArgsTuple>
			void printPlainSorted(std::ostream & outStream, std::vector<std::pair<bool, std::size_t> > positions, ArgsTuple const & argsTuple) const
			{
				PIRANHA_ASSERT(argsTuple.template get<Ancestor::position>().size() == this->size());
				(void)argsTuple;

				this->printElementsSorted(outStream, positions);

			}



			// index in argsTuple and ExpoVector are already assumed to be the same i.e.
			// they correspond!
			// what about negative exponents
			// rational are handled by expoVectorPrintElementPretty
			template <class ArgsTuple>
			void printPretty(std::ostream &outStream, const ArgsTuple &argsTuple) const
			{
				PIRANHA_ASSERT(argsTuple.template get<Ancestor::position>().size() == this->size());

				bool printedSomething = false;

				for (size_type i = 0; i < this->size(); ++i) 
                {
					const value_type &n = (*this)[i];
					// Don't print anything if n is zero.
					if (n != 0) 
                    {
						// Prepend the multiplication operator only if we already printed something.
						if (printedSomething) 
                        {
							outStream << '*';
						}

						// Take care of printing the name of the exponent.
						outStream << argsTuple.template get<Ancestor::position>()[i].getName();

						// Print the pow operator only if exponent is not unitary.
						if (n != 1) 
                        {
							outStream << "**";
							expoVectorPrintElementPretty(outStream, n);
						}
						printedSomething = true;
					}
				}
			}


			template <class ArgsTuple>
			void printTex(std::ostream &outStream, const ArgsTuple &argsTuple) const
			{
				PIRANHA_ASSERT(argsTuple.template get<Ancestor::position>().size() == this->size());

				for (size_type i = 0; i < this->size(); ++i) 
                {
					const value_type &n = (*this)[i];
					// Don't print anything if n is zero.
					if (n != 0) 
                    {
						// Take care of printing the name of the exponent.
						outStream << ' ' << argsTuple.template get<Ancestor::position>()[i].getName() << ' ';
						// Print the pow operator only if exponent is not unitary.
						if (n != 1) 
                        {
							outStream << "^{";
							expoVectorPrintElementTEX(outStream, n);
							outStream << '}';
						}
					}
				}
			}

			// just what we expect. The numerical evaluation of the exponential vector 
			// respecting the parameters time dependency
			// the value is 1 if no exponent present
			// the index of the exponent vector and the argsTuple have to be coordinated externally
			// the number of parameters in the ExpoVector has to be smaller or equal to the the number in the argsTuple

			template <class ArgsTuple>
			double eval(const double t, const ArgsTuple &argsTuple) const
			{
				const size_type w = this->size();
				PIRANHA_ASSERT(w <= argsTuple.template get<Ancestor::position>().size());

				double retval = 1.0;
				for (size_type i = 0; i < w; ++i) 
                {
                    //get the numerical value for the parameter (may be time dependent) and calculate its power according to this ExpoVector
					retval *= std::pow(argsTuple.template get<Ancestor::position>()[i].eval(t), (*this)[i]);
				}

				return retval;
			}


			/// Am I ignorable?
			/**
			 * By construction an array of exponents is never ignorable.
			 */
			template <class ArgsTuple>
			bool isIgnorable(const ArgsTuple &) const
			{
				return false;
			}


			bool isUnity() const
			{
                //its a 1 if all exponante are 0
				return this->elementsAreZero();
			}


			/// Comparison operator.
			bool operator<(const ExpoVector &expoVector) const
			{
				return this->lexComparison(expoVector);
			}


			/// Norm.
			/**
			 * The norm of an exponent array is defined as the absolute value of the evaluation at t = 0.
			 */
			template <class ArgsTuple>
			double norm(const ArgsTuple &argsTuple) const
			{
				return std::abs(eval(0, argsTuple));
			}


			/// Calculate hash value.
            /// Don't rename. Might be needed by boost
			std::size_t hash_value() const
			{
				return this->elementsHasher();
			}


			/// Return the total degree of the exponents array, i.e. sum over all the exponent values
			DegreeType degree(VectorPsym const & symbols) const
			{
                // the length of the symbols can be >= the size of the current vector.
                // this happens during merge of arguments in multiplication and addition.
                // The positions have to correspond to the correct symbols, but this can not be verified here
                PIRANHA_ASSERT(symbols.size() >= this->size())
				DegreeType degree(0);
				for (VectorPsym::size_type i = 0; i < this->size(); ++i)
				{
					degree += symbols[i].order() * (*this)[i];
				}
				return degree;

				//return std::accumulate(this->begin(), this->end(), DegreeType(0));
			}


			/// Total degree of the variables at specified positions pos.
			/**
			 * pos_tuple must be a tuple of vectors of (bool,std::size_t) pairs.
			 */
			template <class PosTuple>
			DegreeType partialDegree(const PosTuple &posTuple) const
			{
				// TODO: check PosTuple type.
				const std::vector<std::pair<bool, std::size_t> > &pos = posTuple.template get<Ancestor::position>();
				const size_type w = this->size(); 
                const size_type posSize = boost::numeric_cast<size_type>(pos.size());
				DegreeType retval(0);

				for (size_type i = 0; i < posSize; ++i) 
                {
					// Add up exponents only if they are present and don't try to read outside boundaries
					// (this last could happen after merging arguments with a second series with smaller
					// number of arguments).
					if (pos[i].first && pos[i].second < w) 
                    {
						retval += (*this)[boost::numeric_cast<size_type>(pos[i].second)];
					}
				}
				return retval;
			}


			/// Minimum total degree of the exponents array.
			/**
			 * Provided for use within the power series toolbox, and defined to be equivalent to degree().
			 */
			DegreeType order(VectorPsym const & symbols) const
			{
				return degree(symbols);
			}

			DegreeType xorder() const
			{
				DegreeType degree(0);
				for (int i = 0; i < this->size(); ++i)
				{
					degree += (*this)[i];
				}
				return degree;

			}
			/// Minimum total degree of the variables at specified positions pos.
			/**
			 * Provided for use within the power series toolbox, and defined to be equivalent to partialDegree().
			 */
			template <class PosTuple>
			DegreeType partialOrder(const PosTuple &posTuple) const
			{
				return partialDegree(posTuple);
			}


			/// Calculate partial derivative.
			// where does Series and PosTuple actually come from.
			// Series is at least a BaseSeries, while PosTuple is only generated in NamedSeries
			// is this a good idea to transport this to this layer?
			// PosTuple may be as a defined class (not a parameter). SHould this tuple not actually be an array?
			// PosTuple pair (bool,int): bool indicates if the originally indicated name is present
			//                           int: indicates the index into the vector of the corresponding element in the ArgsTuple
			// Note : PosTuple is provided from the outside, the resolution from name/psym is already done

			// SHould that be 0 if postuple contains no found element??
			// mathematically yes and I haven't found surrounding usage that would would correct for that
			// Guess it just hasn't happened or somehwere is a check to interpret an empty series as 0.
			// That is not nice becasue we would not be homogenous!! and consistent i.e. two ways of representing 0.
			template <class Series, class PosTuple, class ArgsTuple>
			Series partial(const PosTuple &posTuple, const ArgsTuple &argsTuple) const
			{
				// TODO: check PosTuple type.
				PIRANHA_ASSERT(posTuple.template get<Ancestor::position>().size() == 1);

				const std::size_t pos = posTuple.template get<Ancestor::position>()[0].second;

				PIRANHA_ASSERT(!posTuple.template get<Ancestor::position>()[0].first || pos < this->size());

				// Do something only if the argument of the partial derivation is present in the exponent array
				// and the interesting exponent is not zero.
				Series retval;
				if (posTuple.template get<Ancestor::position>()[0].first && (*this)[boost::numeric_cast<size_type>(pos)] != 0) 
                {
					ExpoVector copy(*this);
					// NOTE: move semantics here?
					copy[pos] -= 1;
					retval = Series::baseSeriesFromKey(copy, argsTuple);
					retval.baseMultBy((*this)[pos], argsTuple);
				}
				return retval;
			}


			//ArgsTuple is nowhere used. Why is it here??? 
			// which combinations do we allow. rational exponents, double exponents?
			//
			// control over the parameter y should probably better
			template <class ArgsTuple>
			ExpoVector pow(const double &y, const ArgsTuple &) const
			{
				if (is_integer(y)) 
                {
					return genericPow((int)y);
				} else 
                {
					return ExpoVectorPowDouble<value_type>::run(*this, y);
				}
			}


			template <class ArgsTuple>
			ExpoVector pow(const mp_rational &q, const ArgsTuple &) const
			{
				return ExpoVectorPowRational<value_type>::run(*this, q);
			}

			// Again a good idea to have that here. Series is a much higher object than an ExpoVector
			// substituting what?? Something from the SubCache
			// ALlows for a single element in position tuple . Position tuple has to come from the outside
			//
			template <class RetSeries, class PosTuple, class SubCaches, class ArgsTuple>
			RetSeries sub(const PosTuple &posTuple, SubCaches &subCaches, const ArgsTuple &argsTuple) const
			{
				// TODO: check on PosTuple.
				RetSeries retval;
				// If the argument is not present here, the return series will have one term consisting
				// of a unitary coefficient and this very ExpoVector.
				// NOTE: for now we can substitute one symbol at a time.
				PIRANHA_ASSERT(posTuple.template get<Ancestor::position>().size() == 1);

				if (!posTuple.template get<Ancestor::position>()[0].first) 
                {
					retval = RetSeries::baseSeriesFromKey(*this, argsTuple);

				} else 
                {
					const std::size_t pos = posTuple.template get<Ancestor::position>()[0].second;

					PIRANHA_ASSERT(pos < this->size());
					
                    ExpoVector tmpEa(*this);
					// Let's turn off the exponent associated to the symbol we are substituting.
					tmpEa[pos] = 0;
					RetSeries orig(RetSeries::baseSeriesFromKey(tmpEa, argsTuple));

					PIRANHA_ASSERT(retval.empty());
					
                    // NOTICE: series multadd here?
					retval.baseAdd(orig, argsTuple);
					retval.baseMultBy(subCaches.template get<Ancestor::position>()[(*this)[pos]], argsTuple);
				}
				return retval;
			}

			// IS that actually complete and or used anywhere ?  It just creates a sereis from the key
			// PosTuple is not used, Subcaches is not used
			template <class RetSeries, class PosTuple, class SubCaches, class ArgsTuple>
			RetSeries eiSubstitute(const PosTuple &, SubCaches &, const ArgsTuple &argsTuple) const
			{
				return RetSeries::baseSeriesFromKey(*this, argsTuple);
			}


		private:

			// Generic exponentiation.
			template <class U>
			ExpoVector genericPow(const U &x) const
			{
				ExpoVector retval(*this);
				const size_type w = this->size();

				for (size_type i = 0; i < w; ++i) 
                {
					retval[i] *= x;
				}

				return retval;
			}
	};


	/// is_ring_exact type trait specialisation for ExpoVector.
	template <__PIRANHA_EXPO_VECTOR_TP_DECL>
	struct is_ring_exact<ExpoVector<__PIRANHA_EXPO_VECTOR_TP> >: boost::true_type {};
}

#undef __PIRANHA_EXPO_VECTOR_TP_DECL
#undef __PIRANHA_EXPO_VECTOR_TP

#endif
