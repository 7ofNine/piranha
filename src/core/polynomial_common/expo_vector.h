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

			typedef VectorKey<__PIRANHA_EXPO_VECTOR_TP, ExpoVector<__PIRANHA_EXPO_VECTOR_TP> > ancestor;

			template <class SubSeries, class ArgsTuple>
			class SubCache: public PowerCache<SubSeries, T, BaseSeriesArithmetics<SubSeries, ArgsTuple> >
			{
					typedef PowerCache<SubSeries, T, BaseSeriesArithmetics<SubSeries, ArgsTuple> > ancestor;

				public:

					SubCache():ancestor() {}

					void setup(const SubSeries &s, const ArgsTuple *argsTuple)
					{
						this->m_arith_functor.argsTuple = argsTuple;
						// NOTE: move semantics here.
						SubSeries tmp;
						tmp.baseAdd(1, *argsTuple);
						this->m_container[T(0)] = tmp;
						this->m_container[T(1)] = s;
					}
			};

			// This is just a stub.
			template <class SubSeries, class ArgsTuple>
			class ei_sub_cache
			{
				public:

					void setup(const SubSeries &, const ArgsTuple *) {}
			};

			template <class, class>
			friend struct expo_vector_pow_double;

			template <class, class>
			friend struct expo_vector_pow_rational;

		public:

			typedef typename ancestor::value_type value_type;
			typedef value_type                    degree_type;
			typedef typename ancestor::size_type  size_type;
			typedef double                        EvalType;

			template <class SubSeries, class SubCachesCons, class ArgsTuple>
			struct sub_cache_selector 
            {
				typedef boost::tuples::cons<SubCache<SubSeries, ArgsTuple>, SubCachesCons> type;
			};

			template <class SubSeries, class SubCachesCons, class ArgsTuple>
			struct ei_sub_cache_selector 
            {
				typedef boost::tuples::cons<ei_sub_cache<SubSeries, ArgsTuple>, SubCachesCons> type;
			};


			// Ctors.
			/// Default ctor.
			ExpoVector(): ancestor() {}


			/// Ctor from string.
			template <class ArgsTuple>
			explicit ExpoVector(const std::string &s, const ArgsTuple &): ancestor()
			{
				std::vector<std::string> sd;
				boost::split(sd, s, boost::is_any_of(std::string(1, this->separator)));
				const size_type w = boost::numeric_cast<size_type>(sd.size());
				for (size_type i = 0; i < w; ++i)
				{
					this->container.push_back(boost::lexical_cast<value_type>(sd[i]));
				}
			}


			/// Ctor from Psym.
			template <class ArgsTuple>
			explicit ExpoVector(const Psym &p, const int &n, const ArgsTuple &a): ancestor(p, n, a) {}


			// Math.
			/// Multiplication.
			template <class ResultType>
			void multiply(const ExpoVector &e2, ResultType &ret) const
			{
				const size_type max_w = this->size();
				const size_type min_w = e2.size();
				// Resize, if needed.
				ret.resize(max_w);

				// Assert widths, *this should always come from a polynomial, and its width should hence be
				// already adjusted my merge_args in multiplication routines.
				PIRANHA_ASSERT(max_w >= min_w);
				PIRANHA_ASSERT(ret.size() == max_w);

				size_type i;
				for (i = 0; i < min_w; ++i) 
                {
					ret[i] = (*this)[i] + e2[i];
				}

				for (; i < max_w; ++i) 
                {
					ret[i] = (*this)[i];
				}
			}


			// I/O.
			template <class ArgsTuple>
			void print_plain(std::ostream &outStream, const ArgsTuple &argsTuple) const
			{
				PIRANHA_ASSERT(argsTuple.template get<ancestor::position>().size() == this->size());

				(void)argsTuple;
				this->printElements(outStream);
			}


			template <class ArgsTuple>
			void print_pretty(std::ostream &outStream, const ArgsTuple &argsTuple) const
			{
				PIRANHA_ASSERT(argsTuple.template get<ancestor::position>().size() == this->size());

				bool printed_something = false;
				for (size_type i = 0; i < this->size(); ++i) 
                {
					const value_type &n = (*this)[i];
					// Don't print anything if n is zero.
					if (n != 0) 
                    {
						// Prepend the multiplication operator only if we already printed something.
						if (printed_something) 
                        {
							outStream << '*';
						}

						// Take care of printing the name of the exponent.
						outStream << argsTuple.template get<ancestor::position>()[i].get_name();
						// Print the pow operator only if exponent is not unitary.
						if (n != 1) 
                        {
							outStream << "**";
							expo_vector_print_element_pretty(outStream, n);
						}
						printed_something = true;
					}
				}
			}


			template <class ArgsTuple>
			void print_tex(std::ostream &outStream, const ArgsTuple &argsTuple) const
			{
				PIRANHA_ASSERT(argsTuple.template get<ancestor::position>().size() == this->size());

				for (size_type i = 0; i < this->size(); ++i) 
                {
					const value_type &n = (*this)[i];
					// Don't print anything if n is zero.
					if (n != 0) 
                    {
						// Take care of printing the name of the exponent.
						outStream << ' ' << argsTuple.template get<ancestor::position>()[i].get_name() << ' ';
						// Print the pow operator only if exponent is not unitary.
						if (n != 1) 
                        {
							outStream << "^{";
							expo_vector_print_element_tex(outStream, n);
							outStream << '}';
						}
					}
				}
			}


			template <class ArgsTuple>
			double eval(const double &t, const ArgsTuple &argsTuple) const
			{
				const size_type w = this->size();
				PIRANHA_ASSERT(w <= argsTuple.template get<ancestor::position>().size());

				double retval = 1.;
				for (size_type i = 0; i < w; ++i) 
                {
					retval *= std::pow(argsTuple.template get<ancestor::position>()[i].eval(t), (*this)[i]);
				}

				return retval;
			}


			/// Am I ignorable?
			/**
			 * By construction an array of exponents is never ignorable.
			 */
			template <class ArgsTuple>
			bool is_ignorable(const ArgsTuple &) const
			{
				return false;
			}


			bool is_unity() const
			{
				return (this->elementsAreZero());
			}


			/// Comparison operator.
			bool operator<(const ExpoVector &e2) const
			{
				return this->lex_comparison(e2);
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
			std::size_t hash_value() const
			{
				return this->elementsHasher();
			}


			/// Return the total degree of the exponents array.
			degree_type degree() const
			{
				return std::accumulate(this->begin(), this->end(), degree_type(0));
			}


			/// Total degree of the variables at specified positions pos.
			/**
			 * pos_tuple must be a tuple of vectors of (bool,std::size_t) pairs.
			 */
			template <class PosTuple>
			degree_type partial_degree(const PosTuple &pos_tuple) const
			{
				// TODO: check PosTuple type.
				const std::vector<std::pair<bool,std::size_t> > &pos = pos_tuple.template get<ancestor::position>();
				const size_type w = this->size(); 
                const size_type pos_size = boost::numeric_cast<size_type>(pos.size());
				degree_type retval(0);
				for (size_type i = 0; i < pos_size; ++i) 
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
			degree_type order() const
			{
				return degree();
			}


			/// Minimum total degree of the variables at specified positions pos.
			/**
			 * Provided for use within the power series toolbox, and defined to be equivalent to partial_degree().
			 */
			template <class PosTuple>
			degree_type partial_order(const PosTuple &pos_tuple) const
			{
				return partial_degree(pos_tuple);
			}


			/// Calculate partial derivative.
			template <class Series, class PosTuple, class ArgsTuple>
			Series partial(const PosTuple &pos_tuple, const ArgsTuple &argsTuple) const
			{
				// TODO: check PosTuple type.
				PIRANHA_ASSERT(pos_tuple.template get<ancestor::position>().size() == 1);

				const std::size_t pos = pos_tuple.template get<ancestor::position>()[0].second;

				PIRANHA_ASSERT(!pos_tuple.template get<ancestor::position>()[0].first || pos < this->size());

				// Do something only if the argument of the partial derivation is present in the exponent array
				// and the interesting exponent is not zero.
				Series retval;
				if (pos_tuple.template get<ancestor::position>()[0].first && (*this)[boost::numeric_cast<size_type>(pos)] != 0) 
                {
					ExpoVector copy(*this);
					// NOTE: move semantics here?
					copy[pos] -= 1;
					retval = Series::baseSeriesFromKey(copy,argsTuple);
					retval.baseMultBy((*this)[pos],argsTuple);
				}
				return retval;
			}


			template <class ArgsTuple>
			ExpoVector pow(const double &y, const ArgsTuple &) const
			{
				if (is_integer(y)) 
                {
					return generic_pow((int)y);
				} else 
                {
					return expo_vector_pow_double<value_type>::run(*this,y);
				}
			}


			template <class ArgsTuple>
			ExpoVector pow(const mp_rational &q, const ArgsTuple &) const
			{
				return expo_vector_pow_rational<value_type>::run(*this,q);
			}


			template <class RetSeries, class PosTuple, class SubCaches, class ArgsTuple>
			RetSeries sub(const PosTuple &posTuple, SubCaches &subCaches, const ArgsTuple &argsTuple) const
			{
				// TODO: check on PosTuple.
				RetSeries retval;
				// If the argument is not present here, the return series will have one term consisting
				// of a unitary coefficient and this very ExpoVector.
				// NOTE: for now we can substitute one symbol at a time.
				PIRANHA_ASSERT(posTuple.template get<ancestor::position>().size() == 1);

				if (!posTuple.template get<ancestor::position>()[0].first) 
                {
					retval = RetSeries::baseSeriesFromKey(*this, argsTuple);

				} else 
                {
					const std::size_t pos = posTuple.template get<ancestor::position>()[0].second;

					PIRANHA_ASSERT(pos < this->size());
					
                    ExpoVector tmp_ea(*this);
					// Let's turn off the exponent associated to the symbol we are substituting.
					tmp_ea[pos] = 0;
					RetSeries orig(RetSeries::baseSeriesFromKey(tmp_ea, argsTuple));

					PIRANHA_ASSERT(retval.empty());
					
                    // NOTICE: series multadd here?
					retval.baseAdd(orig, argsTuple);
					retval.baseMultBy(subCaches.template get<ancestor::position>()[(*this)[pos]], argsTuple);
				}
				return retval;
			}


			template <class RetSeries, class PosTuple, class SubCaches, class ArgsTuple>
			RetSeries ei_sub(const PosTuple &, SubCaches &, const ArgsTuple &argsTuple) const
			{
				return RetSeries::baseSeriesFromKey(*this, argsTuple);
			}


		private:

			// Generic exponentiation.
			template <class U>
			ExpoVector generic_pow(const U &x) const
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
