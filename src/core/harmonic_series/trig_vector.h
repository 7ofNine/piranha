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

#ifndef PIRANHA_TRIG_VECTOR_H
#define PIRANHA_TRIG_VECTOR_H

#include <algorithm>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/functional/hash.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/type_traits/integral_constant.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iostream>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "../base_classes/vector_key.h"
#include "../common_functors.h"
#include "../exceptions.h"
#include "../mp.h"
#include "../power_cache.h"
#include "../type_traits.h"
#include "trig_vector_mp.h"

#define __PIRANHA_TRIG_VECTOR_TP_DECL class T, int Pos
#define __PIRANHA_TRIG_VECTOR_TP T,Pos

namespace piranha
{
	/// Trigonometric vector.
	template < __PIRANHA_TRIG_VECTOR_TP_DECL >
	class TrigVector: public VectorKey<__PIRANHA_TRIG_VECTOR_TP, TrigVector<__PIRANHA_TRIG_VECTOR_TP> >
	{
			typedef VectorKey<__PIRANHA_TRIG_VECTOR_TP, TrigVector<__PIRANHA_TRIG_VECTOR_TP> > Ancestor;

			template <class SubSeries, class ArgsTuple>
			class SubCache: public PowerCache<std::complex<SubSeries>, T, BaseSeriesArithmetics<std::complex<SubSeries>, ArgsTuple> >
			{
					typedef PowerCache<std::complex<SubSeries>, T, BaseSeriesArithmetics<std::complex<SubSeries>, ArgsTuple> > Ancestor;

                    //what is that good for?
					enum Status {
						zero,
						one,
						full
					};

				public:

					SubCache():Ancestor(), status(zero), errmsg() {}


					void setup(const SubSeries &s, const ArgsTuple *argsTuple)
					{
						this->arithmeticFunctor.argsTuple = argsTuple;
						this->container[T(0)] = std::complex<SubSeries>().baseAdd(1, *argsTuple);
						try {
							std::complex<SubSeries> tmp1(s.base_ei(*argsTuple));
							this->container[T(1)] = tmp1;
							status = one;

							SubSeries tmp2(s);
							tmp2.baseMultBy(-1, *argsTuple);

							std::complex<SubSeries> tmp3(tmp2.base_ei(*argsTuple));
							this->container[T(-1)] = tmp3;
							status = full;

						} catch (const value_error &ve) 
						{
							errmsg = ve.what();
						}
					}


					const std::complex<SubSeries> &operator[](const T &n)
					{
						switch (status) {
							case zero:
								if (n != 0) 
								{
									PIRANHA_THROW(value_error,std::string("the substitution cache was unable to "
										"compute the complex exponential of the series used for substitution. "
										"The reported error was:\n") + errmsg);
								}
							case one:
								if (n < 0) 
								{
									PIRANHA_THROW(value_error,std::string("the substitution cache was unable to "
										"compute the inverse complex exponential of the series used for substitution. "
										"The reported error was:\n") + errmsg);
								}
							default:
								;
						}
						return Ancestor::operator[](n);
					}

				private:

					Status		status;
					std::string	errmsg;
			};


			template <class SubSeries, class ArgsTuple>
			class EiSubCache: public PowerCache<SubSeries, T, BaseSeriesArithmetics<SubSeries, ArgsTuple> >
			{
					typedef PowerCache<SubSeries, T, BaseSeriesArithmetics<SubSeries, ArgsTuple> > Ancestor;

				public:

					EiSubCache():Ancestor() {}


					// NOTE: here we assume that s has absolute value equal to one, which lets us calculate its
					// inverse as conjugate. Note it into the documentation.
					void setup(const SubSeries &s, const ArgsTuple *argsTuple)
					{
						this->arithmeticFunctor.argsTuple = argsTuple;
						this->container[T(0)]  = SubSeries().baseAdd(1, *argsTuple);
						this->container[T(1)]  = s;
						this->container[T(-1)] = s.baseConjugate(*argsTuple);
					}
			};

		public:

			typedef typename Ancestor::value_type value_type;
			typedef value_type                    HarmonicDegreeType;
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
			TrigVector(): Ancestor(), flavour(true) {}


			/// Copy ctor from different position.
			template <int Position2>
			TrigVector(const TrigVector<T, Position2> &trigVector): Ancestor(), flavour(trigVector.getFlavour())
			{
				this->resize(trigVector.size());
				std::copy(trigVector.begin(), trigVector.end() ,this->begin());
			}


			/// Ctor from string.
			template <class ArgsTuple>
			explicit TrigVector(const std::string &s, const ArgsTuple &): Ancestor(), flavour(true)
			{
				std::vector<std::string> sd;
				boost::split(sd, s, boost::is_any_of(std::string(1, this->separator)));
				const size_type w = boost::numeric_cast<size_type>(sd.size());
				if (w == 0) 
				{
					// Flavour is already set to true.
					return;
				}

				// Now we know  w >= 1.
				this->resize(w - 1);
				for (size_type i = 0; i < w - 1; ++i) 
				{
					(*this)[i] = boost::lexical_cast<value_type>(sd[i]);
				}

				// Take care of flavour.
				if (*sd.back().c_str() == 's') 
				{
					flavour = false;

				} else if (*sd.back().c_str() != 'c') 
				{
					PIRANHA_THROW(value_error,"unknown flavour");
				}
			}


			template <class ArgsTuple>
			explicit TrigVector(const Psym &p, const int n, const ArgsTuple &a): Ancestor(p, n, a), flavour(true) {}


			// Math.
			/// Multiplication.
			/**
			 * Used in poisson_series_term multiplication.
			 * TODO: update docs below.
			 * Multiplication of two trigonometric functions using Werner's formulas, i.e.
			 * \f[
			 * C\cos\alpha\cdot\cos\beta=
			 * \frac{C}{2} \cos \left( \alpha - \beta \right) + \frac{C}{2} \cos \left( \alpha + \beta \right)
			 * \f]
			 * and the likes. Notice that in the first return value always goes the \f$ \alpha - \beta \f$ term
			 * and in the second one always goes \f$ \alpha + \beta \f$ one.
			 * Please also note that no assumptions are made with respect to return values' content
			 * (e.g., it is not guaranteed that return values are empty).
			 * @param[in] t2 factor.
			 * @param[out] ret1 first return value.
			 * @param[out] ret2 second return value.
			 */
			void multiply(const TrigVector &trigVector, TrigVector &result1, TrigVector &result2) const
			// NOTE: we are not using here a general version of vector addition/subtraction
			// because this way we can do two operations (+ and -) every cycle. This is a performance
			// critical part, so the optimization should be worth the hassle.
			{
				const size_type maxw = this->size();
                const size_type minw = trigVector.size();
				
                // Assert widths, *this should always come from a regular Poisson series, and its width should hence be
				// already adjusted my mergeArgs in multiplication routines.
				PIRANHA_ASSERT(maxw >= minw);

				// Adjust the width of retvals, if needed.
				result1.resize(maxw);
				result2.resize(maxw);
				PIRANHA_ASSERT(result1.size() == maxw);
				PIRANHA_ASSERT(result2.size() == maxw);
				size_type i;
				// TODO: improve speed here.
				for (i = 0; i < minw; ++i) 
				{
					result1[i] = (*this)[i] - trigVector[i];
					result2[i] = (*this)[i] + trigVector[i];
				}
				for (; i < maxw; ++i) 
				{
					result1[i] = (*this)[i];
					result2[i] = (*this)[i];
				}
			}


			/// Get flavour.
			bool getFlavour() const
			{
				return flavour;
			}


			/// Set flavour.
            // true: cos; false: sin
			void setFlavour(bool const f)
			{
				flavour = f;
			}


			// I/O.
			template <class ArgsTuple>
			void printPlain(std::ostream &outStream, const ArgsTuple &argsTuple) const
			{
				PIRANHA_ASSERT(argsTuple.template get<Ancestor::position>().size() == this->size());
				(void)argsTuple;

				this->printElements(outStream);
				// Print the separator before flavour only if we actually printed something above.
				if (this->size() != 0) 
				{
					outStream << this->separator;
				}

				if (flavour) 
				{
					outStream << 'c';

				} else 
                {
					outStream << 's';
				}
			}


			template <class ArgsTuple>
			void printPretty(std::ostream &outStream, const ArgsTuple &argsTuple) const
			{
				PIRANHA_ASSERT(argsTuple.template get<Ancestor::position>().size() == this->size());

				if (flavour) 
				{
					outStream << "cos(";

				} else 
				{
					outStream << "sin(";
				}

				bool printedSomething = false;
				for (size_type i = 0; i < this->size(); ++i) 
				{
					const value_type &n = (*this)[i];
					// Don't print anything if n is zero.
					if (n != 0) 
					{
						// If we already printed something and n is positive we are going to print the sign too.
						if (printedSomething && n > 0) 
						{
							outStream << '+';
						}
						// Take care of printing the multiplier.
						if (n == 1) 
						{
							;
						} else if (n == -1) 
						{
							outStream << '-';
						} else
						{
							outStream << n << '*';
						}
						outStream << argsTuple.template get<Ancestor::position>()[i].getName();
						printedSomething = true;
					}
				}

				outStream << ')';
			}


			template <class ArgsTuple>
			void printTex(std::ostream &outStream, const ArgsTuple &argsTuple) const
			{
				PIRANHA_ASSERT(argsTuple.template get<Ancestor::position>().size() == this->size());
				if (flavour) 
				{
					outStream << "\\cos\\left(";

				} else 
				{
					outStream << "\\sin\\left(";
				}

				bool printedSomething = false;
				for (size_type i = 0; i < this->size(); ++i)
				{
					const value_type &n = (*this)[i];
					// Don't print anything if n is zero.
					if (n != 0) 
					{
						// If we already printed something and n is positive we are going to print the sign too.
						if (printedSomething && n > 0) 
						{
							outStream << '+';
						}
						// Take care of printing the multiplier.
						if (n == 1) 
						{
							;
						} else if (n == -1) 
						{
							outStream << '-';
						} else 
						{
							trigVectorPrintElementTEX(outStream,n);
						}
						outStream << argsTuple.template get<Ancestor::position>()[i].getName();
						printedSomething = true;
					}
				}
				outStream << "\\right)";
			}


			bool isUnity() const
			{
                //vlavour == true for cos
				return this->elementsAreZero() && flavour;
			}


			/// Total harmonic degree.
			/**
			 * The total harmonic degree is defined as the summation of the values of the trigonometric multipliers.
			 */
			HarmonicDegreeType harmonicDegree() const
			{
				return std::accumulate(this->begin(), this->end(), HarmonicDegreeType(0));
			}


			/// Harmonic degree of the variables at specified positions pos.
			/**
			 * pos_tuple must be a tuple of vectors of (bool,std::size_t) pairs.
			 */
			template <class PosTuple>
			HarmonicDegreeType partialHarmonicDegree(const PosTuple &posTuple) const
			{
				const std::vector<std::pair<bool, std::size_t> > &pos = posTuple.template get<Ancestor::position>();
				const size_type w = this->size();
                const size_type posSize = boost::numeric_cast<size_type>(pos.size());
				HarmonicDegreeType retval(0);
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


			/// Minimum total harmonic degree.
			/**
			 * Provided for use within the harmonic series toolbox, and defined to be equivalent to harmonicDegree().
			 */
			HarmonicDegreeType harmonicOrder() const
			{
				return harmonicDegree();
			}


			/// Minimum total harmonic degree of the variables at specified positions pos.
			/**
			 * Provided for use within the harmonic series toolbox, and defined to be equivalent to partialHarmonicDegree().
			 */
			template <class PosTuple>
			HarmonicDegreeType partialHarmonicOrder(const PosTuple &posSize) const
			{
				return partialHarmonicDegree(posSize);
			}


			/// Norm.
			/**
			 * The norm of a trigonometric part is always one.
			 */
			template <class ArgsTuple>
			double norm(const ArgsTuple &argsTuple) const
			{
				PIRANHA_ASSERT(argsTuple.template get<Ancestor::position>().size() >= this->size());
				(void)argsTuple;
				return 1.;
			}


			/// Time evaluation of arguments.
			/**
			 * Returns the value assumed by the linear combination of arguments at time t.
			 * @param[in] t double time of the evaluation.
			 * @param[in] v vector of piranha::Psym pointers.
			 */
			template <class ArgsTuple>
			double eval(const double t, const ArgsTuple &argsTuple) const
			{
				const size_type w = this->size();
				PIRANHA_ASSERT(w <= argsTuple.template get<Ancestor::position>().size());

				double retval = 0.;
				for (size_type i = 0; i < w; ++i) 
				{
					if ((*this)[i] != 0) 
					{
						retval += trigVectorEvalElement((*this)[i]) * argsTuple.template get<Ancestor::position>()[i].eval(t);
					}
				}

				if (flavour) 
				{
					return std::cos(retval);

				} else 
				{
					return std::sin(retval);
				}
			}


			/// Sign.
			/**
			 * Return the sign of the first non-zero element of the combination. Zero is considered positive.
			 * This function is used to test for canonical form in piranha::poisson_series_term.
			 */
			short int sign() const
			{
				const size_type w = this->size();
				for (size_type i = 0; i < w; ++i) 
				{
					if ((*this)[i] > 0) 
					{
						return 1;
					}
					if ((*this)[i] < 0) 
					{
						return -1;
					}
				}
				return 1;
			}


			// Re-implement swap and trim to take into account the flavour.
			void swap(TrigVector &trigVector)
			{
				Ancestor::swap(trigVector);
				std::swap(flavour, trigVector.flavour);
			}


			template <class TrimFlags, class ArgsTuple>
			TrigVector trim(const TrimFlags &trimFlags, const ArgsTuple &argsTuple) const
			{
				TrigVector retval(Ancestor::trim(trimFlags, argsTuple));
				retval.flavour = flavour;
				return retval;
			}


			/// All multipliers are zero and flavour is sine.
			template <class ArgsTuple>
			bool isIgnorable(const ArgsTuple &) const
			{
				return !flavour && this->elementsAreZero();
			}


			/// Equality test.
			bool operator==(const TrigVector &trigVector) const
			{
				return (flavour == trigVector.flavour && this->elementsEqualTo(trigVector));
			}


			/// Less than.
			bool operator<(const TrigVector &trigVector) const
			{
                //sin before cos
				if (flavour < trigVector.flavour)
                {
					return true;

				} else if (flavour > trigVector.flavour)
                {
					return false;
				}

                // both are either sin or cos
				return this->lexComparison(trigVector);
			}


			/// Calculate hash_value.
			/**
			 * Used by the hash_value overload for piranha::BaseTerm.
			 */
			std::size_t hash_value() const
			{
				std::size_t retval = this->elementsHasher();
				boost::hash_combine(retval, flavour);
				return retval;
			}


			/// Partial derivative.
			template <class Series, class PosTuple, class ArgsTuple>
			Series partial(const PosTuple &posTuple, const ArgsTuple &argsTuple) const
			{
				Series retval;
				// Do something only if the argument of the partial derivation is present in the trigonometric vector.
				// Otherwise the above retval will return, and it will deliver a zero integer multiplier to be
				// multiplied by the coefficient in the partial derivation of the whole term.
				PIRANHA_ASSERT(posTuple.template get<Ancestor::position>().size() == 1);

				if (posTuple.template get<Ancestor::position>()[0].first) 
				{
					TrigVector copy(*this);
					const size_type pos = boost::numeric_cast<size_type>(posTuple.template get<Ancestor::position>()[0].second);
					// Change the flavour of the resulting key.
					copy.flavour = !flavour;

					PIRANHA_ASSERT(pos < this->size());
					
                    retval = Series::baseSeriesFromKey(copy, argsTuple);
					if (flavour) 
					{
						retval.baseMultBy(-1, argsTuple);
					}
					retval.baseMultBy((*this)[pos], argsTuple);
				}
				return retval;
			}


			/// Exponentiation.
			template <class ArgsTuple>
			TrigVector pow(const double y, const ArgsTuple &) const
			{
				return powNumber(y);
			}


			template <class ArgsTuple>
			TrigVector pow(const mp_rational &q, const ArgsTuple &) const
			{
				return powNumber(q);
			}


			// NOTE: here argsTuple must be the merge of the series undergoing the substitution and
			// the series used for the substitution.
			template <class ResultSeries, class PosTuple, class SubCaches, class ArgsTuple>
			ResultSeries sub(const PosTuple &posTuple, SubCaches &subCaches, const ArgsTuple &argsTuple) const
			{
				ResultSeries retval;
				// If the argument is not present here, the return series will have one term consisting
				// of a unitary coefficient and this very TrigVector.
				// NOTE: for now we can substitute one symbol at a time.
				PIRANHA_ASSERT(posTuple.template get<Ancestor::position>().size() == 1);

				if (!posTuple.template get<Ancestor::position>()[0].first) 
				{
					retval = ResultSeries::baseSeriesFromKey(*this, argsTuple);

				} else 
				{
					const size_type pos = boost::numeric_cast<size_type>(posTuple.template get<Ancestor::position>()[0].second);
					const value_type &power = (*this)[pos];

					PIRANHA_ASSERT(pos < this->size());

					TrigVector tmpTa(*this);
					// Let's turn off the multiplier associated to the symbol we are substituting.
					tmpTa[pos] = 0;
					// NOTE: important: we need key builders here because we may be building RetSeries
					// whose key is _not_ a TrigVector, in principle, so we cannot build a term consisting
					// of a TrigVector and unity coefficient and simply insert it.
					// Build the orig_cos series.
					tmpTa.setFlavour(true);
					ResultSeries originalCos = ResultSeries::baseSeriesFromKey(tmpTa, argsTuple);
					// Build the orig_sin series.
					tmpTa.setFlavour(false);
					ResultSeries originalSin = ResultSeries::baseSeriesFromKey(tmpTa, argsTuple);

					PIRANHA_ASSERT(retval.empty());
					
                    if (this->getFlavour()) 
					{
						retval.baseAdd(originalCos, argsTuple);
						retval.baseMultBy(subCaches.template get<Ancestor::position>()[power].baseReal(argsTuple), argsTuple);
						originalSin.baseMultBy(subCaches.template get<Ancestor::position>()[power].baseImag(argsTuple), argsTuple);
						retval.baseSubtract(originalSin, argsTuple);

					} else 
					{
						retval.baseAdd(originalSin, argsTuple);
						retval.baseMultBy(subCaches.template get<Ancestor::position>()[power].baseReal(argsTuple), argsTuple);
						originalCos.baseMultBy(subCaches.template get<Ancestor::position>()[power].baseImag(argsTuple), argsTuple);
						// NOTE: series multadd here (and multiply by -1 to do subtraction too)?
						// Below too...
						retval.baseAdd(originalCos, argsTuple);
					}
				}
				return retval;
			}


			template <class ResultSeries, class PosTuple, class SubCaches, class ArgsTuple>
			ResultSeries eiSub(const PosTuple &pos_tuple, SubCaches &subCaches, const ArgsTuple &argsTuple) const
			{
				ResultSeries retval;
				
                PIRANHA_ASSERT(pos_tuple.template get<Ancestor::position>().size() == 1);

				if (!pos_tuple.template get<Ancestor::position>()[0].first) 
				{
					retval = ResultSeries::baseSeriesFromKey(*this, argsTuple);

				} else 
				{
					const size_type pos = boost::numeric_cast<size_type>(pos_tuple.template get<Ancestor::position>()[0].second);
					const value_type &power = (*this)[pos];

					PIRANHA_ASSERT(pos < this->size());
					
                    TrigVector tmpTa(*this);
					tmpTa[pos] = 0;
					tmpTa.setFlavour(true);
					ResultSeries originalCos = ResultSeries::baseSeriesFromKey(tmpTa, argsTuple);
					tmpTa.setFlavour(false);
					ResultSeries originalSin = ResultSeries::baseSeriesFromKey(tmpTa, argsTuple);

					PIRANHA_ASSERT(retval.empty());
					
                    if (flavour) 
					{
						retval.baseAdd(originalCos, argsTuple);
						retval.baseMultBy(subCaches.template get<Ancestor::position>()[power].baseReal(argsTuple), argsTuple);
						originalSin.baseMultBy(subCaches.template get<Ancestor::position>()[power].baseImag(argsTuple), argsTuple);
						retval.baseSubtract(originalSin, argsTuple);
					} else 
					{
						retval.baseAdd(originalSin, argsTuple);
						retval.baseMultBy(subCaches.template get<Ancestor::position>()[power].baseReal(argsTuple), argsTuple);
						originalCos.baseMultBy(subCaches.template get<Ancestor::position>()[power].baseImag(argsTuple), argsTuple);
						retval.baseAdd(originalCos, argsTuple);
					}
				}

				return retval;
			}

		private:
			template <class Number>
			TrigVector powNumber(const Number &y) const 
			{
				const bool intZero = this->elementsAreZero();
				TrigVector retval;
				if (y < 0) 
				{
					if (intZero && !flavour) 
					{
						// 0**-y.
						PIRANHA_THROW(zero_division_error, "cannot divide by zero");
					} else if (intZero && flavour) 
					{
						// 1**-y == 1. Don't do anything because retval is already initialized properly.
						;
					} else 
					{
						// x**-y -> no go.
						PIRANHA_THROW(value_error, "non-unity trigonometric vector is not suitable for negative exponentiation");
					}
				} else if (y == 0) 
				{
					// x**0 == 1. Don't do nothing because retval is already initialized properly.
					;
				} else 
				{
					if (intZero && !flavour) 
					{
						// 0**y == 0.
						retval.flavour = false;
					} else if (intZero && flavour) 
					{
						// 1**y == 1. Don't do anything because retval is already initialized properly.
						;
					} else 
					{
						// x**y --> no go.
						PIRANHA_THROW(value_error, "non-unity trigonometric vector is not suitable for positive exponentiation");
					}
				}
				return retval;
			}

		private:
			//true: cos; false: sin  
			bool flavour;
	};


	/// is_ring_exact type trait specialisation for TrigVector.
	template <__PIRANHA_TRIG_VECTOR_TP_DECL>
	struct is_ring_exact<TrigVector<__PIRANHA_TRIG_VECTOR_TP> >: boost::true_type {};


	/// is_trig_exact type trait specialisation for TrigVector.
	template <__PIRANHA_TRIG_VECTOR_TP_DECL>
	struct is_trig_exact<TrigVector<__PIRANHA_TRIG_VECTOR_TP> >: boost::true_type {};
}

#undef __PIRANHA_TRIG_VECTOR_TP_DECL
#undef __PIRANHA_TRIG_VECTOR_TP

#endif
