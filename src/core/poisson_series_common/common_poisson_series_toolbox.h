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

#ifndef PIRANHA_COMMON_POISSON_SERIES_TOOLBOX_H
#define PIRANHA_COMMON_POISSON_SERIES_TOOLBOX_H

#include <algorithm> // For sorting.
#include <boost/lambda/lambda.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/type_traits/is_same.hpp>
#include <complex>
#include <cstddef>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

#include "../config.h"
#include "../exceptions.h"
#include "../math.h"
#include "../mp.h"
#include "../ntuple.h"
#include "../Psym.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	// NOTE: this seems to be valid only for specific position of arguments. Document. Maybe put static const int "trig_index" or something like
	// that, as we do for power series?
	template <class Derived>
	class common_poisson_series
	{
		public:
			std::complex<Derived> ei() const
			{
				// In order to account for a potential integer linear combination of arguments
				// we must merge in as trigonometric arguments the polynomial arguments. The safe
				// way to do this is by using NamedSeries::mergeArgs with a phony series having zero
				// polynomial arguments and as trigonometric arguments the polynomial arguments of this.
				Derived copy(*derived_const_cast), tmp;
				typename Derived::ArgsTupleType tmp_args;
				tmp_args.template get<1>() = derived_const_cast->arguments().template get<0>();
				tmp.setArguments(tmp_args);
				copy.mergeArgs(tmp);
				// Now we can call in the ei method from above.
				std::complex<Derived> retval(copy.base_ei(copy.arguments()));
				retval.setArguments(copy.arguments());
				retval.trim();
				return retval;
			}

			Derived cos() const
			{
				return ei().real();
			}

			Derived sin() const
			{
				return ei().imag();
			}

			// We have to specialise this in order to prepare the arguments tuple for the fact
			// that poly args of s may be added as trig args (as s will be used as argument for sines
			// and/or cosines).
			template <class SubSeries>
			Derived sub(const std::string &name, const SubSeries &series) const
			{
				typedef typename Derived::ArgsTupleType ArgsTupleType;
				typedef typename NTuple<std::vector<std::pair<bool, std::size_t> >, Derived::echelonLevel + 1>::Type pos_tuple_type;
				typedef typename Derived::TermType::CfType::
					template SubstitutionCacheSelector<SubSeries, typename Derived::TermType::KeyType::
					template SubstitutionCacheSelector<SubSeries, boost::tuples::null_type, ArgsTupleType>
					::Type, ArgsTupleType>::Type sub_caches_type;

				PIRANHA_STATIC_CHECK(boost::tuples::length<sub_caches_type>::value == boost::tuples::length<pos_tuple_type>::value,
					                 "Size mismatch for position and cache tuples in Poisson series substitution.");
				const Psym p(name);
				Derived this_copy(*derived_const_cast);
				SubSeries s_copy(series), tmp;
				// Assign as tmp's trig arguments series's polynomial arguments.
				ArgsTupleType tmp_args;
				tmp_args.template get<1>() = series.arguments().template get<0>();
				tmp.setArguments(tmp_args);
				// After the next line, s_copy's args layout is compatible with tmp's.
				s_copy.mergeArgs(tmp);
				// After the next line, this_copy's args layout is compatible with s_copy's
				this_copy.mergeArgs(s_copy);
				// Finally, have s_copy have compatible arguments with this_copy. This is needed because
				// we will be using this_copy's arguments as argsTuple in all base functions used from now
				// on, including functions taking s_copy as arguments and which do not know anything about
				// this_copy.
				s_copy.mergeArgs(this_copy);
				// Init sub caches using s and this_copy.m_arguments.
				sub_caches_type sub_caches;
				InitSubstitutionCaches<sub_caches_type, SubSeries, ArgsTupleType>::run(sub_caches, s_copy, &this_copy.arguments());
				const pos_tuple_type pos_tuple = psyms2pos(std::vector<Psym>(1, p), this_copy.arguments());

				Derived retval(this_copy.template baseSub<Derived, typename Derived::SubstitutionFunctor>(pos_tuple, sub_caches, this_copy.arguments()));
				retval.setArguments(this_copy.arguments());
				retval.trim();

				return retval;
			}


			template <class SubSeries>
			Derived eiSubstitute(const std::string &name, const SubSeries &series) const
			{
				typedef typename Derived::ArgsTupleType ArgsTupleType;
				typedef typename NTuple<std::vector<std::pair<bool, std::size_t> >, Derived::echelonLevel + 1>::Type pos_tuple_type;
				typedef typename Derived::TermType::CfType::
					template EiSubstitutionCacheSelector<SubSeries, typename Derived::TermType::KeyType::
					template EiSubstitutionCacheSelector<SubSeries, boost::tuples::null_type, ArgsTupleType>
					::Type, ArgsTupleType>::Type sub_caches_type;

				PIRANHA_STATIC_CHECK(boost::tuples::length<sub_caches_type>::value == boost::tuples::length<pos_tuple_type>::value,
					                 "Size mismatch for position and cache tuples in Poisson series ei substitution.");
				const Psym p(name);
				Derived this_copy(*derived_const_cast);
				SubSeries s_copy(series);
				this_copy.mergeArgs(s_copy);
				s_copy.mergeArgs(this_copy);
				sub_caches_type sub_caches;
				InitSubstitutionCaches<sub_caches_type, SubSeries, ArgsTupleType>::run(sub_caches, s_copy, &this_copy.arguments());

				const pos_tuple_type pos_tuple = psyms2pos(std::vector<Psym>(1,p), this_copy.arguments());

				Derived retval(this_copy.template baseSub<Derived, EiSubstitutionFunctor>(pos_tuple, sub_caches, this_copy.arguments()));
				retval.setArguments(this_copy.arguments());
				retval.trim();
				
                return retval;
			}


			// NOTE: either this must be generalised somehow, or dropped completely.
			template <class FourierSeries>
			FourierSeries to_fs() const
			{
				typedef typename Derived::const_iterator const_iterator;
				typedef typename FourierSeries::TermType fourier_term;
				typename NTuple<VectorPsym,1>::Type argsTuple(derived_const_cast->arguments().template get<1>());
				FourierSeries retval;
				retval.setArguments(argsTuple);
				const const_iterator it_f = derived_const_cast->end();

				for (const_iterator it = derived_const_cast->begin(); it != it_f; ++it)
                {
					if (!it->cf.isSingleCf())
                    {
						PIRANHA_THROW(value_error,"polynomial coefficient cannot be converted to numerical coefficient");
					}

					retval.insert(fourier_term(typename fourier_term::CfType(it->cf.begin()->cf), typename fourier_term::KeyType(it->key)), argsTuple);
				}

				return retval;
			}


			Derived integrate(const std::string &name) const
			{
				typedef typename NTuple<std::vector<std::pair<bool, std::size_t> >, 2>::Type pos_tuple_type;
				const Psym p(name);
				const pos_tuple_type pos_tuple = psyms2pos(VectorPsym(1,p), derived_const_cast->arguments());
				Derived retval;

				if (pos_tuple.template get<0>()[0].first || pos_tuple.template get<1>()[0].first) 
                {
					// If the symbol is present either as a poly arg or a trig arg, invoke the base integration routine.
					// We must prepare arguments list to take an extra polynomial argument, for cases such as 1 + cos(x)
					// where integration adds a symbol to the list of polynomial symbols.
					Derived this_copy(*derived_const_cast);
					this_copy.mergeArgs(Derived(p));
					const pos_tuple_type new_pos_tuple = psyms2pos(VectorPsym(1,p),this_copy.arguments());

					retval = this_copy.baseIntegrate(new_pos_tuple, this_copy.arguments());
					retval.setArguments(this_copy.arguments());
					retval.trim();
				} else 
                {
					// If the symbol is not present at all in the series, just multiply this by the series generated
					// by the input symbol.
					retval = *derived_const_cast;
					retval *= Derived(p);
				}

				return retval;
			}


		//protected:
			// Integrate supposing that the symbol is present in the Poisson series.
			template <typename PosTuple, typename ArgsTuple>
			Derived baseIntegrate(const PosTuple &pos_tuple, const ArgsTuple &argsTuple) const
			{
				PIRANHA_STATIC_CHECK(boost::tuples::length<PosTuple>::value == boost::tuples::length<ArgsTuple>::value,
					"Size mismatch between args tuple and pos tuple in Poisson series integration.");

				typedef typename Derived::const_iterator                     const_iterator;
				typedef typename Derived::TermType::CfType::DegreeType       degree_type;
				typedef typename Derived::TermType::KeyType::HarmonicDegreeType HarmonicDegreeType;

				// Make sure that the position tuple contains just one symbol in each element of the tuple,
				// and that the symbol is present in the series.
				PIRANHA_ASSERT(pos_tuple.template get<0>().size() == 1 && pos_tuple.template get<1>().size() == 1);
				PIRANHA_ASSERT(pos_tuple.template get<0>()[0].first || pos_tuple.template get<1>()[0].first);
				
                Derived retval;
				const const_iterator it_f = derived_const_cast->end();
				for (const_iterator it = derived_const_cast->begin(); it != it_f; ++it) 
                {
					if (pos_tuple.template get<1>()[0].first && it->key[pos_tuple.template get<1>()[0].second] != 0) 
                    {
						// Integrand argument appears as trigonometric argument: try to integrate recursively by parts.
						typedef typename Derived::TermType::CfType::DegreeType degree_type;
						const degree_type degree(it->cf.partialDegree(pos_tuple));
						if (degree < 0 || !is_integer(degree)) 
                        {
							PIRANHA_THROW(value_error,"cannot integrate Poisson series term if the polynomial degree of the integrand argument "
								"is negative or not an integer");
						}

						const HarmonicDegreeType trig_mult(it->key[pos_tuple.template get<1>()[0].second]);
						typename Derived::TermType tmp(*it);
						typename Derived::TermType::CfType tmp_cf(it->cf);

						tmp_cf.divideBy(trig_mult, argsTuple);
						
                        for (mp_integer i(0); i < degree + 1; i += 1) 
                        {
							// Flip the flavour of the trigonometric part.
							tmp.key.setFlavour(!tmp.key.getFlavour());
							const mp_integer rem(i % 4);
							int mult;
							if (rem == 0) 
                            {
								if (it->key.getFlavour()) 
                                {
									mult = 1;
								} else 
                                {
									mult = -1;
								}
							} else if (rem == 1) 
                            {
								mult = -1;
							} else if (rem == 2) 
                            {
								if (it->key.getFlavour()) 
                                {
									mult = -1;
								} else 
                                {
									mult = 1;
								}
							} else 
                            {
								PIRANHA_ASSERT(rem == 3);
								mult = 1;
							}

							tmp.cf = tmp_cf;
							tmp.cf.multBy(mult * cs_phase(i),argsTuple);

							retval.insert(tmp,argsTuple);
							
                            // Prepare tmp's cf for next step.
							tmp_cf = tmp_cf.template partial<typename Derived::TermType::CfType>(pos_tuple, argsTuple);
							tmp_cf.divideBy(trig_mult, argsTuple);
						}
					} else 
                    {
						typename Derived::TermType tmp(*it);
						tmp.cf = tmp.cf.integrate(pos_tuple, argsTuple);

						retval.insert(tmp, argsTuple);
					}
				}

				return retval;
			}


			// NOTICE: move this into private?
			// NOTICE: this method assumes that the input args tuple already hase merged in as
			// trig arguments the poly arguments (see also below).
			template <class ArgsTuple>
			std::complex<Derived> base_ei(const ArgsTuple &argsTuple) const
			{
				typedef typename std::complex<Derived>::TermType                 complex_term_type;
				typedef typename Derived::TermType                               term_type;
				typedef typename term_type::CfType::TermType::CfType           poly_cf_type;
				typedef typename std::complex<Derived>::TermType::CfType        complex_cf_type;
				typedef typename std::vector<term_type const *>::const_iterator  const_iterator;

				// Cache and sort the terms according to the criterion defined in the truncator.
				std::vector<term_type const *> cache;
				try {
					cache = derived_const_cast->template get_sorted_series<Derived>(argsTuple);
					// Reverse the series, we want to start multiplication from the least significant terms.
					std::reverse(cache.begin(), cache.end());
				} catch (const value_error &) {
					// Cache term pointers.
					std::transform(derived_const_cast->begin(),derived_const_cast->end(),
						std::insert_iterator<std::vector<term_type const *> >(cache,cache.begin()),
						&(boost::lambda::_1));
				}

				// Determine the iterator to be avoided - parts of which we might be able to treat exactly.
				const_iterator it_avoid = cache.end();
				for (const_iterator it = cache.begin(); it != cache.end(); ++it) 
                {
					if ((*it)->key.isUnity()) 
                    {
						it_avoid = it;
						break;
					}
				}

				// Expand using Jacobi-Anger's identity.
				std::complex<Derived> retval;
				derived_const_cast->jacang(cache, it_avoid, retval, argsTuple);
				// Now handle the avoided iterator, if any.
				if (it_avoid != cache.end()) 
                {
					// Split the polynomial coefficient in two parts: exactly treatable and not.
					typename term_type::CfType exact;
                    typename term_type::CfType residual;
					for (typename term_type::CfType::const_iterator it = (*it_avoid)->cf.begin(); it != (*it_avoid)->cf.end(); ++it) 
                    {
						// Exact part has the following requisites: exactly one "active" variable with unitary exponent and a coefficient
						// that is convertible to the type representing the harmonic degree.
						typename term_type::CfType::TermType::KeyType::size_type n_active = 0;
						bool has_unitary = false;
						for (typename term_type::CfType::TermType::KeyType::size_type i = 0; i < it->key.size(); ++i) 
                        {
							if (it->key[i] == 1) 
                            {
								has_unitary = true;
							}
							if (it->key[i] != 0) 
                            {
								++n_active;
							}
						}

						bool is_exact = (has_unitary && n_active == 1);
						if (is_exact) 
                        {
							// Try to convert the coefficient.
							try {
								boost::lexical_cast<typename term_type::KeyType::HarmonicDegreeType>(it->cf.get_value());
							} catch (const boost::bad_lexical_cast &) {
								is_exact = false;
							}
						}

						if (is_exact) 
                        {
							exact.insert(*it,argsTuple);
						} else 
                        {
							residual.insert(*it,argsTuple);
						}
					}

					// Let's deal with the exact part.
					if (exact.length()) 
                    {
						std::complex<Derived> tmp_series;
						// Prepare the terms to be inserted.
						complex_term_type tmp_term1;
						tmp_term1.cf = complex_cf_type(std::complex<double>(1, 0), argsTuple);
						tmp_term1.key.resize(boost::numeric_cast<typename complex_term_type::KeyType::size_type>(argsTuple.template get<1>().size()));
						tmp_term1.key.setFlavour(true);
					
                        complex_term_type tmp_term2;
						tmp_term2.cf = complex_cf_type(std::complex<double>(0, 1), argsTuple);
						tmp_term2.key.resize(boost::numeric_cast<typename complex_term_type::KeyType::size_type>(argsTuple.template get<1>().size()));
						tmp_term2.key.setFlavour(false);

						// Copy over the exact part from polynomial into trigonometric.
						for (typename term_type::CfType::const_iterator it = exact.begin(); it != exact.end(); ++it)
                        {
							// NOTE: we use just 1 size type here, but we should be covered by prior arguments merging.
							typedef typename complex_term_type::KeyType::size_type size_type;
							const size_type pos_poly = boost::numeric_cast<size_type>(
								std::distance(it->key.begin(),
								std::find_if(it->key.begin(),it->key.end(),boost::lambda::_1 == 1)
								)
							);

							PIRANHA_ASSERT(pos_poly < argsTuple.template get<0>().size());

							const size_type pos_trig = boost::numeric_cast<size_type>(
								std::distance(
								argsTuple.template get<1>().begin(),
								std::find_if(argsTuple.template get<1>().begin(),
								argsTuple.template get<1>().end(),
								boost::lambda::_1 ==
								argsTuple.template get<0>()[pos_poly])
								)
							);

							PIRANHA_ASSERT(pos_trig < argsTuple.template get<1>().size());

							tmp_term1.key[pos_trig] = boost::lexical_cast<typename complex_term_type::KeyType::value_type>(it->cf.get_value());
							tmp_term2.key[pos_trig] = boost::lexical_cast<typename complex_term_type::KeyType::value_type>(it->cf.get_value());
						}

						tmp_series.insert(tmp_term1, argsTuple);
						tmp_series.insert(tmp_term2, argsTuple);
						retval.baseMultBy(tmp_series, argsTuple);
					}

					// Finally, the residual.
					if (residual.length()) 
                    {
						term_type tmp;
						tmp.cf       = residual;
						tmp.key      = (*it_avoid)->key;
						term_type *ptr = &tmp;

						retval.baseMultBy(derived_const_cast->jacang_term(&ptr, argsTuple), argsTuple);
					}
				}

				return retval;
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
