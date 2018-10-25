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

#ifndef PIRANHA_BASE_SERIES_MULTIPLIER_H
#define PIRANHA_BASE_SERIES_MULTIPLIER_H

#include <algorithm>
#include <boost/lambda/lambda.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/thread/thread.hpp>
#include <boost/type_traits/is_same.hpp> // For key type detection.
#include <boost/tuple/tuple.hpp>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <vector>

#include "../config.h"
#include "../exceptions.h"
#include "../settings.h"
#include "../stats.h"
#include "base_series_multiplier_mp.h"
#include "null_truncator.h"


#define DC static_cast<Derived *>(this)

namespace piranha
{
	/// Base series multiplier.
	/**
	 * This class is meant to be extended to build specific multipliers.
	 */
	template <class Series1, class Series2, class ArgsTuple, class Truncator, class Derived>
	class BaseSeriesMultiplier
	{
	//		friend class base_insert_multiplication_result;

		protected:
			// Alias for term type of first input series and return value series.
			typedef typename Series1::TermType TermType1;
			// Alias for term type of second input series.
			typedef typename Series2::TermType TermType2;

            static_assert((boost::is_same<typename TermType1::KeyType, typename TermType2::KeyType>::value), "Key type mismatch in base multiplier.");

			/// Compute block size for multiplication.
			/**
			 * Resulting block size depends on the cache memory available, and will always be in the [16,8192] range. N is
			 * the size of the base storage unit type used to store the result of the multiplication (e.g., coefficient type in vector coded, coded hash table
			 * term in the sparse hashed multiplication, etc.).
			 */
			template <std::size_t N>
			static std::size_t computeBlockSize()
			{
				// NOTE: this function is used typically considering only the output storage requirements, since storage of input series
				//       will be a small fraction of storage for the output series.
                static_assert(N > 0, "");

				const std::size_t shift = boost::numeric_cast<std::size_t>(
					        std::log(std::max<double>(16., std::sqrt(static_cast<double>((settings::cache_size * 1024) / N)))) / std::log(2.) - 1
				);

				return (std::size_t(2) << std::min<std::size_t>(std::size_t(12), shift));
			}


			template <class MulitplicationFunctor>
			static void blockedMultiplication(std::size_t const blockSize, std::size_t const size1, std::size_t const size2, MulitplicationFunctor &m)
			{
				PIRANHA_ASSERT(blockSize > 0);

				const std::size_t nblocks1 = size1 / blockSize;
                const std::size_t nblocks2 = size2 / blockSize;

				for (std::size_t n1 = 0; n1 < nblocks1; ++n1)
                {
					const std::size_t iStart = n1 * blockSize; 
                    const std::size_t iEnd   = iStart + blockSize;

					// regulars1 * regulars2
					for (std::size_t n2 = 0; n2 < nblocks2; ++n2)
                    {
						const std::size_t jStart = n2 * blockSize; 
                        const std::size_t jEnd   = jStart + blockSize;

						for (std::size_t i = iStart; i < iEnd; ++i)
                        {
							for (std::size_t j = jStart; j < jEnd; ++j)
                            {
								if (!m(i, j))
                                {
									break;
								}
							}
						}
					}

					// regulars1 * rem2
					for (std::size_t i = iStart; i < iEnd; ++i)
                    {
						for (std::size_t j = nblocks2 * blockSize; j < size2; ++j)
                        {
							if (!m(i, j))
                            {
								break;
							}
						}
					}
				}

				// rem1 * regulars2
				for (std::size_t n2 = 0; n2 < nblocks2; ++n2)
                {
					const std::size_t jStart = n2 * blockSize; 
                    const std::size_t jEnd   = jStart + blockSize;

					for (std::size_t i = nblocks1 * blockSize; i < size1; ++i)
                    {
						for (std::size_t j = jStart; j < jEnd; ++j)
                        {
							if (!m(i, j))
                            {
								break;
							}
						}
					}
				}

				// rem1 * rem2.
				for (std::size_t i = nblocks1 * blockSize; i < size1; ++i)
                {
					for (std::size_t j = nblocks2 * blockSize; j < size2; ++j)
                    {
						if (!m(i, j))
                        {
							break;
						}
					}
				}
			}


			/// Cache pointers to series' terms in the internal storage.
			template <class Container1, class Container2>
			void cacheTermsPointers(const Container1 &container1, const Container2 &container2)
			{
				PIRANHA_ASSERT(terms1.empty() && terms2.empty());

				std::transform(container1.begin(), container1.end(), std::insert_iterator<std::vector<typename Series1::TermType const *> >(terms1, terms1.begin()), &(boost::lambda::_1));
				std::transform(container2.begin(), container2.end(), std::insert_iterator<std::vector<typename Series2::TermType const *> >(terms2, terms2.begin()), &(boost::lambda::_1));
			}


		public:

			BaseSeriesMultiplier(const Series1 &series1, const Series2 &series2, Series1 &result, const ArgsTuple &argsTuple)
              :	series1(series1), series2(series2), argsTuple(argsTuple), retval(result)
			{
				PIRANHA_ASSERT(series1.length() > 0 && series2.length() > 0);
			}


			// Plain multiplication.
			void performPlainMultiplication()
			{
				performPlainThreadedMultiplication();
			}

		private:

			template <class GenericTruncator>
			struct PlainFunctor
            {
				typedef typename TermType1::multiplication_result ResultType;

				PlainFunctor(ResultType &result, const TermType1 **term1, const TermType2 **term2, const GenericTruncator &truncator, Series1 &retval, const ArgsTuple &argsTuple)
                             : result(result), term1(term1), term2(term2), truncator(truncator), retval(retval), argsTuple(argsTuple)
				{}

				bool operator()(const std::size_t i, const std::size_t j)
				{
					if (truncator.skip(&term1[i], &term2[j]))
                    {
						return false;
					}

					TermType1::multiply(*term1[i], *term2[j], result, argsTuple);
					InsertMultiplicationResult<ResultType>::run(result, retval, argsTuple);

					return true;
				}


				ResultType		        &result;
				const TermType1	        **term1;
				const TermType2	        **term2;
				const GenericTruncator	&truncator;
				Series1			        &retval;
				const ArgsTuple		    &argsTuple;
			};



			struct PlainWorker {

				PlainWorker(BaseSeriesMultiplier &multiplier, Series1 &retval)
                    : multiplier(multiplier), retval(retval), terms1(multiplier.terms1)
				{}

				PlainWorker(BaseSeriesMultiplier &multiplier, Series1 &retval, std::vector<std::vector<TermType1 const *> > &split1, const std::size_t idx)
                    : multiplier(multiplier), retval(retval), terms1(split1[idx])
				{}

				void operator()()
				{
					// Build the truncator.
					const typename Truncator::template GetType<Series1, Series2, ArgsTuple> truncator(terms1, multiplier.terms2, multiplier.argsTuple);
					// Use the selected truncator only if it really truncates, otherwise use the
					// null truncator.
					if (truncator.isEffective())
                    {
						plainImplementation(truncator);

					} else
                    {
						plainImplementation(NullTruncator::template GetType<Series1, Series2, ArgsTuple>(terms1, multiplier.terms2, multiplier.argsTuple) );
					}
				}


				template <class GenericTruncator>
				void plainImplementation(const GenericTruncator &truncator)
				{
					typedef typename TermType1::multiplication_result ResultType;
					ResultType res;
					const std::size_t size1 = terms1.size();
                    const std::size_t size2 = multiplier.terms2.size();

					PIRANHA_ASSERT(size1 && size2);

					const TermType1 **t1 = &terms1[0];
					const TermType2 **t2 = &multiplier.terms2[0];
					PlainFunctor<GenericTruncator> plainFunctor(res, t1, t2, truncator, retval, multiplier.argsTuple);
					const std::size_t blockSize = computeBlockSize<boost::tuples::length<ResultType>::value * sizeof(TermType1)>();

					blockedMultiplication(blockSize, size1, size2, plainFunctor);
				}


				BaseSeriesMultiplier		    &multiplier;
				Series1				            &retval;
				std::vector<TermType1 const *>	terms1;
			};


			// Threaded multiplication.
			void performPlainThreadedMultiplication()
			{
				// Effective number of threads to use. If the two series are small, we want to use one single thread.
				// NOTE: here the number 2500 is a kind of rule-of thumb. Basically multiplications of series < 50 elements
				// will use just one thread.
				if (double(terms1.size()) * double(terms2.size()) < 2500)
                {
					stats::trace_stat("mult_st", std::size_t(0), increment);

					PlainWorker w(*DC, retval);
					w();

				} else
                {
					stats::trace_stat("mult_mt", std::size_t(0), increment);
					// TODO: fix numeric casting here.
					// If size1 is less than the number of desired threads,
					// use size1 as number of threads.
					const std::size_t n = std::min(boost::numeric_cast<typename std::vector<TermType1 const *>::size_type > (settings::get_nthread()), terms1.size());
					std::vector<std::vector<TermType1 const *> > split1(n);
					// m is the number of terms per thread for regular blocks.
					const std::size_t m = terms1.size() / n;

					// Iterate up to n - 1 because that's the number up to which we can divide series1 into
					// regular blocks.
					for (std::size_t i = 0; i < n - 1; ++i)
                    {
						split1[i].insert(split1[i].end(), terms1.begin() + i * m, terms1.begin() + (i + 1) * m);
					}


					// Last iteration.
					split1[n - 1].insert(split1[n - 1].end(), terms1.begin() + (n - 1) * m, terms1.end());
					boost::thread_group threadGroup;
					std::vector<Series1> retvals(n, Series1());
					for (std::size_t i = 0; i < n; ++i)
                    {
						threadGroup.create_thread(PlainWorker(*DC, retvals[i], split1, i));
					}
					threadGroup.join_all();

					// Take the retvals and insert them into final retval.
					for (std::size_t i = 0; i < n; ++i)
                    {
						retval.insertRange(retvals[i].begin(), retvals[i].end(), argsTuple);
					}
				}
			}


		public:
			// TODO: make these protected?
			// References to the series.
			const Series1					&series1;
			const Series2					&series2;
			// Reference to the arguments tuple.
			const ArgsTuple					&argsTuple;
			// Reference to the result.
			Series1						    &retval;
			// Vectors of pointers to the input terms.
			std::vector<TermType1 const *>	terms1;
			std::vector<TermType2 const *>	terms2;
	};
}

#undef DC

#endif
