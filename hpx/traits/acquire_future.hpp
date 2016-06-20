//  Copyright (c) 2007-2016 Hartmut Kaiser
//  Copyright (c) 2016 Agustin Berge
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef HPX_TRAITS_ACQUIRE_FUTURE_HPP
#define HPX_TRAITS_ACQUIRE_FUTURE_HPP

#include <hpx/config.hpp>
#include <hpx/util/range.hpp>
#include <hpx/traits/detail/reserve.hpp>
#include <hpx/traits/is_future.hpp>
#include <hpx/traits/is_future_range.hpp>
#include <hpx/traits/is_range.hpp>

#include <algorithm>
#include <cstddef>
#include <cstddef>
#include <iterator>
#include <type_traits>
#include <utility>

namespace hpx { namespace traits
{
    namespace detail
    {
        template <typename T, typename Enable = void>
        struct acquire_future_impl;
    }

    template <typename T, typename Enable = void>
    struct acquire_future
      : detail::acquire_future_impl<typename std::decay<T>::type>
    {};

    struct acquire_future_disp
    {
        template <typename T>
        HPX_FORCEINLINE typename acquire_future<T>::type
        operator()(T && t) const
        {
            return acquire_future<T>()(std::forward<T>(t));
        }
    };

    namespace detail
    {
        ///////////////////////////////////////////////////////////////////////
        template <typename T, typename Enable>
        struct acquire_future_impl
        {
            static_assert(!is_future_or_future_range<T>::value, "");

            typedef T type;

            template <typename T_>
            HPX_FORCEINLINE
            T operator()(T_ && value) const
            {
                return std::forward<T_>(value);
            }
        };

        ///////////////////////////////////////////////////////////////////////
        template <typename R>
        struct acquire_future_impl<hpx::lcos::future<R> >
        {
            typedef hpx::lcos::future<R> type;

            HPX_FORCEINLINE hpx::lcos::future<R>
            operator()(hpx::lcos::future<R>& future) const
            {
                return std::move(future);
            }

            HPX_FORCEINLINE hpx::lcos::future<R>
            operator()(hpx::lcos::future<R>&& future) const
            {
                return std::move(future);
            }
        };

        template <typename R>
        struct acquire_future_impl<hpx::lcos::shared_future<R> >
        {
            typedef hpx::lcos::shared_future<R> type;

            HPX_FORCEINLINE hpx::lcos::shared_future<R>
            operator()(hpx::lcos::shared_future<R> const& future) const
            {
                return future;
            }

            HPX_FORCEINLINE hpx::lcos::shared_future<R>
            operator()(hpx::lcos::shared_future<R>&& future) const
            {
                return std::move(future);
            }
        };

        ///////////////////////////////////////////////////////////////////////
        template <typename Range>
        struct acquire_future_impl<
            Range,
            typename std::enable_if<
                traits::is_future_range<Range>::value
            >::type
        >
        {
            typedef typename traits::future_range_traits<Range>::future_type
                future_type;

            typedef Range type;

            template <typename Range_>
            HPX_FORCEINLINE Range
            operator()(Range_&& futures) const
            {
                Range values;
                detail::reserve_if_random_access(values, futures);

                std::transform(
                    util::begin(futures), util::end(futures),
                    std::back_inserter(values), acquire_future_disp());

                return values;
            }
        };
    }
}}

#endif /*HPX_TRAITS_ACQUIRE_FUTURE_HPP*/
