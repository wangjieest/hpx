//  Copyright (c) 2014 Grant Mercer
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>
#include <hpx/include/parallel_count.hpp>
#include <hpx/util/lightweight_test.hpp>

#include "test_utils.hpp"

////////////////////////////////////////////////////////////////////////////
template <typename ExPolicy, typename IteratorTag>
void test_count_if(ExPolicy const& policy, IteratorTag)
{
    BOOST_STATIC_ASSERT(hpx::parallel::is_execution_policy<ExPolicy>::value);

    typedef std::vector<std::size_t>::iterator base_iterator;
    typedef test::test_iterator<base_iterator, IteratorTag> iterator;

    std::vector<std::size_t> c(100007);
    std::iota(boost::begin(c), boost::begin(c) + 50, 0);
    std::iota(boost::begin(c) + 50, boost::end(c), std::rand() + 50);

    std::size_t num_items = hpx::parallel::count_if(policy,
        iterator(boost::begin(c)), iterator(boost::end(c)),
        [](std::size_t x) { return x < 50; });

    HPX_TEST_EQ(num_items, 50u);
}

template <typename ExPolicy, typename IteratorTag>
void test_count_if_async(ExPolicy const& p, IteratorTag)
{
    typedef std::vector<std::size_t>::iterator base_iterator;
    typedef std::vector<std::size_t>::difference_type diff_type;
    typedef test::test_iterator<base_iterator, IteratorTag> iterator;

    std::vector<std::size_t> c(10007);
    std::iota(boost::begin(c), boost::begin(c) + 50, 0);
    std::iota(boost::begin(c) + 50, boost::end(c), std::rand() + 50);

    hpx::future<diff_type> f =
        hpx::parallel::count_if(p,
            iterator(boost::begin(c)), iterator(boost::end(c)),
            [](std::size_t x) { return x < 50; });

    HPX_TEST_EQ(f.get(), 50);
}

template <typename IteratorTag>
void test_count_if()
{
    using namespace hpx::parallel;
    test_count_if(seq, IteratorTag());
    test_count_if(par, IteratorTag());
    test_count_if(par_vec, IteratorTag());

    test_count_if_async(seq(task), IteratorTag());
    test_count_if_async(par(task), IteratorTag());

    test_count_if(execution_policy(seq), IteratorTag());
    test_count_if(execution_policy(par), IteratorTag());
    test_count_if(execution_policy(par_vec), IteratorTag());

    test_count_if(execution_policy(seq(task)), IteratorTag());
    test_count_if(execution_policy(par(task)), IteratorTag());
}

void count_if_test()
{
    test_count_if<std::random_access_iterator_tag>();
    test_count_if<std::forward_iterator_tag>();
    test_count_if<std::input_iterator_tag>();
}

////////////////////////////////////////////////////////////////////////////
template <typename ExPolicy, typename IteratorTag>
void test_count_if_exception(ExPolicy const& policy, IteratorTag)
{
    BOOST_STATIC_ASSERT(hpx::parallel::is_execution_policy<ExPolicy>::value);

    typedef std::vector<std::size_t>::iterator base_iterator;
    typedef test::decorated_iterator<base_iterator, IteratorTag>
        decorated_iterator;
    std::vector<std::size_t> c(10007);
    std::iota(boost::begin(c), boost::end(c), std::rand());

    bool caught_exception = false;
    try {
        //pred should never proc, so simple 'returns true'
        hpx::parallel::count_if(policy,
            decorated_iterator(
                boost::begin(c),
                [](){ throw std::runtime_error("test"); }),
            decorated_iterator(boost::end(c)),
            [](std::size_t val) { return true; });
        HPX_TEST(false);
    }
    catch(hpx::exception_list const& e) {
        caught_exception = true;
        test::test_num_exceptions<ExPolicy, IteratorTag>::call(policy, e);
    }
    catch(...) {
        HPX_TEST(false);
    }

    HPX_TEST(caught_exception);
}

template <typename ExPolicy, typename IteratorTag>
void test_count_if_exception_async(ExPolicy const& p, IteratorTag)
{
    typedef std::vector<std::size_t>::iterator base_iterator;
    typedef std::vector<std::size_t>::difference_type diff_type;
    typedef test::decorated_iterator<base_iterator, IteratorTag>
        decorated_iterator;

    std::vector<std::size_t> c(10007);
    std::fill(boost::begin(c), boost::end(c), 10);

    bool caught_exception = false;
    try {
        hpx::future<diff_type> f =
            hpx::parallel::count_if(p,
                decorated_iterator(
                    boost::begin(c),
                    [](){ throw std::runtime_error("test"); }),
                decorated_iterator(boost::end(c)),
                [](std::size_t val) { return true; });
        f.get();

        HPX_TEST(false);
    }
    catch(hpx::exception_list const& e) {
        caught_exception = true;
        test::test_num_exceptions<ExPolicy, IteratorTag>::call(p, e);
    }
    catch(...) {
        HPX_TEST(false);
    }

    HPX_TEST(caught_exception);
}

template <typename IteratorTag>
void test_count_if_exception()
{
    using namespace hpx::parallel;

    // If the execution policy object is of type vector_execution_policy,
    // std::terminate shall be called. therefore we do not test exceptions
    // with a vector execution policy
    test_count_if_exception(seq, IteratorTag());
    test_count_if_exception(par, IteratorTag());

    test_count_if_exception_async(seq(task), IteratorTag());
    test_count_if_exception_async(par(task), IteratorTag());

    test_count_if_exception(execution_policy(seq), IteratorTag());
    test_count_if_exception(execution_policy(par), IteratorTag());

    test_count_if_exception(execution_policy(seq(task)), IteratorTag());
    test_count_if_exception(execution_policy(par(task)), IteratorTag());
}

void count_if_exception_test()
{
    test_count_if_exception<std::random_access_iterator_tag>();
    test_count_if_exception<std::forward_iterator_tag>();
}

//////////////////////////////////////////////////////////////////////////////
template <typename ExPolicy, typename IteratorTag>
void test_count_if_bad_alloc(ExPolicy const& policy, IteratorTag)
{
    BOOST_STATIC_ASSERT(hpx::parallel::is_execution_policy<ExPolicy>::value);

    typedef std::vector<std::size_t>::iterator base_iterator;
    typedef test::decorated_iterator<base_iterator, IteratorTag>
        decorated_iterator;

    std::vector<std::size_t> c(10007);
    std::iota(boost::begin(c), boost::end(c), std::rand());

    bool caught_bad_alloc = false;
    try {
        hpx::parallel::count_if(policy,
            decorated_iterator(
                boost::begin(c),
                [](){ throw std::bad_alloc(); }),
            decorated_iterator(boost::end(c)),
            [](std::size_t v) { return true; });
        HPX_TEST(false);
    }
    catch (std::bad_alloc const&) {
        caught_bad_alloc = true;
    }
    catch (...) {
        HPX_TEST(false);
    }

    HPX_TEST(caught_bad_alloc);
}

template <typename ExPolicy, typename IteratorTag>
void test_count_if_bad_alloc_async(ExPolicy const& p, IteratorTag)
{
    typedef std::vector<std::size_t>::iterator base_iterator;
    typedef std::vector<std::size_t>::difference_type diff_type;
    typedef test::decorated_iterator<base_iterator, IteratorTag>
        decorated_iterator;

    std::vector<std::size_t> c(10007);
    std::iota(boost::begin(c), boost::end(c), std::rand());

    bool caught_bad_alloc = false;
    try {
        hpx::future<diff_type> f =
            hpx::parallel::count_if(p,
                decorated_iterator(
                    boost::begin(c),
                    [](){ throw std::bad_alloc(); }),
                decorated_iterator(boost::end(c)),
                [](std::size_t v) { return true; });

        f.get();

        HPX_TEST(false);
    }
    catch(std::bad_alloc const&) {
        caught_bad_alloc = true;
    }
    catch(...) {
        HPX_TEST(false);
    }

    HPX_TEST(caught_bad_alloc);
}

template <typename IteratorTag>
void test_count_if_bad_alloc()
{
    using namespace hpx::parallel;

    // If the execution policy object is of type vector_execution_policy,
    // std::terminate shall be called. therefore we do not test exceptions
    // with a vector execution policy
    test_count_if_bad_alloc(seq, IteratorTag());
    test_count_if_bad_alloc(par, IteratorTag());

    test_count_if_bad_alloc_async(seq(task), IteratorTag());
    test_count_if_bad_alloc_async(par(task), IteratorTag());

    test_count_if_bad_alloc(execution_policy(seq), IteratorTag());
    test_count_if_bad_alloc(execution_policy(par), IteratorTag());

    test_count_if_bad_alloc(execution_policy(seq(task)), IteratorTag());
    test_count_if_bad_alloc(execution_policy(par(task)), IteratorTag());
}

void count_if_bad_alloc_test()
{
    test_count_if_bad_alloc<std::random_access_iterator_tag>();
    test_count_if_bad_alloc<std::forward_iterator_tag>();
    test_count_if_bad_alloc<std::input_iterator_tag>();
}

int hpx_main()
{
    count_if_test();
    count_if_exception_test();
    count_if_bad_alloc_test();
    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    std::vector<std::string> cfg;
    cfg.push_back("hpx.os_threads=" +
        boost::lexical_cast<std::string>(hpx::threads::hardware_concurrency()));

    HPX_TEST_EQ_MSG(hpx::init(argc, argv, cfg), 0,
        "HPX main exited with non-zero status");

    return hpx::util::report_errors();
}
