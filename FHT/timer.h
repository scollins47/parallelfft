/**
 * timer.h
 *
 * Custom timing macros for serial/OpenMP programs. Uses omp_get_wtime() if
 * _OPENMP is defined and gettimeofday() otherwise.
 *
 * Example:
 *
 *      START_TIMER(tag1)
 *      do_stuff();
 *      STOP_TIMER(tag1)
 *
 *      START_TIMER(tag2)
 *      do_more_stuff();
 *      STOP_TIMER(tag2)
 *
 *      printf("tag1: %8.4fs  tag2: %8.4fs\n",
 *          GET_TIMER(tag1), GET_TIMER(tag2));
 */

#ifdef _OPENMP
#   include <omp.h>
#   define START_TIMER(X) double _timer_ ## X = omp_get_wtime();
#   define STOP_TIMER(X)  _timer_ ## X = omp_get_wtime() - (_timer_ ## X);
#   define GET_TIMER(X)   (_timer_ ## X)
#else
#   include <sys/time.h>
    struct timeval _tv;
#   define START_TIMER(X) gettimeofday(&_tv, NULL); \
        double _timer_ ## X = _tv.tv_sec+(_tv.tv_usec/1000000.0);
#   define STOP_TIMER(X)  gettimeofday(&_tv, NULL); \
        _timer_ ## X = _tv.tv_sec+(_tv.tv_usec/1000000.0) - (_timer_ ## X);
#   define GET_TIMER(X)   (_timer_ ## X)
#endif

