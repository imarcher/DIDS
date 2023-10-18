#ifndef DIDS_CONFIGURE_H
#define DIDS_CONFIGURE_H

#include <ctime>
#include <sys/time.h>

// Partial parallelism, incomplete
#define DIDS_PARALLELISM 1
#define DIDS_PRINT_INFO 1
#define DIDS_HDD 0


// All unimportant parameters
struct Const {
    static constexpr uint32_t thread_num = 32;
    static constexpr uint32_t btree_max_leaf_size = 1000;
    static constexpr uint32_t btree_write_ts_batch_size = 800;
    static constexpr uint32_t btree_read_sax_batch_size = 2560;
};

#if DIDS_PRINT_INFO
#define EXACT_SEARCH_START gettimeofday(&exact_search_time_start, NULL);
#define EXACT_SEARCH_END  gettimeofday(&current_time, NULL); \
                                      start_time = exact_search_time_start.tv_sec*1000000 + (exact_search_time_start.tv_usec); \
                                      end_time  = current_time.tv_sec*1000000 + (current_time.tv_usec); \
                                      exact_search_time = (end_time - start_time) / 1000;

#define APPROXIMATE_SEARCH_START gettimeofday(&approximate_search_time_start, NULL);
#define APPROXIMATE_SEARCH_END  gettimeofday(&current_time, NULL); \
                                      start_time = approximate_search_time_start.tv_sec*1000000 + (approximate_search_time_start.tv_usec); \
                                      end_time  = current_time.tv_sec*1000000 + (current_time.tv_usec); \
                                      approximate_search_time = (end_time - start_time) / 1000;

#define COUNT_READ_SAX(num) read_sax_count += num;
#define COUNT_READ_TS(num) read_ts_count += num;

#if DIDS_PARALLELISM
#define CLEAR_INFO read_sax_count = 0, read_ts_count = 0, exact_search_time = 0, approximate_search_time = 0; \
                    memset(read_sax_count_i, 0, sizeof(uint64_t) * Const::thread_num);                                        \
                    memset(read_ts_count_i, 0, sizeof(uint64_t) * Const::thread_num);
#define COUNT_READ_SAX_I(num, i) read_sax_count_i[i] += num;
#define COUNT_READ_TS_I(num, i) read_ts_count_i[i] += num;

#else
#define CLEAR_INFO read_sax_count = 0, read_ts_count = 0, exact_search_time = 0, approximate_search_time = 0;
#endif


#define PRINT_INFO std::cout<<"a query:"<<std::endl;\
std::cout<<"    approximate search time: "<<  approximate_search_time <<"ms"<<std::endl; \
std::cout<<"    exact search time: "<<  exact_search_time <<"ms"<<std::endl; \
std::cout<<"    read sax: "<<  (double)read_sax_count * 100 / total_ts_num <<"%"<<std::endl; \
std::cout<<"    read ts: "<<  (double)read_ts_count * 100 / total_ts_num<<"%"<<std::endl;

#else
#define EXACT_SEARCH_START
#define EXACT_SEARCH_END

#define APPROXIMATE_SEARCH_START
#define APPROXIMATE_SEARCH_END

#define COUNT_READ_SAX(num)
#define COUNT_READ_TS(num)
#define CLEAR_INFO

#define PRINT_INFO

#endif


#endif

