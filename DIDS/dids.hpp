#ifndef DIDS_DIDS_H
#define DIDS_DIDS_H

#include <cstdint>
#include "vector"

namespace dids {
    class DIDS {
    public:
        virtual ~DIDS() = default;

        /**
         * exact search
         * @param search_ts_vec the query ts
         * @param k the k of k-NN
         * @param approximate_search_node_num the graph nodes to search in the approximate search
         * @return the distances and positions in files of answers
         */
        virtual std::vector<std::pair<float, uint64_t>>
        search(void *search_ts_vec, uint32_t k, uint32_t approximate_search_node_num) = 0;

        /**
         * approximate search
         * @param search_ts_vec the query ts
         * @param k the k of k-NN
         * @param approximate_search_node_num the graph nodes to search in the approximate search
         * @return the distances and positions in files of answers
         */
        virtual std::vector<std::pair<float, uint64_t>>
        approximateSearch(void *search_ts_vec, uint32_t k, uint32_t approximate_search_node_num) = 0;
    };
}


#endif
