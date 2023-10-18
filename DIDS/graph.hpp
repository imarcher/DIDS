#ifndef DIDS_GRAPH_H
#define DIDS_GRAPH_H

#include <cstdint>
#include "types.hpp"
#include "space_hnsw.hpp"


namespace dids {

    template<typename Helper, typename PointType, typename SpaceType>
    class GraphImpl : public Graph {
        using graph_type = hnswlib::HierarchicalNSW<float>;
        static constexpr int ef_construction = 40;
        static constexpr int m = 16;

    public:
        explicit GraphImpl(const uint32_t point_num, const uint32_t point_length) : space(new SpaceType(point_length)),
                                                                                    graph_(new graph_type(space,
                                                                                                          point_num, m,
                                                                                                          ef_construction)) {}

        explicit GraphImpl(const uint32_t point_length, const std::string &filename) : space(
                new SpaceType(point_length)),
                                                                                       graph_(new graph_type(space,
                                                                                                             filename)) {}

        ~GraphImpl() override {
            delete space;
            delete graph_;
        }

        inline void insert(void *point, uint32_t id) override {
            graph_->addPoint(point, id);
        }

        inline std::vector<std::pair<float, uint32_t>>
        search(void *search_point, uint32_t k) override {
            return graph_->searchKnn(search_point, k);
        }

        inline void saveIndex(const std::string &filename) override {
            graph_->saveIndex(filename);
        }

    private:
        SpaceType *space;
        graph_type *graph_;
    };

    template<typename Helper, typename PointType>
    struct GraphFactory {
        static Graph *createTSVec(const uint32_t point_num) {
            using space_type = EDSpace<Helper>;
            return new GraphImpl<Helper, PointType, space_type>(point_num, PointType::ts_length);
        }

        static Graph *createTSVec(const std::string &filename) {
            using space_type = EDSpace<Helper>;
            return new GraphImpl<Helper, PointType, space_type>(PointType::ts_length, filename);
        }

    };

}


#endif
