#ifndef DIDS_SPACE_HNSW_H
#define DIDS_SPACE_HNSW_H

#include "types.hpp"
#include "hnswlib/hnswalg.h"


namespace dids {


    template<typename Helper>
    static float getDisBetweenTSAndTS(const void *pVect1v, const void *pVect2v, const void *qty_ptr) {
        return Helper::getDisBetweenTSAndTS(*(typename Helper::ts_vec_type *) pVect1v,
                                            *(typename Helper::ts_vec_type *) pVect2v);
    }


    template<typename Helper>
    class EDSpace : public hnswlib::SpaceInterface<float> {
        hnswlib::DISTFUNC<float> fstdistfunc_;
        size_t data_size_;
        size_t dim_;

    public:
        EDSpace(size_t dim) {
            fstdistfunc_ = getDisBetweenTSAndTS<Helper>;
            dim_ = dim;
            data_size_ = sizeof(typename Helper::ts_vec_type);
        }

        size_t get_data_size() {
            return data_size_;
        }

        hnswlib::DISTFUNC<float> get_dist_func() {
            return fstdistfunc_;
        }

        void *get_dist_func_param() {
            return &dim_;
        }

        ~EDSpace() {}
    };
}


#endif
