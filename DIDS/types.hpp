#ifndef DIDS_TYPES_H
#define DIDS_TYPES_H

#include <cstdint>
#include <algorithm>
#include <vector>
#include <cstring>
#include "breakpoints.hpp"
#include <thread>
#include "immintrin.h"
#include "configure.hpp"

namespace dids {

    class Graph {
    public:
        virtual ~Graph() = default;

        virtual inline void insert(void *point, uint32_t id) = 0;

        virtual inline std::vector<std::pair<float, uint32_t>>
        search(void *search_point, uint32_t k) = 0;

        virtual inline void saveIndex(const std::string &filename) = 0;
    };

    template<uint32_t TSLength>
    struct TSVec {
        static constexpr uint32_t ts_length = TSLength;
        float val[TSLength];
    };

    template<uint32_t PAALength>
    struct PAAVec {
        static constexpr uint32_t paa_length = PAALength;
        float val[PAALength];
    };

    template<uint32_t SAXLength>
    struct SAXVec {
        using sax_type = unsigned char;
        static constexpr uint32_t sax_length = SAXLength;
        sax_type val[SAXLength];
    };


    struct SketchUnit {
        float upper_bound;
        float lower_bound;
    };

    template<typename TSVecType>
    struct Sketch {
        static constexpr uint32_t sketch_length = TSVecType::ts_length;

        inline float getDisI(uint32_t id) {
            return val[id].upper_bound - val[id].lower_bound;
        }

        inline void updateI(float new_val, uint32_t id) {
            val[id].upper_bound = std::max(new_val, val[id].upper_bound);
            val[id].lower_bound = std::min(new_val, val[id].lower_bound);
        }

        inline void init() {
            for (int i = 0; i < sketch_length; i++) {
                val[i].upper_bound = -MAXFLOAT;
                val[i].lower_bound = MAXFLOAT;
            }
        }

        SketchUnit val[sketch_length];
    };


    template<typename TSVecType, typename PAAVecType, typename SAXVecType>
    struct TypesHelper {
        static constexpr uint32_t num_per_segment = TSVecType::ts_length / SAXVecType::sax_length;
        using ts_vec_type = TSVecType;

        static inline void getPAAFromTS(TSVecType &ts_vec, PAAVecType &paa_vec) {
            for (int i = 0, s = 0; s < SAXVecType::sax_length; i += num_per_segment, s++) {
                paa_vec.val[s] = 0;
                for (int j = i; j < i + num_per_segment; j++) {
                    paa_vec.val[s] += ts_vec.val[j];
                }
                paa_vec.val[s] /= num_per_segment;
            }
        }

        static inline void getSAXFromTS(TSVecType &ts_vec, SAXVecType &sax_vec) {
            PAAVecType tmp_paa;
            for (int i = 0, s = 0; s < SAXVecType::sax_length; i += num_per_segment, s++) {
                tmp_paa.val[s] = 0;
                for (int j = i; j < i + num_per_segment; j++) {
                    tmp_paa.val[s] += ts_vec.val[j];
                }
                tmp_paa.val[s] /= num_per_segment;
            }

            for (int s = 0; s < SAXVecType::sax_length; s++) {
                int l = 0, r = 256;
                while (l < r) {
                    int mid = (l + r + 1) >> 1;
                    if (sax_breakpoints[mid] < tmp_paa.val[s]) l = mid;
                    else r = mid - 1;
                }
                sax_vec.val[s] = l;
            }
        }

        // not sqrt
        static inline float getDisBetweenTSAndTS(const TSVecType &ts_vec_1, const TSVecType &ts_vec_2) {
#if DIDS_PARALLELISM == 0
            float dis = 0;
            for (int i = 0; i < TSVecType::ts_length; i++) {
                float tmp = ts_vec_1.val[i] - ts_vec_2.val[i];
                dis += tmp * tmp;
            }
            return dis;
#else

            float dis = 0;

            int len = TSVecType::ts_length / 8 * 8;


            __m256 xfsSum = _mm256_setzero_ps();
            __m256 a;    // 加载.
            __m256 b;


            for (int i = 0; i < len; i += 8) {
                a = _mm256_loadu_ps(&ts_vec_1.val[i]);
                b = _mm256_loadu_ps(&ts_vec_2.val[i]);
                __m256 c = _mm256_sub_ps(a, b);
                __m256 d = _mm256_mul_ps(c, c);
                xfsSum = _mm256_add_ps(xfsSum, d);
            }
            const auto *q = (const float *) &xfsSum;
            dis = q[0] + q[1] + q[2] + q[3] + q[4] + q[5] + q[6] + q[7];

            for (int i = len; i < TSVecType::ts_length; i++) {
                dis += (ts_vec_1.val[i] - ts_vec_2.val[i]) * (ts_vec_1.val[i] - ts_vec_2.val[i]);
            }
            return dis;
#endif
        }

        // not sqrt
        static inline float getDisBetweenSAXAndPAA(const SAXVecType &sax_vec, const PAAVecType &paa_vec) {
#if DIDS_PARALLELISM
            float dis = 0;
            int len = SAXVecType::sax_length / 8 * 8;


            float breakpoints_lower[SAXVecType::sax_length];
            float breakpoints_upper[SAXVecType::sax_length];

            for (int i = 0; i < SAXVecType::sax_length; i++) {
                auto region = sax_vec.val[i];
                breakpoints_lower[i] = sax_breakpoints[region];
                breakpoints_upper[i] = sax_breakpoints[region + 1];
            }

            __m256 dis8 = _mm256_setzero_ps();

            for (int i = 0; i < len; i += 8) {
                __m256 paa8 = _mm256_loadu_ps(&paa_vec.val[i]);

                __m256 lower8 = _mm256_loadu_ps(breakpoints_lower + i);
                __m256 to_sub = _mm256_sub_ps(lower8, paa8);
                __m256 f = _mm256_mul_ps(to_sub, to_sub);
                __m256 mask = _mm256_cmp_ps(lower8, paa8, _CMP_GT_OQ);
                dis8 = _mm256_add_ps(dis8, _mm256_and_ps(f, mask));

                __m256 upper8 = _mm256_loadu_ps(breakpoints_upper + i);
                to_sub = _mm256_sub_ps(upper8, paa8);
                f = _mm256_mul_ps(to_sub, to_sub);
                mask = _mm256_cmp_ps(upper8, paa8, _CMP_LT_OQ);
                dis8 = _mm256_add_ps(dis8, _mm256_and_ps(f, mask));

            }

            const auto *q = (const float *) &dis8;
            dis = q[0] + q[1] + q[2] + q[3] + q[4] + q[5] + q[6] + q[7];

            for (int i = len; i < SAXVecType::sax_length; i++) {
                if (breakpoints_lower[i] > paa_vec.val[i]) {
                    dis += (breakpoints_lower[i] - paa_vec.val[i]) * (breakpoints_lower[i] - paa_vec.val[i]);
                } else if (breakpoints_upper[i] < paa_vec.val[i]) {
                    dis += (breakpoints_upper[i] - paa_vec.val[i]) * (breakpoints_upper[i] - paa_vec.val[i]);
                }
            }

            return dis * num_per_segment;

#else
            float dis = 0;
            for (int i = 0; i < SAXVecType::sax_length; i++) {
                auto region = sax_vec.val[i];
                float breakpoint_lower = sax_breakpoints[region];
                float breakpoint_upper = sax_breakpoints[region + 1];

                if (breakpoint_lower > paa_vec.val[i]) {
                    dis += (breakpoint_lower - paa_vec.val[i]) * (breakpoint_lower - paa_vec.val[i]);
                } else if (breakpoint_upper < paa_vec.val[i]) {
                    dis += (breakpoint_upper - paa_vec.val[i]) * (breakpoint_upper - paa_vec.val[i]);
                }
            }
            return dis * num_per_segment;
#endif
        }

    };
}

#endif
