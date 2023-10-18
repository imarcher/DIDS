#ifndef DIDS_RANDOMDATA_H
#define DIDS_RANDOMDATA_H


#include <cstdint>
#include <string>
#include "vector"

static int getRecallNum(const float real_ans, const std::vector<std::pair<float, uint64_t>> &find_ans) {
    int recall_num = 0;
    for (const auto & find_an : find_ans) {
        if (find_an.first < real_ans + 1e-4) {
            recall_num ++;
        }
    }
    return recall_num;
}

struct RandomData {
    static void generate_data(uint64_t ts_length, uint64_t ts_num, const std::string &filename_);
};


#endif //DIDS_RANDOMDATA_H
