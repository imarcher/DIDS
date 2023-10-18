#include "DIDS/dids_factory.hpp"
#include "random_data.h"

#include <iostream>

using namespace std;

int main() {

    // generate data
    const string input_filename = "./ts.bin";
    const string output_directory = "./random_index/";
    const string data_name = "random";
    const uint64_t ts_length = 128;
    const uint64_t sax_length = 16;
    const uint64_t ts_num = 1000000;
    const uint64_t k = 10;

    RandomData::generate_data(ts_length, ts_num, input_filename);
    auto file = fopen(input_filename.c_str(), "r");
    float query[ts_length];
    fread(query, sizeof(float), ts_length, file);
    fclose(file);

    // build the DIDS
    dids::DIDSFactory<ts_length, sax_length>::create(ts_num, 1000, 500, 10000, 100, data_name, input_filename,
                                                     output_directory);
    // load the DIDS
    auto dids_index = dids::DIDSFactory<ts_length, sax_length>::createFromIndex(data_name, output_directory);
    // approximate search
    auto approximate_ans = dids_index->approximateSearch(query, k, 5);
    // exact search
    auto exact_ans = dids_index->search(query, k, 5);

    cout << "recall: " << (float) getRecallNum(exact_ans[k - 1].first, approximate_ans) / k * 100 << "%" << endl;
    delete dids_index;
}