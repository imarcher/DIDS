# DIDS

DIDS: Double Indices and Double Summarizations for Fast Similarity Search

DIDS is a disk-based index supporting both approximate and exact searches

## Prerequisites

* GCC 9.4+ with OpenMP
* CMake 3.16+
* TCMalloc

## Datasets

## API

DIDS consists entirely of `.hpp` files. 
You can integrate it by including `DIDS/dids_factory.hpp`.
DIDS has a total of 4 APIs, two for building and two for querying.

* `DIDSFactory<TSLength, SAXLength>::create(ts_num, ref_objs_size, approximate_leaf_size, ts_buffer_size_for_read, ts_buffer_size_per_ref_obj, data_name, input_filename, output_directory)` create DIDS index
  * `ts_num` the tuples of dataset
  * `ref_objs_size` the number of reference objects
  * `approximate_leaf_size` the leaf size of the binary tree
  * `ts_buffer_size_for_read` the size of read buffer for building
  * `ts_buffer_size_per_ref_obj` the size of write buffer per reference object for building
  * `data_name` the name of dataset
  * `input_filename` the file path of dataset
  * `output_directory` the directory of the index to store
* `DIDSFactory<TSLength, SAXLength>::createFromIndex(data_name, output_directory)` create DIDS index from disk files
* `approximateSearch(search_ts_vec, k, approximate_search_node_num)` approximate search
  * `search_ts_vec` the query data series
  * `k` the k of k-NN
  * `approximate_search_node_num` the graph nodes to search in the approximate search
* `search(search_ts_vec, k, approximate_search_node_num)` exact search

## A running example

Here's an example in `example.cpp` that generates 1,000,000 data sequences, builds a DIDS index, loads the DIDS index from disk, and performs approximate and exact queries.

```c++
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
```
You can run it with:
```
mkdir build
cd build
cmake ..
make
./example
```

## Multi-threading

DIDS has a multi-threading version (partial), which can be enabled or disabled in the `DIDS/configure.hpp` file.















