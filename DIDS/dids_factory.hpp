#ifndef DIDS_DIDSFACTORY_H
#define DIDS_DIDSFACTORY_H


#include "dids_impl.hpp"

namespace dids {

    template<uint32_t TSLength, uint32_t SAXLength>
    class DIDSFactory {
        using ts_vec_type = TSVec<TSLength>;
        using paa_vec_type = PAAVec<SAXLength>;
        using sax_vec_type = SAXVec<SAXLength>;
        using tree_type = DIDSImpl<ts_vec_type, paa_vec_type, sax_vec_type>;
    public:
        /**
         * firstly create DIDS index
         * @param ts_num the tuples of dataset
         * @param ref_objs_size the number of reference objects
         * @param approximate_leaf_size the leaf size of the binary tree
         * @param ts_buffer_size_for_read the size of read buffer for building
         * @param ts_buffer_size_per_ref_obj the size of write buffer per reference object for building
         * @param data_name the name of dataset
         * @param input_filename the file path of dataset
         * @param output_directory the directory of the index to store
         */
        static void create(const uint64_t ts_num, const uint32_t ref_objs_size, const uint32_t approximate_leaf_size,
                           const uint64_t ts_buffer_size_for_read,
                           const uint64_t ts_buffer_size_per_ref_obj,
                           const std::string &data_name,
                           const std::string &input_filename, const std::string &output_directory) {
            tree_type(ts_num, ref_objs_size, approximate_leaf_size, ts_buffer_size_for_read, ts_buffer_size_per_ref_obj,
                      data_name, input_filename, output_directory);
        }

        /**
         * create DIDS index from the disk files
         * @param data_name the name of dataset
         * @param output_directory the directory of the index
         * @return a DIDS index
         */
        static DIDS *createFromIndex(const std::string &data_name, const std::string &output_directory) {
            return new tree_type(data_name, output_directory);
        }
    };
}


#endif
