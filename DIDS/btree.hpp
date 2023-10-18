#ifndef DIDS_BTREE_H
#define DIDS_BTREE_H

#include "string"
#include "types.hpp"


namespace dids {
    namespace btree {

        class LeafNode {
        public:
            uint32_t len;
            uint64_t file_pos;
        };

        struct InternalKey {
            float dis;
            bool is_leaf;
            void *p;
        };

        struct InternalNode {
            InternalNode() : len(0) {}

            uint32_t len;
            std::vector<InternalKey> keys;
        };


        class BTree {
        public:

            BTree() {}

            ~BTree() {
                delete[] leaf_nodes;
                for (auto i: internal_nodes_vec) {
                    delete[] i;
                }
            }

            void build(const float *inputs, uint64_t input_size, const uint64_t now_file_pos) {
                btree_file_pos = now_file_pos;
                all_size = input_size;
                radius = inputs[input_size - 1];
                uint64_t n = all_size;
                uint64_t num_leafs = (n - 1) / btree_max_leaf_size + 1;
                leaf_nums = num_leafs;
                leaf_nodes = new LeafNode[num_leafs];
                internal_keys = new InternalKey[num_leafs];
                auto tmp_internal_keys = internal_keys;

                uint64_t file_pos = now_file_pos;
                for (int i = 0; i < num_leafs - 1; i++) {
                    leaf_nodes[i].len = btree_max_leaf_size;
                    leaf_nodes[i].file_pos = file_pos;
                    tmp_internal_keys[i].dis = inputs[file_pos - now_file_pos];
                    tmp_internal_keys[i].is_leaf = true;
                    tmp_internal_keys[i].p = &leaf_nodes[i];
                    file_pos += btree_max_leaf_size;
                }
                leaf_nodes[num_leafs - 1].len = n - btree_max_leaf_size * (num_leafs - 1);
                leaf_nodes[num_leafs - 1].file_pos = file_pos;
                tmp_internal_keys[num_leafs - 1].dis = inputs[file_pos - now_file_pos];
                tmp_internal_keys[num_leafs - 1].is_leaf = true;
                tmp_internal_keys[num_leafs - 1].p = &leaf_nodes[num_leafs - 1];


                uint64_t num_keys = num_leafs;
                uint64_t num_output_nodes = (num_keys - 1) / btree_max_leaf_size + 1;
                InternalNode *new_internal_nodes = new InternalNode[num_output_nodes];
                internal_nodes_vec.push_back(new_internal_nodes);
                BuildInternalNode(tmp_internal_keys, num_keys, new_internal_nodes, num_output_nodes);
//                delete[] tmp_internal_keys;

                while (num_output_nodes > 1) {

                    num_keys = num_output_nodes;
                    num_output_nodes = (num_keys - 1) / btree_max_leaf_size + 1;
                    tmp_internal_keys = new InternalKey[num_keys];
                    for (int i = 0; i < num_keys; i++) {
                        tmp_internal_keys[i].dis = new_internal_nodes[i].keys[0].dis;
                        tmp_internal_keys[i].is_leaf = false;
                        tmp_internal_keys[i].p = &new_internal_nodes[i];
                    }
                    new_internal_nodes = new InternalNode[num_output_nodes];
                    internal_nodes_vec.push_back(new_internal_nodes);
                    BuildInternalNode(tmp_internal_keys, num_keys, new_internal_nodes, num_output_nodes);
                    delete[] tmp_internal_keys;
                }

                root = new_internal_nodes;

            }

            void build(FILE *file) {

                uint64_t num_leafs = leaf_nums;
                leaf_nodes = new LeafNode[num_leafs];
                internal_keys = new InternalKey[num_leafs];
                auto tmp_internal_keys = internal_keys;
                fread(leaf_nodes, sizeof(LeafNode), num_leafs, file);
                fread(internal_keys, sizeof(InternalKey), num_leafs, file);
                for (int i = 0; i < num_leafs; i++) {
                    internal_keys[i].p = leaf_nodes + i;
                }

                uint64_t num_keys = num_leafs;
                uint64_t num_output_nodes = (num_keys - 1) / btree_max_leaf_size + 1;
                InternalNode *new_internal_nodes = new InternalNode[num_output_nodes];
                internal_nodes_vec.push_back(new_internal_nodes);
                BuildInternalNode(tmp_internal_keys, num_keys, new_internal_nodes, num_output_nodes);
                delete[] tmp_internal_keys;

                while (num_output_nodes > 1) {

                    num_keys = num_output_nodes;
                    num_output_nodes = (num_keys - 1) / btree_max_leaf_size + 1;
                    tmp_internal_keys = new InternalKey[num_keys];
                    for (int i = 0; i < num_keys; i++) {
                        tmp_internal_keys[i].dis = new_internal_nodes[i].keys[0].dis;
                        tmp_internal_keys[i].is_leaf = false;
                        tmp_internal_keys[i].p = &new_internal_nodes[i];
                    }
                    new_internal_nodes = new InternalNode[num_output_nodes];
                    internal_nodes_vec.push_back(new_internal_nodes);
                    BuildInternalNode(tmp_internal_keys, num_keys, new_internal_nodes, num_output_nodes);
                    delete[] tmp_internal_keys;
                }

                root = new_internal_nodes;

            }

            // return the begin and end pos that cannot be pruned
            std::pair<uint64_t, uint64_t>
            search(const float dis_to_ref_obj, const float top_sqrt, float *tmp_leaf_for_search_, FILE *btree_file) {
                //not include
                btree_file_for_search = btree_file;
                tmp_leaf_for_search = tmp_leaf_for_search_;
                const float begin_dis = dis_to_ref_obj - top_sqrt;
                const float end_dis = dis_to_ref_obj + top_sqrt;
                auto begin_pos = dfs_begin(begin_dis, root, false);
                auto end_pos = dfs_end(end_dis, root, false);
                return {begin_pos, end_pos};
            }


        private:


            static void BuildInternalNode(InternalKey *input_internal_keys, uint64_t num_input_keys,
                                          InternalNode *output_internal_nodes, uint64_t num_output_nodes) {

                for (int i = 0; i < num_output_nodes - 1; i++) {
                    output_internal_nodes[i].keys.resize(btree_max_leaf_size);
                    output_internal_nodes[i].len = btree_max_leaf_size;
                    memcpy(output_internal_nodes[i].keys.data(), input_internal_keys + i * btree_max_leaf_size,
                           sizeof(InternalKey) * btree_max_leaf_size);
                }

                uint64_t len = num_input_keys - (num_output_nodes - 1) * btree_max_leaf_size;
                output_internal_nodes[num_output_nodes - 1].keys.resize(len);
                output_internal_nodes[num_output_nodes - 1].len = len;
                memcpy(output_internal_nodes[num_output_nodes - 1].keys.data(),
                       input_internal_keys + (num_output_nodes - 1) * btree_max_leaf_size, sizeof(InternalKey) * len);
            }

            // find first > dis
            uint64_t dfs_begin(const float dis, const void *node, const bool is_leaf) const {
                if (is_leaf) {
                    auto leaf_node = (LeafNode *) node;
                    fseek(btree_file_for_search, leaf_node->file_pos * sizeof(float), SEEK_SET);
                    fread(tmp_leaf_for_search, sizeof(float), leaf_node->len, btree_file_for_search);

                    uint64_t l = 0, r = leaf_node->len - 1;
                    while (l < r) {
                        uint64_t mid = (l + r) >> 1;
                        if (tmp_leaf_for_search[mid] <= dis) {
                            l = mid + 1;
                        } else {
                            r = mid;
                        }
                    }
                    if (tmp_leaf_for_search[l] <= dis) return leaf_node->file_pos + leaf_node->len;
                    else return leaf_node->file_pos + l;
                } else {
                    auto internal_node = (InternalNode *) node;
                    uint64_t l = 0, r = internal_node->len - 1;
                    while (l < r) {
                        uint64_t mid = (l + r + 1) >> 1;
                        if (internal_node->keys[mid].dis <= dis) {
                            l = mid;
                        } else {
                            r = mid - 1;
                        }
                    }
                    return dfs_begin(dis, internal_node->keys[l].p, internal_node->keys[l].is_leaf);
                }
            }

            uint64_t dfs_end(const float dis, const void *node, const bool is_leaf) const {
                if (is_leaf) {
                    auto leaf_node = (LeafNode *) node;
                    fseek(btree_file_for_search, leaf_node->file_pos * sizeof(float), SEEK_SET);
                    fread(tmp_leaf_for_search, sizeof(float), leaf_node->len, btree_file_for_search);

                    uint64_t l = 0, r = leaf_node->len - 1;
                    while (l < r) {
                        uint64_t mid = (l + r + 1) >> 1;
                        if (tmp_leaf_for_search[mid] < dis) {
                            l = mid;
                        } else {
                            r = mid - 1;
                        }
                    }
                    if (tmp_leaf_for_search[l] >= dis) return leaf_node->file_pos - 1;
                    else return leaf_node->file_pos + l;
                } else {
                    auto internal_node = (InternalNode *) node;
                    uint64_t l = 0, r = internal_node->len - 1;
                    while (l < r) {
                        uint64_t mid = (l + r + 1) >> 1;
                        if (internal_node->keys[mid].dis < dis) {
                            l = mid;
                        } else {
                            r = mid - 1;
                        }
                    }
                    return dfs_end(dis, internal_node->keys[l].p, internal_node->keys[l].is_leaf);
                }
            }


            const std::string btree_filename;
            InternalNode *root;
            float *tmp_leaf_for_search;
            FILE *btree_file_for_search;
            std::vector<InternalNode *> internal_nodes_vec;

        public:
            static constexpr int btree_max_leaf_size = Const::btree_max_leaf_size;
            uint64_t btree_file_pos;
            uint64_t all_size;
            uint64_t leaf_nums;
            float radius;
            LeafNode *leaf_nodes;
            InternalKey *internal_keys;

        };
    }
}

#endif
