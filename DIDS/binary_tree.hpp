#ifndef DIDS_BINARY_TREE_H
#define DIDS_BINARY_TREE_H

#include <cstdio>
#include <vector>
#include "types.hpp"
#include "iostream"
#include "set"
#include "btree.hpp"
#include <functional>


namespace dids {
    namespace binary_tree {


        template<typename TSVecType>
        struct InternalNode {
            InternalNode(uint16_t split_segment_id, float split_val,
                         void *left_child, void *right_child) :
                    left_is_leaf(true), right_is_leaf(true),
                    split_segment_id(split_segment_id), split_val(split_val),
                    left_child(left_child), right_child(right_child) {}

            bool left_is_leaf;
            bool right_is_leaf;
            uint16_t split_segment_id;
            // left <= split_val, right > split_val
            float split_val;
            void *left_child;
            void *right_child;
        };


        template<typename TSVecType>
        class LeafNodeForBuild {
        public:
            using sketch_type = Sketch<TSVecType>;
            using internal_node_type = InternalNode<TSVecType>;
            using sort_type = std::pair<float, uint64_t>;

        private:


            internal_node_type *split(TSVecType *ts_buffer) {
                std::pair<float, uint32_t> max_dis_id = {0, 1};
                for (int i = 0; i < sketch_type::sketch_length; i++) {
                    float dis = sketch.getDisI(i);
                    if (dis > max_dis_id.first) {
                        max_dis_id = {dis, i};
                    }
                }
                if (max_dis_id.first < 1e-6) {
                    std::cout << "too many same data series, can not split the node";
                    exit(2);
                }

                uint32_t split_id = max_dis_id.second;

                auto to_find_mid_mem = new sort_type[size_];
                for (int i = 0; i < size_; i++) {
                    to_find_mid_mem[i] = {ts_buffer[ts_vec_ids[i]].val[split_id], ts_vec_ids[i]};
                }
                std::sort(to_find_mid_mem, to_find_mid_mem + size_);
                float mid_val = to_find_mid_mem[size_ / 2].first;

                if (sketch.val[split_id].upper_bound - mid_val < 1e-6) {
                    uint64_t l = 0, r = size_ - 1;
                    while (l < r) {
                        uint64_t mid = (l + r + 1) >> 1;
                        if (to_find_mid_mem[mid].first < mid_val) {
                            l = mid;
                        } else {
                            r = mid - 1;
                        }
                    }
                    mid_val = to_find_mid_mem[l].first;
                }

                auto *left_child = this;
                auto *right_child = new LeafNodeForBuild(leaf_size);
                auto *new_internal_node = new internal_node_type(split_id, mid_val, left_child, right_child);

                //move data to child nodes and update the sketch of child nodes
                uint32_t old_size = size_;
                clear();

                for (int i = 0; i < old_size; i++) {
                    if (to_find_mid_mem[i].first <= mid_val) {
                        insert(to_find_mid_mem[i].second, ts_buffer);
                    } else {
                        right_child->insert(to_find_mid_mem[i].second, ts_buffer);
                    }
                }

                delete[] to_find_mid_mem;
                return new_internal_node;
            }

            void updateSketch(const uint64_t ts_vec_id, TSVecType *ts_buffer) {
                auto &ts_vec = ts_buffer[ts_vec_id];
                for (int i = 0; i < sketch_type::sketch_length; i++) {
                    sketch.updateI(ts_vec.val[i], i);
                }
            }

            void clear() {
                sketch.init();
                size_ = 0;
                ts_vec_ids.clear();
            }

        public:


            LeafNodeForBuild(const uint32_t leaf_size) : leaf_size(leaf_size) {
                ts_vec_ids.reserve(leaf_size);
                clear();
            }

            internal_node_type *insert(const uint64_t ts_vec_id, TSVecType *ts_buffer) {
                ts_vec_ids.push_back(ts_vec_id);
                size_++;
                updateSketch(ts_vec_id, ts_buffer);
                if (size_ >= leaf_size) {
                    return split(ts_buffer);
                }
                return nullptr;
            }

            TSVecType getCenter(TSVecType *ts_buffer) {
                TSVecType center;
                for (int j = 0; j < sketch_type::sketch_length; j++) {
                    center.val[j] = 0;
                }
                for (int i = 0; i < size_; i++) {
                    for (int j = 0; j < sketch_type::sketch_length; j++) {
                        center.val[j] += ts_buffer[ts_vec_ids[i]].val[j];
                    }
                }
                for (int j = 0; j < sketch_type::sketch_length; j++) {
                    center.val[j] /= size_;
                }
                return center;
            }

            sketch_type sketch;
            uint32_t leaf_size;
            uint32_t size_;
            std::vector<uint64_t> ts_vec_ids;
        };


        struct LeafNodeForSearch {
            uint64_t file_pos;
            uint64_t file_pos_end;

            bool operator<(const LeafNodeForSearch &t) const {
                return file_pos < t.file_pos;
            }
        };

        template<typename TSVecType>
        class BinaryTree {
            using leaf_node_for_build_type = LeafNodeForBuild<TSVecType>;
            using internal_node_type = InternalNode<TSVecType>;
        public:

            explicit BinaryTree(const uint32_t leaf_size) : root(new leaf_node_for_build_type(leaf_size)),
                                                            root_is_leaf(true) {
                leaf_nodes.push_back((leaf_node_for_build_type *) root);
            }

            ~BinaryTree() {
                for (auto node: internal_nodes) {
                    delete node;
                }
                for (auto node: leaf_nodes) {
                    delete node;
                }
            }

            void insert(const uint64_t ts_vec_id, TSVecType *ts_buffer) {
                insert(ts_vec_id, ts_buffer, root, root_is_leaf);
            }

            // write into file and build the search array
            void write(uint64_t file_pos, std::vector<LeafNodeForSearch> &search_leaf_nodes) {
                uint32_t leaf_node_num = leaf_nodes.size();
                for (int i = 0; i < leaf_node_num; i++) {
                    uint32_t node_len = leaf_nodes[i]->size_;
                    uint64_t max_id = 0;
                    uint64_t min_id = UINT64_MAX;
                    for (int j = 0; j < node_len; j++) {
                        max_id = std::max(max_id, leaf_nodes[i]->ts_vec_ids[j]);
                        min_id = std::min(min_id, leaf_nodes[i]->ts_vec_ids[j]);
                    }
                    search_leaf_nodes.push_back({min_id + file_pos, max_id + file_pos});
                }
            }

            void insertGraph(Graph *graph, uint64_t &leaf_id, TSVecType *ts_buffer) {
                uint32_t leaf_node_num = leaf_nodes.size();
                for (int i = 0; i < leaf_node_num; i++) {
                    auto center = leaf_nodes[i]->getCenter(ts_buffer);
                    graph->insert(&center, leaf_id + i);
                }
                leaf_id += leaf_node_num;
            }

        private:

            void insert(const uint64_t ts_vec_id, TSVecType *ts_buffer, void *&to_insert_node, bool &is_leaf) {
                if (is_leaf) {
                    auto real_to_insert_node = (leaf_node_for_build_type *) to_insert_node;
                    auto new_internal_node = real_to_insert_node->insert(ts_vec_id, ts_buffer);
                    if (new_internal_node) {
                        is_leaf = false;
                        to_insert_node = new_internal_node;
                        internal_nodes.push_back(new_internal_node);
                        leaf_nodes.push_back((leaf_node_for_build_type *) new_internal_node->right_child);
                    }
                } else {
                    auto real_to_insert_node = (internal_node_type *) to_insert_node;
                    if (ts_buffer[ts_vec_id].val[real_to_insert_node->split_segment_id] <=
                        real_to_insert_node->split_val) {
                        insert(ts_vec_id, ts_buffer, real_to_insert_node->left_child,
                               real_to_insert_node->left_is_leaf);
                    } else {
                        insert(ts_vec_id, ts_buffer, real_to_insert_node->right_child,
                               real_to_insert_node->right_is_leaf);
                    }
                }
            }


            void *root;
            bool root_is_leaf;
            std::vector<internal_node_type *> internal_nodes;
        public:
            std::vector<leaf_node_for_build_type *> leaf_nodes;
        };
    }
}
#endif
