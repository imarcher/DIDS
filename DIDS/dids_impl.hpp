#ifndef DIDS_DIDS_IMPL_H
#define DIDS_DIDS_IMPL_H

#include "btree.hpp"
#include "binary_tree.hpp"
#include "graph.hpp"
#include "dids.hpp"
#include "dkm/dkm_parallel.hpp"

namespace dids {

    class TopKHeap {
    public:
        TopKHeap(const uint32_t k) : k(k) {}

        void pushNotKnowK(float dis, uint64_t pos) {
            if (sizeIsK()) {
                if (dis < top) {
                    heap_.pop();
                    heap_.push({dis, pos});
                    top = heap_.top().first;
                    top_sqrt = sqrtf(top);
                }
            } else {
                heap_.push({dis, pos});
                if (sizeIsK()) {
                    top = heap_.top().first;
                    top_sqrt = sqrtf(top);
                }
            }
        }

        void push(float dis, uint64_t pos) {
            heap_.pop();
            heap_.push({dis, pos});
            top = heap_.top().first;
            top_sqrt = sqrtf(top);

        }

        bool sizeIsK() {
            return heap_.size() == k;
        }

        std::vector<std::pair<float, uint64_t>> getRes() {
            std::vector<std::pair<float, uint64_t>> res(k);
            for (int i = 0; i < k; i++) {
                res[i] = heap_.top();
                heap_.pop();
            }
            std::sort(res.begin(), res.end());
            return res;
        }

    private:
        std::priority_queue<std::pair<float, uint64_t>, std::vector<std::pair<float, uint64_t>>> heap_;
        const uint32_t k;
    public:
        float top;
        float top_sqrt;
    };


    template<typename TSVecType, typename PAAVecType, typename SAXVecType>
    class DIDSImpl : public DIDS {
        using types_helper = TypesHelper<TSVecType, PAAVecType, SAXVecType>;
        using binary_tree_type = binary_tree::BinaryTree<TSVecType>;
        using approximate_node_type = binary_tree::LeafNodeForSearch;

    public:


        DIDSImpl(const uint64_t ts_num, const uint32_t ref_objs_size, const uint32_t approximate_leaf_size,
                 const uint64_t ts_buffer_size_for_read,
                 const uint64_t ts_buffer_size_for_one_ref_obj,
                 const std::string &data_name,
                 const std::string &input_filename, const std::string &output_directory) :
                ts_filename(output_directory + data_name + "_ts.bin"),
                sax_filename(output_directory + data_name + "_sax.bin"),
                btree_filename(output_directory + data_name + "_btree.bin"),
                btree_leaf_nodes_filename(output_directory + data_name + "_btree_leaf_nodes.bin"),
                ref_objs_filename(output_directory + data_name + "_ref_objs.bin"),
                approximate_leaf_nodes_filename(output_directory + data_name + "_approximate_leaf_nodes.bin"),
                approximate_graph_filename(output_directory + data_name + "_approximate_graph.bin"),
                tmp_filename(output_directory + "_tmp_"),
                tmp_dis_filename(output_directory + "_tmp_dis_"),
                total_ts_num(ts_num),
                ref_objs_size(ref_objs_size), ref_objs(ref_objs_size),
                btrees(ref_objs_size) {

            system(("mkdir -p " + output_directory).c_str());

            fprintf(stderr, ">> select ref objs\n");

            selectRefObjs(ts_num, input_filename);

            fprintf(stderr, ">> group ts by ref objs\n");

            auto btree_all_size_vec = groupTSByRefObjs(ts_num, ts_buffer_size_for_read, ts_buffer_size_for_one_ref_obj,
                                                       input_filename);

            fprintf(stderr, ">> load ts within same ref objs, partition by binary tree and build the graph\n");
            buildApproximateStructs(ts_num, approximate_leaf_size, btree_all_size_vec);

            fprintf(stderr, ">> build finished\n");
        }


        DIDSImpl(const std::string &data_name, const std::string &output_directory) :
                ts_filename(output_directory + data_name + "_ts.bin"),
                sax_filename(output_directory + data_name + "_sax.bin"),
                btree_filename(output_directory + data_name + "_btree.bin"),
                btree_leaf_nodes_filename(output_directory + data_name + "_btree_leaf_nodes.bin"),
                ref_objs_filename(output_directory + data_name + "_ref_objs.bin"),
                approximate_leaf_nodes_filename(output_directory + data_name + "_approximate_leaf_nodes.bin"),
                approximate_graph_filename(output_directory + data_name + "_approximate_graph.bin"),
                tmp_filename(output_directory + "_tmp_"),
                tmp_dis_filename(output_directory + "_tmp_dis_") {

            fprintf(stderr, ">> load ref objs\n");
            loadRefObjs();

            fprintf(stderr, ">> load BTrees\n");
            loadBtrees();

            fprintf(stderr, ">> load approximate structures\n");
            loadApproximateStructs();

            fprintf(stderr, ">> build finished\n");
        }

        ~DIDSImpl() override {
            delete approximate_graph;
        }

        std::vector<std::pair<float, uint64_t>>
        search(void *search_ts_vec, uint32_t k, uint32_t approximate_search_node_num) override {
            CLEAR_INFO
            return search(*(TSVecType *) search_ts_vec, k, approximate_search_node_num);
        }


        std::vector<std::pair<float, uint64_t>>
        approximateSearch(void *search_ts_vec, uint32_t k, uint32_t approximate_search_node_num) override {
            CLEAR_INFO
            return approximateSearch(*(TSVecType *) search_ts_vec, k, approximate_search_node_num);
        }


    private:

        std::vector<std::pair<float, uint64_t>>
        search(TSVecType &search_ts_vec, uint32_t k, uint32_t approximate_search_node_num) {

            FILE *ts_file = getFileForRead(ts_filename);
            FILE *sax_file = getFileForRead(sax_filename);
            FILE *btree_file = getFileForRead(btree_filename);
            SAXVecType *tmp_sax_vec = new SAXVecType[btree_read_sax_batch_size];

            EXACT_SEARCH_START
            // search graph
            auto found_leaf_nodes_ids = approximate_graph->search(&search_ts_vec, approximate_search_node_num);


            auto to_search_ranges = getApproximateSearchRanges(found_leaf_nodes_ids);

            // search node from graph
            TopKHeap top_k_heap(k);
            PAAVecType search_paa_vec;
            types_helper::getPAAFromTS(search_ts_vec, search_paa_vec);

            APPROXIMATE_SEARCH_START
            for (int i = 0; i < to_search_ranges.size(); i++) {
                if (top_k_heap.sizeIsK()) {
                    approximateSearchNodes(search_ts_vec, search_paa_vec, to_search_ranges[i].first,
                                           to_search_ranges[i].second,
                                           tmp_sax_vec, top_k_heap, ts_file, sax_file);
                } else {
                    approximateSearchNodesNotK(search_ts_vec, to_search_ranges[i].first, to_search_ranges[i].second,
                                               top_k_heap, ts_file);
                }
            }

            // search all BTree
            // make sure heap.size is k
            if (!top_k_heap.sizeIsK()) {
                std::cout << "not k" << std::endl;
                exit(1);
            }

#if DIDS_PARALLELISM
            fclose(ts_file);
            fclose(sax_file);
            fclose(btree_file);
            delete[] tmp_sax_vec;

            APPROXIMATE_SEARCH_END


            const uint32_t thread_num = Const::thread_num;
            std::vector<TopKHeap> heaps(thread_num, top_k_heap);

            {
#pragma omp parallel for
                for (int ii = 0; ii < thread_num; ii++) {

                    std::vector<std::pair<uint64_t, uint64_t>> pos_pair_vec;
                    std::vector<uint32_t> ts_pos_vec;
                    auto tmp_leaf_node = new float[btree::BTree::btree_max_leaf_size];

                    FILE *btree_file_i = getFileForRead(btree_filename);
                    if (ii == thread_num - 1) {
                        for (int i = ii * (ref_objs_size / thread_num); i < ref_objs_size; i++) {
                            if (!btrees[i].all_size) continue;
                            float dis_b = sqrt(types_helper::getDisBetweenTSAndTS(search_ts_vec, ref_objs[i]));
                            auto pos_pair = btrees[i].search(dis_b, heaps[ii].top_sqrt, tmp_leaf_node, btree_file_i);
                            pos_pair_vec.push_back(pos_pair);
                        }
                    } else {
                        for (int i = ii * (ref_objs_size / thread_num);
                             i < ii * (ref_objs_size / thread_num) + ref_objs_size / thread_num; i++) {
                            if (!btrees[i].all_size) continue;
                            float dis_b = sqrt(types_helper::getDisBetweenTSAndTS(search_ts_vec, ref_objs[i]));
                            auto pos_pair = btrees[i].search(dis_b, heaps[ii].top_sqrt, tmp_leaf_node, btree_file_i);
                            pos_pair_vec.push_back(pos_pair);
                        }

                    }
                    fclose(btree_file_i);

                    FILE *sax_file_i = getFileForRead(sax_filename);
                    SAXVecType *tmp_sax_vec_i = new SAXVecType[btree_read_sax_batch_size];
                    for (auto &i: pos_pair_vec) {

                        auto pos_pairs = getRealRange(i, to_search_ranges);
                        for (int u = 0; u < pos_pairs.size(); u++) {
                            uint32_t to_search_size = pos_pairs[u].second - pos_pairs[u].first + 1;
                            uint64_t begin_pos = pos_pairs[u].first;
                            COUNT_READ_SAX_I(to_search_size, ii)

                            fseek(sax_file_i, begin_pos * sizeof(SAXVecType), SEEK_SET);


                            for (int kk = 0; kk < (to_search_size - 1) / btree_read_sax_batch_size + 1; kk++) {
                                uint64_t read_batch = btree_read_sax_batch_size;
                                if (to_search_size % btree_read_sax_batch_size &&
                                    kk == (to_search_size - 1) / btree_read_sax_batch_size) {
                                    read_batch = to_search_size % btree_read_sax_batch_size;
                                }
                                fread(tmp_sax_vec_i, sizeof(SAXVecType), read_batch, sax_file_i);
                                for (int j = 0; j < read_batch; j++) {
                                    float sax_dis = types_helper::getDisBetweenSAXAndPAA(tmp_sax_vec_i[j],
                                                                                         search_paa_vec);
                                    if (sax_dis >= heaps[ii].top) continue;
                                    ts_pos_vec.push_back(begin_pos + j);
                                }
                                begin_pos += read_batch;
                            }
                        }
                    }
                    fclose(sax_file_i);
                    delete[] tmp_sax_vec_i;


                    uint64_t search_ts_size = ts_pos_vec.size();
                    COUNT_READ_TS_I(search_ts_size, ii)

                    FILE *ts_file_i = getFileForRead(ts_filename);
                    TSVecType tmp_ts;
                    for (int i = 0; i < search_ts_size; i++) {
                        fseek(ts_file_i, ts_pos_vec[i] * sizeof(TSVecType), SEEK_SET);
                        fread(&tmp_ts, sizeof(TSVecType), 1, ts_file_i);
                        float dis = types_helper::getDisBetweenTSAndTS(search_ts_vec, tmp_ts);
                        if (dis >= heaps[ii].top) continue;
                        heaps[ii].push(dis, ts_pos_vec[i]);
                    }
                    fclose(ts_file_i);
                }
            }

#if DIDS_PRINT_INFO
            for (int i = 0; i < thread_num; i++) {
                COUNT_READ_SAX(read_sax_count_i[i])
                COUNT_READ_TS(read_ts_count_i[i])
            }
#endif

            std::vector<std::pair<float, uint64_t>> all_ans(k * thread_num);
            for (int i = 0; i < thread_num; i++) {
                auto ans_i = heaps[i].getRes();
                memcpy(&all_ans[i * k], ans_i.data(), sizeof(std::pair<float, uint64_t>) * k);
            }
            sort(all_ans.begin(), all_ans.end());
            all_ans.erase(unique(all_ans.begin(), all_ans.end()), all_ans.end());
            all_ans.resize(k);
            EXACT_SEARCH_END
            PRINT_INFO
            return all_ans;
#else

#if DIDS_HDD

            TSVecType tmp_ts;

            APPROXIMATE_SEARCH_END
            std::vector<std::pair<uint64_t, uint64_t>> pos_pair_vec;
            std::vector<uint32_t> ts_pos_vec;
            ts_pos_vec.reserve(100000);
            auto tmp_leaf_node = new float[btree::BTree::btree_max_leaf_size];

            for (int i = 0; i < ref_objs_size; i++) {
                if (!btrees[i].all_size) continue;
                float dis_b = sqrt(types_helper::getDisBetweenTSAndTS(search_ts_vec, ref_objs[i]));
                auto pos_pair = btrees[i].search(dis_b, top_k_heap.top_sqrt, tmp_leaf_node, btree_file);
                pos_pair_vec.push_back(pos_pair);
            }
            fclose(btree_file);

            for (auto & i : pos_pair_vec) {

                auto pos_pairs = getRealRange(i, to_search_ranges);
                for (int u=0;u<pos_pairs.size();u++) {
                    uint32_t to_search_size = pos_pairs[u].second - pos_pairs[u].first + 1;
                    uint64_t begin_pos = pos_pairs[u].first;
                    COUNT_READ_SAX(to_search_size)

                    fseek(sax_file, begin_pos * sizeof(SAXVecType), SEEK_SET);


                    for (int kk = 0; kk < (to_search_size - 1) / btree_read_sax_batch_size + 1; kk++) {
                        uint64_t read_batch = btree_read_sax_batch_size;
                        if (to_search_size % btree_read_sax_batch_size &&
                            kk == (to_search_size - 1) / btree_read_sax_batch_size) {
                            read_batch = to_search_size % btree_read_sax_batch_size;
                        }
                        fread(tmp_sax_vec, sizeof(SAXVecType), read_batch, sax_file);
                        for (int j = 0; j < read_batch; j++) {
                            float sax_dis = types_helper::getDisBetweenSAXAndPAA(tmp_sax_vec[j], search_paa_vec);
                            if (sax_dis >= top_k_heap.top) continue;
                            ts_pos_vec.push_back(begin_pos + j);
                        }
                        begin_pos += read_batch;
                    }
                }
            }

            fclose(sax_file);
            delete[] tmp_sax_vec;


            uint64_t search_ts_size = ts_pos_vec.size();
            COUNT_READ_TS(search_ts_size)


            for(int i=0;i<search_ts_size;i++) {
                fseek(ts_file, ts_pos_vec[i] * sizeof(TSVecType), SEEK_SET);
                fread(&tmp_ts, sizeof(TSVecType), 1, ts_file);
                float dis = types_helper::getDisBetweenTSAndTS(search_ts_vec, tmp_ts);
                if (dis >= top_k_heap.top) continue;
                top_k_heap.push(dis, ts_pos_vec[i]);
            }
            fclose(ts_file);

            EXACT_SEARCH_END
            PRINT_INFO
            return top_k_heap.getRes();

#else
            TSVecType tmp_ts;

            APPROXIMATE_SEARCH_END



            auto tmp_leaf_node = new float[btree::BTree::btree_max_leaf_size];
            for (int i = 0; i < ref_objs_size; i++) {
                if (!btrees[i].all_size) continue;
                float dis_b = sqrt(types_helper::getDisBetweenTSAndTS(search_ts_vec, ref_objs[i]));
                auto pos_pair = btrees[i].search(dis_b, top_k_heap.top_sqrt, tmp_leaf_node, btree_file);

                auto pos_pairs = getRealRange(pos_pair, to_search_ranges);
                for (int u=0;u<pos_pairs.size();u++) {
                    uint32_t to_search_size = pos_pairs[u].second - pos_pairs[u].first + 1;
                    uint64_t begin_pos = pos_pairs[u].first;
                    COUNT_READ_SAX(to_search_size)

                    fseek(sax_file, begin_pos * sizeof(SAXVecType), SEEK_SET);


                    for (int kk = 0; kk < (to_search_size - 1) / btree_read_sax_batch_size + 1; kk++) {
                        uint64_t read_batch = btree_read_sax_batch_size;
                        if (to_search_size % btree_read_sax_batch_size &&
                            kk == (to_search_size - 1) / btree_read_sax_batch_size) {
                            read_batch = to_search_size % btree_read_sax_batch_size;
                        }
                        fread(tmp_sax_vec, sizeof(SAXVecType), read_batch, sax_file);
                        for (int j = 0; j < read_batch; j++) {
                            float sax_dis = types_helper::getDisBetweenSAXAndPAA(tmp_sax_vec[j], search_paa_vec);
                            if (sax_dis >= top_k_heap.top) continue;
                            COUNT_READ_TS(1)
                            uint64_t read_ts_pos = begin_pos + j;
                            fseek(ts_file, read_ts_pos * sizeof(TSVecType), SEEK_SET);
                            fread(&tmp_ts, sizeof(TSVecType), 1, ts_file);
                            float dis = types_helper::getDisBetweenTSAndTS(search_ts_vec, tmp_ts);
                            if (dis >= top_k_heap.top) continue;

                            top_k_heap.push(dis, read_ts_pos);
                        }
                        begin_pos += read_batch;
                    }
                }
            }

            delete[] tmp_sax_vec;
            delete[] tmp_leaf_node;
            fclose(sax_file);
            fclose(btree_file);
            fclose(ts_file);
            EXACT_SEARCH_END
            PRINT_INFO
            return top_k_heap.getRes();

#endif
#endif
        }

        std::vector<std::pair<float, uint64_t>>
        approximateSearch(TSVecType &search_ts_vec, uint32_t k, uint32_t approximate_search_node_num) {
            FILE *ts_file = getFileForRead(ts_filename);
            FILE *sax_file = getFileForRead(sax_filename);
            SAXVecType *tmp_sax_vec = new SAXVecType[btree_read_sax_batch_size];

            // search graph
            auto found_leaf_nodes_ids = approximate_graph->search(&search_ts_vec, approximate_search_node_num);
            auto to_search_ranges = getApproximateSearchRanges(found_leaf_nodes_ids);

            // search node from graph
            TopKHeap top_k_heap(k);
            PAAVecType search_paa_vec;
            types_helper::getPAAFromTS(search_ts_vec, search_paa_vec);


            APPROXIMATE_SEARCH_START
            for (int i = 0; i < to_search_ranges.size(); i++) {
                if (top_k_heap.sizeIsK()) {
                    approximateSearchNodes(search_ts_vec, search_paa_vec, to_search_ranges[i].first,
                                           to_search_ranges[i].second,
                                           tmp_sax_vec, top_k_heap, ts_file, sax_file);
                } else {
                    approximateSearchNodesNotK(search_ts_vec, to_search_ranges[i].first, to_search_ranges[i].second,
                                               top_k_heap, ts_file);
                }
            }

            fclose(sax_file);
            fclose(ts_file);
            delete[] tmp_sax_vec;
            APPROXIMATE_SEARCH_END
            PRINT_INFO


            return top_k_heap.getRes();
        }

        std::vector<std::pair<uint64_t, uint64_t>>
        getRealRange(std::pair<uint64_t, uint64_t> pos_pair,
                     std::vector<std::pair<uint64_t, uint32_t>> &to_search_ranges) {
            std::vector<std::pair<uint64_t, uint64_t>> new_pos_pairs;
            if (pos_pair.second < pos_pair.first) return new_pos_pairs;
            for (auto &to_search_range: to_search_ranges) {
                uint64_t start = to_search_range.first;
                uint64_t end = to_search_range.first + to_search_range.second - 1;
                if (end < pos_pair.first) continue;
                if (start > pos_pair.second) break;

                if (start > pos_pair.first) {
                    new_pos_pairs.template emplace_back(pos_pair.first, start - 1);
                }
                pos_pair.first = end + 1;
                if (pos_pair.first > pos_pair.second) break;
            }
            if (pos_pair.first <= pos_pair.second) {
                new_pos_pairs.template emplace_back(pos_pair.first, pos_pair.second);
            }
            return new_pos_pairs;
        }

        std::vector<std::pair<uint64_t, uint32_t>>
        getApproximateSearchRanges(std::vector<std::pair<float, uint32_t>> &found_leaf_nodes_ids) {
            std::vector<std::pair<uint64_t, uint32_t>> to_search_ranges;

            uint64_t node_num = found_leaf_nodes_ids.size();
            std::vector<approximate_node_type> to_sort_nodes;
            to_sort_nodes.reserve(node_num);
            for (int i = 0; i < node_num; i++) {
                to_sort_nodes.push_back(approximate_leaf_nodes[found_leaf_nodes_ids[i].second]);
            }
            std::sort(to_sort_nodes.begin(), to_sort_nodes.end());
            uint64_t new_start = to_sort_nodes[0].file_pos;
            uint64_t new_start_end = to_sort_nodes[0].file_pos_end;
            for (int i = 1; i < node_num; i++) {
                if (to_sort_nodes[i].file_pos < new_start_end) {
                    new_start_end = std::max(new_start_end, to_sort_nodes[i].file_pos_end);
                } else {
                    to_search_ranges.emplace_back(new_start, new_start_end - new_start + 1);
                    new_start = to_sort_nodes[i].file_pos;
                    new_start_end = to_sort_nodes[i].file_pos_end;
                }
            }
            to_search_ranges.emplace_back(new_start, new_start_end - new_start + 1);
            return to_search_ranges;
        }

        void approximateSearchNodes(const TSVecType &search_ts_vec, const PAAVecType &search_paa_vec,
                                    const uint64_t read_start,
                                    const uint32_t read_size, SAXVecType *tmp_sax_vec,
                                    TopKHeap &top_k_heap,
                                    FILE *ts_file, FILE *sax_file) {

            COUNT_READ_SAX(read_size)
            uint64_t begin_pos = read_start;
            TSVecType tmp_ts;
            fseek(sax_file, read_start * sizeof(SAXVecType), SEEK_SET);
            for (int kk = 0; kk < (read_size - 1) / btree_read_sax_batch_size + 1; kk++) {
                uint64_t read_batch = btree_read_sax_batch_size;
                if (read_size % btree_read_sax_batch_size &&
                    kk == (read_size - 1) / btree_read_sax_batch_size) {
                    read_batch = read_size % btree_read_sax_batch_size;
                }
                fread(tmp_sax_vec, sizeof(SAXVecType), read_batch, sax_file);
                for (int j = 0; j < read_batch; j++) {
                    float sax_dis = types_helper::getDisBetweenSAXAndPAA(tmp_sax_vec[j], search_paa_vec);
                    if (sax_dis >= top_k_heap.top) continue;
                    COUNT_READ_TS(1)
                    uint64_t read_ts_pos = begin_pos + j;
                    fseek(ts_file, read_ts_pos * sizeof(TSVecType), SEEK_SET);
                    fread(&tmp_ts, sizeof(TSVecType), 1, ts_file);
                    float dis = types_helper::getDisBetweenTSAndTS(search_ts_vec, tmp_ts);
                    if (dis >= top_k_heap.top) continue;

                    top_k_heap.push(dis, read_ts_pos);

                }
                begin_pos += read_batch;
            }
        }

        void approximateSearchNodesNotK(const TSVecType &search_ts_vec, const uint64_t read_start,
                                        const uint32_t read_size,
                                        TopKHeap &top_k_heap,
                                        FILE *ts_file) {
            COUNT_READ_TS(read_size)
            fseek(ts_file, read_start * sizeof(TSVecType), SEEK_SET);
            TSVecType tmp_ts;
            for (int i = 0; i < read_size; i++) {
                fread(&tmp_ts, sizeof(TSVecType), 1, ts_file);
                float dis = types_helper::getDisBetweenTSAndTS(search_ts_vec, tmp_ts);
                top_k_heap.pushNotKnowK(dis, read_start + i);
            }
        }

    private:

        void loadRefObjs() {
            FILE *ref_objs_file = getFileForRead(ref_objs_filename);
            fread(&total_ts_num, sizeof(uint64_t), 1, ref_objs_file);
            fread(&ref_objs_size, sizeof(uint32_t), 1, ref_objs_file);
            ref_objs.resize(ref_objs_size);
            btrees.resize(ref_objs_size);
            fread(ref_objs.data(), sizeof(TSVecType), ref_objs_size, ref_objs_file);
            fclose(ref_objs_file);
        }

        void loadBtrees() {
            FILE *btree_leaf_nodes_file = getFileForRead(btree_leaf_nodes_filename);
            uint64_t file_pos = 0;
            for (int i = 0; i < ref_objs_size; i++) {
                uint64_t all_size;
                fread(&btrees[i].all_size, sizeof(uint64_t), 1, btree_leaf_nodes_file);
                all_size = btrees[i].all_size;
                if (!all_size) {
                    continue;
                }
                fread(&btrees[i].leaf_nums, sizeof(uint64_t), 1, btree_leaf_nodes_file);
                fread(&btrees[i].radius, sizeof(float), 1, btree_leaf_nodes_file);
                btrees[i].build(btree_leaf_nodes_file);
                btrees[i].btree_file_pos = file_pos;
                file_pos += all_size;
            }
            fclose(btree_leaf_nodes_file);
        }

        void loadApproximateStructs() {
            FILE *approximate_leaf_nodes_file = getFileForRead(approximate_leaf_nodes_filename);
            uint64_t approximate_leaf_nodes_size;
            fread(&approximate_leaf_nodes_size, sizeof(uint64_t), 1, approximate_leaf_nodes_file);
            approximate_leaf_nodes.resize(approximate_leaf_nodes_size);
            fread(approximate_leaf_nodes.data(), sizeof(approximate_node_type), approximate_leaf_nodes_size,
                  approximate_leaf_nodes_file);
            fclose(approximate_leaf_nodes_file);
            approximate_graph = GraphFactory<types_helper, TSVecType>::createTSVec(approximate_graph_filename);

        }


    private:

        void buildApproximateStructs(const uint64_t ts_num, const uint32_t approximate_leaf_size,
                                     std::vector<uint64_t> &btree_all_size_vec) {
            uint64_t max_approximate_leaf_num = ts_num / approximate_leaf_size * 2;
            approximate_graph = GraphFactory<types_helper, TSVecType>::createTSVec(max_approximate_leaf_num);


            uint64_t max_size = 0;
            for (int i = 0; i < ref_objs_size; i++) {
                max_size = std::max(btree_all_size_vec[i], max_size);
            }

            FILE *ts_file = getFileForWrite(ts_filename);
            fclose(ts_file);
            FILE *sax_file = getFileForWrite(sax_filename);
            fclose(sax_file);
            FILE *btree_file = getFileForWrite(btree_filename);
            FILE *btree_leaf_nodes_file = getFileForWrite(btree_leaf_nodes_filename);
            auto buffer_for_read = new TSVecType[max_size];

            auto btree_buffer_for_sort = new std::pair<float, uint32_t>[max_size];
            auto btree_buffer_for_write = new float[max_size];
            auto ts_buffer_for_write = new TSVecType[btree_write_ts_batch_size];
            auto sax_buffer_for_write = new SAXVecType[btree_write_ts_batch_size];

            uint64_t file_pos = 0;
            uint64_t leaf_id = 0;
            for (int i = 0; i < ref_objs_size; i++) {
                uint64_t all_size = btree_all_size_vec[i];
                if (!all_size) {
                    btrees[i].all_size = 0;
                    fwrite(&btrees[i].all_size, sizeof(uint64_t), 1, btree_leaf_nodes_file);
                    continue;
                }

                FILE *tmp_file = getFileForRead(tmp_filename + std::to_string(i));
                FILE *tmp_dis_file = getFileForRead(tmp_dis_filename + std::to_string(i));
                fread(buffer_for_read, sizeof(TSVecType), all_size, tmp_file);
                fread(btree_buffer_for_sort, sizeof(std::pair<float, uint32_t>), all_size, tmp_dis_file);
                fclose(tmp_file);
                fclose(tmp_dis_file);

                ts_file = getFileForWriteAppend(ts_filename);
                sax_file = getFileForWriteAppend(sax_filename);
                // do BTree
                // write ts and build BTree
                std::sort(btree_buffer_for_sort, btree_buffer_for_sort + all_size);
                uint64_t now_pos = 0;
                for (int k = 0; k < (all_size - 1) / btree_write_ts_batch_size + 1; k++) {
                    uint64_t read_batch = btree_write_ts_batch_size;
                    if (all_size % btree_write_ts_batch_size && k == (all_size - 1) / btree_write_ts_batch_size) {
                        read_batch = all_size % btree_write_ts_batch_size;
                    }
                    for (int j = 0; j < read_batch; j++) {
                        // write ts
                        ts_buffer_for_write[j] = buffer_for_read[btree_buffer_for_sort[now_pos].second];
                        types_helper::getSAXFromTS(ts_buffer_for_write[j], sax_buffer_for_write[j]);
                        // write BTree dis
                        btree_buffer_for_write[now_pos] = btree_buffer_for_sort[now_pos].first;
                        now_pos++;
                    }
                    fwrite(ts_buffer_for_write, sizeof(TSVecType), read_batch, ts_file);
                    fwrite(sax_buffer_for_write, sizeof(SAXVecType), read_batch, sax_file);
                }
                fclose(ts_file);
                fclose(sax_file);
                fwrite(btree_buffer_for_write, sizeof(float), all_size, btree_file);
                btrees[i].build(btree_buffer_for_write, all_size, file_pos);

                // save BTree leaf node and key
                fwrite(&btrees[i].all_size, sizeof(uint64_t), 1, btree_leaf_nodes_file);
                fwrite(&btrees[i].leaf_nums, sizeof(uint64_t), 1, btree_leaf_nodes_file);
                fwrite(&btrees[i].radius, sizeof(float), 1, btree_leaf_nodes_file);
                fwrite(btrees[i].leaf_nodes, sizeof(btree::LeafNode), btrees[i].leaf_nums, btree_leaf_nodes_file);
                fwrite(btrees[i].internal_keys, sizeof(btree::InternalKey), btrees[i].leaf_nums, btree_leaf_nodes_file);
                delete[] btrees[i].internal_keys;
                // read the saved ts
                ts_file = getFileForRead(ts_filename);
                fseek(ts_file, file_pos * sizeof(TSVecType), SEEK_SET);
                fread(buffer_for_read, sizeof(TSVecType), all_size, ts_file);
                fclose(ts_file);
                // build binary tree
                binary_tree_type *tmp_binary_tree = new binary_tree_type(approximate_leaf_size);
                for (int j = 0; j < all_size; j++) {
                    tmp_binary_tree->insert(j, buffer_for_read);
                }
                tmp_binary_tree->insertGraph(approximate_graph, leaf_id, buffer_for_read);

                tmp_binary_tree->write(file_pos, approximate_leaf_nodes);

                if (leaf_id + tmp_binary_tree->leaf_nodes.size() > max_approximate_leaf_num) {
                    std::cout << "graph leaf num " << max_approximate_leaf_num << " not enough " << leaf_id
                              << std::endl;
                    exit(1);
                }
                delete tmp_binary_tree;
                file_pos += all_size;

                remove((tmp_dis_filename + std::to_string(i)).c_str());
                remove((tmp_filename + std::to_string(i)).c_str());

                fprintf(stderr, "\r\x1b[m>> build: \x1b[36m%2.2lf%%\x1b[0m",
                        (float) ((float) (i + 1) / ref_objs_size) * 100);
            }
            fprintf(stderr, "\r\x1b[m>> build: \x1b[36m%2.2lf%%\x1b[0m", (float) 100);
            fprintf(stderr, "\n");


            delete[] buffer_for_read;
            delete[] btree_buffer_for_sort;
            delete[] btree_buffer_for_write;
            delete[] ts_buffer_for_write;
            fclose(btree_file);
            fclose(btree_leaf_nodes_file);

            // save approximate nodes and graph
            FILE *approximate_leaf_nodes_file = getFileForWrite(approximate_leaf_nodes_filename);
            fwrite(&leaf_id, sizeof(uint64_t), 1, approximate_leaf_nodes_file);
            fwrite(approximate_leaf_nodes.data(), sizeof(approximate_node_type), leaf_id, approximate_leaf_nodes_file);
            fclose(approximate_leaf_nodes_file);

            approximate_graph->saveIndex(approximate_graph_filename);

        }


        static void run(TSVecType *ts_buffer_for_read, Graph *ref_objs_graph, uint64_t begin, uint64_t end,
                        std::vector<std::pair<float, uint32_t>> &ts_ref_obj) {
            for (uint64_t i = begin; i < end; i++) {
                ts_ref_obj[i] = (ref_objs_graph->search(&ts_buffer_for_read[i], 1))[0];
            }
        }


        std::vector<uint64_t> groupTSByRefObjs(const uint64_t ts_num, const uint64_t ts_buffer_size_for_read,
                                               const uint64_t ts_buffer_size_for_one_ref_obj,
                                               const std::string &input_filename) {

            std::vector<uint64_t> btree_all_size_vec(ref_objs_size, 0);

            // build the ts graph
            auto ref_objs_graph = GraphFactory<types_helper, TSVecType>::createTSVec(ref_objs_size);
            for (int i = 0; i < ref_objs_size; i++) {
                ref_objs_graph->insert(&ref_objs[i], i);
            }

            // malloc the read and write buffer
            auto ts_buffer_for_read = new TSVecType[ts_buffer_size_for_read];
            std::vector<std::pair<TSVecType *, uint64_t>> ts_buffer_for_write(ref_objs_size);
            std::vector<std::pair<float, uint32_t> *> ts_dis_buffer_for_write(ref_objs_size);
            auto tmp_p = new TSVecType[ts_buffer_size_for_one_ref_obj * ref_objs_size];
            auto tmp_dis_p = new std::pair<float, uint32_t>[ts_buffer_size_for_one_ref_obj * ref_objs_size];

            for (int i = 0; i < ref_objs_size; i++) {
                ts_buffer_for_write[i] = {tmp_p + ts_buffer_size_for_one_ref_obj * i, 0};
                ts_dis_buffer_for_write[i] = tmp_dis_p + ts_buffer_size_for_one_ref_obj * i;
                // clean all files
                FILE *tmp_file = getFileForWrite(tmp_filename + std::to_string(i));
                FILE *tmp_dis_file = getFileForWrite(tmp_dis_filename + std::to_string(i));
                fclose(tmp_file);
                fclose(tmp_dis_file);
            }

//            std::cout << "read ts and group, write the full buffer" << std::endl;
            // read ts and group, write the full buffer
            FILE *input_file = getFileForRead(input_filename);

            uint64_t now_pos = 0;
            for (int i = 0; i < (ts_num - 1) / ts_buffer_size_for_read + 1; i++) {
                uint64_t read_batch = ts_buffer_size_for_read;
                if (ts_num % ts_buffer_size_for_read && i == (ts_num - 1) / ts_buffer_size_for_read) {
                    read_batch = ts_num % ts_buffer_size_for_read;
                }
                fread(ts_buffer_for_read, sizeof(TSVecType), read_batch, input_file);

#if DIDS_PARALLELISM
                std::vector<std::pair<float, uint32_t>> ts_ref_obj(read_batch);
                const uint64_t thread_num = Const::thread_num;
                std::thread thread_vec[thread_num];

                uint64_t batch_thread = read_batch / thread_num;
                for (int j = 0; j < thread_num - 1; j++) {
                    thread_vec[j] = std::thread(run, ts_buffer_for_read, ref_objs_graph, j * batch_thread,
                                                j * batch_thread + batch_thread, std::ref(ts_ref_obj));
                }
                thread_vec[thread_num - 1] = std::thread(run, ts_buffer_for_read, ref_objs_graph,
                                                         (thread_num - 1) * batch_thread, read_batch,
                                                         std::ref(ts_ref_obj));
                for (int j = 0; j < thread_num; j++) {
                    thread_vec[j].join();
                }

                for (int j = 0; j < read_batch; j++) {
                    // find the ref obj by graph
                    auto ts_ref_obj_id = ts_ref_obj[j].second;
                    auto ts_dis = std::sqrt(ts_ref_obj[j].first);
                    auto &now_pos_in_buffer = ts_buffer_for_write[ts_ref_obj_id].second;
                    ts_buffer_for_write[ts_ref_obj_id].first[now_pos_in_buffer] = ts_buffer_for_read[j];
                    ts_dis_buffer_for_write[ts_ref_obj_id][now_pos_in_buffer].first = ts_dis;
                    ts_dis_buffer_for_write[ts_ref_obj_id][now_pos_in_buffer].second = btree_all_size_vec[ts_ref_obj_id]++;
                    now_pos_in_buffer++;
                    if (now_pos_in_buffer == ts_buffer_size_for_one_ref_obj) {
                        // full append write
                        FILE *tmp_file = getFileForWriteAppend(tmp_filename + std::to_string(ts_ref_obj_id));
                        FILE *tmp_dis_file = getFileForWriteAppend(tmp_dis_filename + std::to_string(ts_ref_obj_id));
                        fwrite(ts_buffer_for_write[ts_ref_obj_id].first, sizeof(TSVecType),
                               ts_buffer_size_for_one_ref_obj, tmp_file);
                        fwrite(ts_dis_buffer_for_write[ts_ref_obj_id], sizeof(std::pair<float, uint32_t>),
                               ts_buffer_size_for_one_ref_obj, tmp_dis_file);
                        fclose(tmp_file);
                        fclose(tmp_dis_file);
                        now_pos_in_buffer = 0;
                    }
                }
#else
                for (int j = 0; j < read_batch; j++) {
                    // find the ref obj by graph
                    auto ts_ref_obj = ref_objs_graph->search(&ts_buffer_for_read[j], 1);
                    auto ts_ref_obj_id = ts_ref_obj[0].second;
                    auto ts_dis = std::sqrt(ts_ref_obj[0].first);
                    auto &now_pos_in_buffer = ts_buffer_for_write[ts_ref_obj_id].second;
                    ts_buffer_for_write[ts_ref_obj_id].first[now_pos_in_buffer] = ts_buffer_for_read[j];
                    ts_dis_buffer_for_write[ts_ref_obj_id][now_pos_in_buffer].first = ts_dis;
                    ts_dis_buffer_for_write[ts_ref_obj_id][now_pos_in_buffer].second = btree_all_size_vec[ts_ref_obj_id]++;
                    now_pos_in_buffer++;
                    if (now_pos_in_buffer == ts_buffer_size_for_one_ref_obj) {
                        // full append write
                        FILE *tmp_file = getFileForWriteAppend(tmp_filename + std::to_string(ts_ref_obj_id));
                        FILE *tmp_dis_file = getFileForWriteAppend(tmp_dis_filename + std::to_string(ts_ref_obj_id));
                        fwrite(ts_buffer_for_write[ts_ref_obj_id].first, sizeof(TSVecType),
                               ts_buffer_size_for_one_ref_obj, tmp_file);
                        fwrite(ts_dis_buffer_for_write[ts_ref_obj_id], sizeof(std::pair<float, uint32_t>),
                               ts_buffer_size_for_one_ref_obj, tmp_dis_file);
                        fclose(tmp_file);
                        fclose(tmp_dis_file);
                        now_pos_in_buffer = 0;
                    }
                }
#endif
                now_pos += read_batch;
                fprintf(stderr, "\r\x1b[m>> group: \x1b[36m%2.2lf%%\x1b[0m",
                        (float) ((float) now_pos / total_ts_num) * 100);
            }

            // do the rest

            for (int i = 0; i < ref_objs_size; i++) {
                if (ts_buffer_for_write[i].second != 0) {
                    FILE *tmp_file = getFileForWriteAppend(tmp_filename + std::to_string(i));
                    FILE *tmp_dis_file = getFileForWriteAppend(tmp_dis_filename + std::to_string(i));
                    fwrite(ts_buffer_for_write[i].first, sizeof(TSVecType), ts_buffer_for_write[i].second, tmp_file);
                    fwrite(ts_dis_buffer_for_write[i], sizeof(std::pair<float, uint32_t>),
                           ts_buffer_for_write[i].second, tmp_dis_file);
                    fclose(tmp_file);
                    fclose(tmp_dis_file);
                }
            }
            fprintf(stderr, "\r\x1b[m>> group: \x1b[36m%2.2lf%%\x1b[0m", (float) 100.);
            fprintf(stderr, "\n");

            fclose(input_file);
            delete[] ts_buffer_for_read;
            delete[] tmp_p;
            delete[] tmp_dis_p;
            delete ref_objs_graph;

            return btree_all_size_vec;
        }


        // random select
        void selectRefObjs(const uint64_t ts_num, const std::string &input_filename) {

            FILE *input_file = getFileForRead(input_filename);

            const int interval = 100;
            int select_ts_num = ts_num / interval;
            std::vector<std::array<float, TSVecType::ts_length>> select_ts_vec(select_ts_num);
            for (int i = 0; i < select_ts_num; i++) {
                fseek(input_file, i * interval * sizeof(TSVecType), SEEK_SET);
                fread(&select_ts_vec[i], sizeof(TSVecType), 1, input_file);
            }
//            std::cout << "kmeans start" << std::endl;
#if DIDS_PARALLELISM
            auto k_res = dkm::kmeans_lloyd_parallel(select_ts_vec, ref_objs_size, 30);
#else
            auto k_res = dkm::kmeans_lloyd<GraphFactory<types_helper, TSVecType>, float, TSVecType::ts_length>(select_ts_vec, ref_objs_size, 30);
#endif
//            std::cout << "kmeans over" << std::endl;
            auto k_res_vec = std::move(std::get<0>(k_res));
            memcpy(ref_objs.data(), k_res_vec.data(), sizeof(TSVecType) * ref_objs_size);

            fclose(input_file);
            FILE *ref_objs_file = getFileForWrite(ref_objs_filename);
            fwrite(&total_ts_num, sizeof(uint64_t), 1, ref_objs_file);
            fwrite(&ref_objs_size, sizeof(uint32_t), 1, ref_objs_file);
            fwrite(ref_objs.data(), sizeof(TSVecType), ref_objs_size, ref_objs_file);
            fclose(ref_objs_file);

        }

    private:


        inline FILE *getFileForRead(const std::string &filename) {
            FILE *file_ = fopen(filename.c_str(), "rb");
            if (!file_) {
                std::cout << "cannot find " << filename << std::endl;
                exit(1);
            }
            return file_;
        }

        inline FILE *getFileForWrite(const std::string &filename) {
            FILE *file_ = fopen(filename.c_str(), "wb");
            if (!file_) {
                std::cout << "cannot find " << filename << std::endl;
                exit(1);
            }
            return file_;
        }

        inline FILE *getFileForWriteAppend(const std::string &filename) {
            FILE *file_ = fopen(filename.c_str(), "ab");
            if (!file_) {
                std::cout << "cannot find " << filename << std::endl;
                exit(1);
            }
            return file_;
        }


        static constexpr uint32_t btree_write_ts_batch_size = Const::btree_write_ts_batch_size;
        static constexpr uint32_t btree_read_sax_batch_size = Const::btree_read_sax_batch_size;

        const std::string ts_filename;
        const std::string sax_filename;
        // float
        const std::string btree_filename;
        // leaf node and key
        const std::string btree_leaf_nodes_filename;
        const std::string ref_objs_filename;
        const std::string approximate_ts_pos_filename;
        const std::string approximate_leaf_nodes_filename;
        const std::string approximate_graph_filename;
        const std::string tmp_filename;
        const std::string tmp_dis_filename;

        uint64_t total_ts_num;
        uint32_t ref_objs_size;
        std::vector<TSVecType> ref_objs;


        std::vector<btree::BTree> btrees;

        // for approximate search
        std::vector<approximate_node_type> approximate_leaf_nodes;
        Graph *approximate_graph;

#if DIDS_PRINT_INFO
        timeval exact_search_time_start;
        timeval approximate_search_time_start;
        timeval current_time;
        double start_time;
        double end_time;
        double exact_search_time;
        double approximate_search_time;

        uint64_t read_ts_count;
        uint64_t read_sax_count;
#if DIDS_PARALLELISM
        uint64_t read_ts_count_i[Const::thread_num];
        uint64_t read_sax_count_i[Const::thread_num];
#endif
#endif

    };


}


#endif
