/* ms_pointers - Computes the matching statistics pointers from BWT and Thresholds
    Copyright (C) 2020 Massimiliano Rossi

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*!
   \file ms_pointers.hpp
   \brief ms_pointers.hpp Computes the matching statistics pointers from BWT and Thresholds.
   \author Massimiliano Rossi
   \date 09/07/2020
*/

#ifndef _AUG_MS_POINTERS_HH
#define _AUG_MS_POINTERS_HH


/** FLAGS **/
// #define MEASURE_TIME 1  //measure the time for LCE and backward search?
//#define NAIVE_LCE_SCHEDULE 1 //stupidly execute two LCEs without heurstics
#include "Common.hpp"
#include <cstddef>
#ifndef NAIVE_LCE_SCHEDULE //apply a heuristic
#define SORT_BY_DISTANCE_HEURISTIC 1 // apply a heuristic to compute the LCE with the closer BWT position first
#endif

#ifndef DCHECK_HPP
#define DCHECK_HPP
#include <string>
#include <sstream>
#include <stdexcept>

#ifndef DCHECK
#ifdef NDEBUG
#define ON_DEBUG(x)
#define DCHECK_(x,y,z)
#define DCHECK(x)
#define DCHECK_EQ(x, y)
#define DCHECK_NE(x, y)
#define DCHECK_LE(x, y)
#define DCHECK_LT(x, y)
#define DCHECK_GE(x, y)
#define DCHECK_GT(x, y)
#else//NDEBUG
#define ON_DEBUG(x) x
#define DCHECK_(x,y,z) \
  if (!(x)) throw std::runtime_error(std::string(" in file ") + __FILE__ + ':' + std::to_string(__LINE__) + (" the check failed: " #x) + ", we got " + std::to_string(y) + " vs " + std::to_string(z))
#define DCHECK(x) \
  if (!(x)) throw std::runtime_error(std::string(" in file ") + __FILE__ + ':' + std::to_string(__LINE__) + (" the check failed: " #x))
#define DCHECK_EQ(x, y) DCHECK_((x) == (y), x,y)
#define DCHECK_NE(x, y) DCHECK_((x) != (y), x,y)
#define DCHECK_LE(x, y) DCHECK_((x) <= (y), x,y)
#define DCHECK_LT(x, y) DCHECK_((x) < (y) ,x,y)
#define DCHECK_GE(x, y) DCHECK_((x) >= (y),x,y )
#define DCHECK_GT(x, y) DCHECK_((x) > (y) ,x,y)
#endif //NDEBUG
#endif //DCHECK
#endif /* DCHECK_HPP */


#include <common.hpp>

#include <malloc_count.h>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>

#include <r_index.hpp>

#include <ms_rle_string.hpp>
#include <thresholds_ds.hpp>
#include <thr_lce_ds.hpp>

#include "PlainSlp.hpp"
#include "PoSlp.hpp"
#include "ShapedSlp_Status.hpp"
#include "ShapedSlp.hpp"
#include "ShapedSlpV2.hpp"
#include "SelfShapedSlp.hpp"
#include "SelfShapedSlpV2.hpp"
#include "DirectAccessibleGammaCode.hpp"
#include "IncBitLenCode.hpp"
#include "FixedBitLenCode.hpp"
#include "SelectType.hpp"
#include "VlcVec.hpp"

#ifdef MEASURE_TIME
struct Stopwatch {
    std::chrono::high_resolution_clock::time_point m_val = std::chrono::high_resolution_clock::now();

    void reset()
    {
        m_val = std::chrono::high_resolution_clock::now();
    }
    double seconds() const
    {
        return std::chrono::duration<double, std::ratio<1>>(std::chrono::high_resolution_clock::now() - m_val).count();
    }

};
#endif//MEASURE_TIME

using var_t = uint32_t;
using Fblc = FixedBitLenCode<>;
using SelSd = SelectSdvec<>;
using SelMcl = SelectMcl<>;
using DagcSd = DirectAccessibleGammaCode<SelSd>;
using DagcMcl = DirectAccessibleGammaCode<SelMcl>;
using Vlc64 = VlcVec<sdsl::coder::elias_delta, 64>;
using Vlc128 = VlcVec<sdsl::coder::elias_delta, 128>;


template <
    class SlpT = PlainSlp<var_t, Fblc, Fblc>,
    class rle_string_t = ms_rle_string_sd,
    class thr_lce_t = thr_lce_plain<rle_string_t>,
    class sparse_bv_type = ri::sparse_sd_vector,
    class thresholds_t = thr_bv<rle_string_t>
    >
class ms_pointers : ri::r_index<sparse_bv_type, rle_string_t>
{
public:
    thresholds_t thresholds;
    thr_lce_t thr_lce;

    SlpT slp;

    int_vector<> samples_start;

    typedef size_t size_type;

    ms_pointers()
        : ri::r_index<sparse_bv_type, rle_string_t>()
    {}

    void build(const std::string& filename, int bytes = 0)
    {
        verbose("Building the r-index from BWT");

        std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

        std::string bwt_fname = filename + ".bwt";

        verbose("RLE encoding BWT and computing SA samples");

        if (true) {
            std::string bwt_heads_fname = bwt_fname + ".heads";
            std::ifstream ifs_heads(bwt_heads_fname);
            std::string bwt_len_fname = bwt_fname + ".len";
            std::ifstream ifs_len(bwt_len_fname);
            this->bwt = rle_string_t(ifs_heads, ifs_len);

            ifs_heads.seekg(0);
            ifs_len.seekg(0);
            this->build_F_(ifs_heads, ifs_len);
        } else {
            std::ifstream ifs(bwt_fname);
            this->bwt = rle_string_t(ifs);

            ifs.seekg(0);
            this->build_F(ifs);
        }

        this->r = this->bwt.number_of_runs();
        ri::ulint n = this->bwt.size();
        int log_n = bitsize(uint64_t(this->bwt.size()));

        verbose("Number of BWT equal-letter runs: r = ", this->r);
        verbose("Rate n/r = ", double(this->bwt.size()) / this->r);
        verbose("log2(r) = ", log2(double(this->r)));
        verbose("log2(n/r) = ", log2(double(this->bwt.size()) / this->r));
        verbose("log2(n) = ", log_n);

        read_samples(filename + ".ssa", this->r, n, samples_start);
        read_samples(filename + ".esa", this->r, n, this->samples_last);

        std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();
        verbose("R-index construction complete");
        verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
        verbose(3);

        //load_grammar(filename);

        verbose("text length: ", slp.getLen());
        verbose("bwt length: ", this->bwt.size());

        thresholds = thresholds_t(filename, &this->bwt);
        thr_lce = thr_lce_t(filename, &this->bwt, bytes);

        verbose("finished augmented thresholds construction");
    }


    void load_grammar(const std::string& filename)
    {
        {
            verbose("Load Grammar");
            std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();
            if (std::is_same<SlpT, PlainSlp<var_t, Fblc, Fblc>>::value) {
                ifstream fs(filename + ".plain.slp");
                slp.load(fs);
                fs.close();
            } else {
                ifstream fs(filename + ".slp");
                slp.load(fs);
                fs.close();
            }

            std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();
            verbose("Memory peak: ", malloc_count_peak());
            verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
        }
        DCHECK_EQ(slp.getLen()+1, this->bwt.size());
    }



    void read_samples(std::string filename, ulint r, ulint n, int_vector<> &samples)
    {
        int log_n = bitsize(uint64_t(n));
        struct stat filestat;
        FILE *fd;

        if ((fd = fopen(filename.c_str(), "r")) == nullptr)
            error("open() file " + filename + " failed");

        int fn = fileno(fd);
        if (fstat(fn, &filestat) < 0)
            error("stat() file " + filename + " failed");

        if (filestat.st_size % SSABYTES != 0)
            error("invilid file " + filename);

        size_t length = filestat.st_size / (2*SSABYTES);
        //Check that the length of the file is 2*r elements of 5 bytes
        assert(length == r);

        // Create the vector
        samples = int_vector<>(r, 0, log_n);

        // Read the vector
        uint64_t left = 0;
        uint64_t right = 0;
        size_t i = 0;
        while (fread((char *)&left, SSABYTES, 1, fd) && fread((char *)&right, SSABYTES, 1,fd)) {
            ulint val = (right ? right - 1 : n - 1);
            assert(bitsize(uint64_t(val)) <= log_n);
            samples[i++] = val;
        }

        fclose(fd);
    }

    void write_int(ostream& os, const size_t& i)
    {
        os.write(reinterpret_cast<const char*>(&i), sizeof(size_t));
    }

    // Computes the matching statistics pointers for the given pattern
    size_t query(const std::string& patternfile, const std::string& len_filename, const std::string& ref_filename) {
        const char* p;
        size_t m;
        map_file(patternfile.c_str(), p, m);

        auto pattern_at = [&] (const size_t pos) {
            return p[pos];
        };

        ofstream len_file(len_filename, std::ios::binary);
        ofstream ref_file(ref_filename, std::ios::binary);

        const size_t n = slp.getLen();
        verbose("pattern length: ", m);

        size_t last_len;
        size_t last_ref;

        size_t lce_freq = 3, lce_cnt = 0; // freq of LCE skipping and a counter for it
        vector<size_t> stored_sample_pos(lce_freq), stored_ptr(lce_freq);
        vector<int> stored_it(lce_freq);
        vector<bool> direction(lce_freq);
        bool lce_is_paused = false;

        auto write_len = [&] (const size_t i, bool lce_is_paused=false) {
            if (lce_is_paused) return;
            last_len = i;
            write_int(len_file, last_len);
        };
        auto write_ref = [&] (const size_t i) {
            last_ref = i;
            write_int(ref_file, last_ref);
        };

        write_len(1, lce_is_paused);

        //! Start with the last character Adrian: do we assume the last char occurs in T??
        auto pos = this->bwt.select(1, pattern_at(m-1));
        {
            const ri::ulint run_of_j = this->bwt.run_of_position(pos);
            write_ref(samples_start[run_of_j]);
        }
        pos = LF(pos, pattern_at(m-1));

#ifdef MEASURE_TIME
        double time_lce = 0;
        double time_backwardstep = 0;
        size_t count_lce_total = 0;
        size_t count_lce_skips = 0;
#endif

        for (size_t i = 1; i < m; ++i) {
            verbose("i= ", i);

            const auto c = pattern_at(m - i - 1);
            const size_t number_of_runs_of_c = this->bwt.number_of_letter(c);
            if (number_of_runs_of_c == 0) {
                write_len(0, lce_is_paused);
                write_ref(1); // why 1???
            } else if (pos < this->bwt.size() && this->bwt[pos] == c) {
                write_len(last_len+1, lce_is_paused);
                write_ref(last_ref-1);
            } else { // we jump
                const ri::ulint rank = this->bwt.rank(pos, c);

                ri::ulint run0, run1;
                size_t sa0, sa1;

                if (rank < number_of_runs_of_c) {
                    sa1 = this->bwt.select(rank, c);
                }
                if (rank >= number_of_runs_of_c) {
                    sa0 = this->bwt.select(rank-1, c);
                }

                struct Triplet {
                    size_t sa, ref, len;
                };

                auto delay_succeeding_lce = [&] () -> Triplet {
                    const size_t textposStart = this->samples_start[run1];

                    stored_it[lce_cnt] = i;
                    stored_sample_pos[lce_cnt] = textposStart;
                    stored_ptr[lce_cnt] = last_ref;
                    direction[lce_cnt] = 1;
                    lce_cnt++;
                    return {sa1, textposStart, 0};

                };

                auto delay_preceding_lce = [&] () -> Triplet {
                    const ri::ulint run0 = this->bwt.run_of_position(sa0);
                    const size_t textposLast = this->samples_last[run0];

                    stored_it[lce_cnt] = i;
                    stored_sample_pos[lce_cnt] = textposLast;
                    stored_ptr[lce_cnt] = last_ref;
                    direction[lce_cnt] = 0;
                    lce_cnt++;
                    return {sa0, textposLast, 0};
                };

                auto compute_succeeding_lce = [&] (const size_t textposStart, size_t last_ref, size_t last_len) -> size_t {
                    const size_t lenStart = textposStart+1 >= n ? 0 : lceToRBounded(slp, textposStart+1, last_ref, last_len);
                    const size_t len1 = std::min(last_len, lenStart);
                    return len1;
                };

                auto compute_preceding_lce = [&] (const size_t textposLast, size_t last_ref, size_t last_len) -> size_t {
                    const size_t lenLast = textposLast+1 >= n ? 0 : lceToRBounded(slp, textposLast+1, last_ref, last_len);
                    const size_t len0 = std::min(last_len, lenLast);
                    return len0;
                };

                const Triplet jump = [&] () -> Triplet {
                    if(rank == 0) {
                        // Check only succeeding -> we ignore thresholds in this case
                        return delay_succeeding_lce();
                    } else if(rank >= number_of_runs_of_c) {
                        // Check only preceding -> we ignore thresholds in this case
                        return delay_preceding_lce();
                    }
                    // Check thresholds and boundary LCEs first
                    run1 = this->bwt.run_of_position(sa1);
                    const size_t thr = thresholds[run1];
                    if (pos < thr) {
                        if (!lce_is_paused && thr_lce.skip_preceding_lce(run1, last_len)) {
                            const ri::ulint run0 = this->bwt.run_of_position(sa0);
                            const size_t ref0 = this->samples_last[run0];
                            const size_t len0 = last_len;
                            return {sa0, ref0, len0};
                        } else {
                            return delay_preceding_lce();
                        }
                    } else {
                        if (!lce_is_paused && thr_lce.skip_succeeding_lce(run1, last_len)) {
                            const size_t ref1 = this->samples_start[run1];;
                            const size_t len1 = last_len;
                            return {sa1, ref1, len1};
                        } else {
                            return delay_succeeding_lce();
                        }
                    }
                }();

                if (lce_cnt == lce_freq) {
                    lce_is_paused = false;

                    size_t lce = 0;
                    if (direction[lce_freq-1] == 0)     lce = compute_preceding_lce(stored_sample_pos[lce_freq-1],last_ref, last_len);
                    else                                lce = compute_succeeding_lce(stored_sample_pos[lce_freq-1], last_ref, last_len);

                    int skipped_steps = (stored_it[lce_freq-1] - stored_it[0]);
                    // if we still process the same MEM, we know the lens right ahead
                    if ((last_len + skipped_steps) == lce) {
                        for (int j = 0; j < skipped_steps; j++) write_len(last_len+1, lce_is_paused); 
                    } else {
                        // do the delayed LCEs
                        for (int j = 0; j < lce_freq; j++) {
                            // TODO: try to use the stored LCE's
                            if (direction[j] == 0)  last_len = compute_preceding_lce(stored_sample_pos[j],last_ref, last_len);
                            else                    last_len = compute_succeeding_lce(stored_sample_pos[j], last_ref, last_len);

                            if (j < (lce_freq - 1)) {
                                for (int k = 0; k < (stored_it[j + 1] - stored_it[j]); k++) {
                                    write_len(last_len+1, lce_is_paused);
                                }
                            }
                        }
                        // the last LCE wasn't written yet
                        write_len(last_len, lce_is_paused);
                        lce_cnt = 0;
                    }
                } else {
                    write_len(1 + jump.len, lce_is_paused);
                }
                write_ref(jump.ref);
                pos = jump.sa;
            }
            pos = LF(pos, c); //! Perform one backward step
        }
#ifdef MEASURE_TIME
        cout << "Time backwardsearch: " << time_backwardstep << std::endl;
        cout << "Time lce: " << time_lce << std::endl;
        cout << "Count lce: " << count_lce_total << std::endl;
        cout << "Count lce skips: " << count_lce_skips << std::endl;
#endif
        return m;
    }
    // Computes the matching statistics pointers for the given pattern
    std::pair<std::vector<size_t>, std::vector<size_t>> query(const std::vector<uint8_t> &pattern)
    {
        size_t m = pattern.size();

        return _query(pattern.data(), m);
    }

    std::pair<std::vector<size_t>, std::vector<size_t>> query(const char *pattern, const size_t m, const size_t lcp_delay)
    {
        return _query(pattern, m, lcp_delay);
    }

    // Computes the matching statistics pointers for the given pattern
    template <typename string_t>
    std::pair<std::vector<size_t>, std::vector<size_t>> _query(const string_t &pattern, const size_t m, const size_t lce_freq)
    {
        auto pattern_at = [&] (const size_t pos) {
            return pattern[pos];
        };

        std::vector<size_t> lens;
        std::vector<size_t> refs;

        lens.reserve(m);
        refs.reserve(m);

        const size_t n = slp.getLen();
        verbose("Pattern length: ", m);
        slp.precompute_pattern(pattern); // KR hashes of the pattern prefixes

        size_t last_len;
        size_t last_ref;

        size_t lce_cnt = 0;
        vector<size_t> stored_sample_pos(lce_freq), stored_ptr(lce_freq), stored_run(lce_freq);
        vector<int> stored_it(lce_freq+1, 0);
        vector<bool> direction(lce_freq);
        bool lce_is_paused = false;

        auto write_len = [&] (const size_t l, bool lce_is_paused) {
            if (lce_is_paused) return;
            last_len = l;
            lens.push_back(last_len);
        };
        auto write_ref = [&] (const size_t p) {
            //verbose("p= ", p);
            last_ref = p;
            refs.push_back(last_ref);
        };

        auto write_len_segment = [&] (const int n_iters) {
            for (int j = 0; j < n_iters; j++) write_len(last_len+1, lce_is_paused);
        };

        write_len(1, lce_is_paused);

        //! Start with the last character
        auto pos = this->bwt.select(1, pattern_at(m-1));
        {
            const ri::ulint run_of_j = this->bwt.run_of_position(pos);
            write_ref(samples_start[run_of_j]);
        }
        pos = LF(pos, pattern_at(m-1));

        struct Triplet {
            size_t sa, ref, len;
        };
        auto delay_preceding_lce = [&] (const size_t run1, const size_t rank, char c, int i) -> Triplet {
            const size_t sa0 = this->bwt.select(rank-1, c);
            const ri::ulint run0 = this->bwt.run_of_position(sa0);
            const size_t textposLast = this->samples_last[run0];
            //verbose("i = ", i, "last_ref = ", last_ref, " lce_is_paused =", lce_is_paused, "UP= ", textposLast);
            lce_is_paused = true;
            stored_it[lce_cnt] = i;
            stored_sample_pos[lce_cnt] = textposLast;
            stored_ptr[lce_cnt] = last_ref;
            stored_run[lce_cnt] = run1;
            direction[lce_cnt] = 0;
            lce_cnt++;
            return {sa0, textposLast, 0};
        };

        auto delay_succeeding_lce = [&] (const size_t sa, const size_t run, int i) -> Triplet {
            const size_t textposStart = this->samples_start[run];
            //verbose("i = ", i, "last_ref = ", last_ref, " lce_is_paused =", lce_is_paused, "DOWN= ", textposStart, " run = ", this->bwt.run_of_position(sa1));
            lce_is_paused = true;
            stored_it[lce_cnt] = i;
            stored_sample_pos[lce_cnt] = textposStart;
            stored_ptr[lce_cnt] = last_ref;
            stored_run[lce_cnt] = run;
            direction[lce_cnt] = 1;
            lce_cnt++;
            return {sa, textposStart, 0};

        };

        auto compute_lce = [&] (const size_t pos_sample, const size_t pos_pattern, const size_t max_len) -> size_t {
            const size_t len = ((pos_sample + 1) >= n) ? 0 : match_length_query(slp, pos_sample + 1, pos_pattern);
            const size_t res_len = std::min(max_len, len);
            return res_len;
        };

        auto empty_stack = [&] () {
            for (int j = 0; j < lce_cnt; j++) {
                if (((direction[j] == 0) && !thr_lce.skip_preceding_lce(stored_run[j], last_len)) ||
                    ((direction[j] == 1) && !thr_lce.skip_succeeding_lce(stored_run[j], last_len))) {
                    last_len = compute_lce(stored_sample_pos[j], m-stored_it[j], last_len);
                }
                //else {
                //    verbose("using stored threshold LCE");
                //}
                write_len(last_len+1, lce_is_paused);
                write_len_segment(stored_it[j+1] - stored_it[j] - 1);
            }
        };

        for (size_t i = 1; i < m; ++i) {
            const auto c = pattern_at(m - i - 1);
            const size_t number_of_runs_of_c = this->bwt.number_of_letter(c);

            if (number_of_runs_of_c == 0) {
                write_len(0, lce_is_paused);
                write_ref(1); // why 1???
            } else if (pos < this->bwt.size() && this->bwt[pos] == c) {
                write_len(last_len+1, lce_is_paused);
                write_ref(last_ref-1);
            } else { // we jump
                const ri::ulint rank = this->bwt.rank(pos, c);

                size_t sa0, sa1; // SA indices of the previous and next entries with BWT[i] = c
                ri::ulint run1;

                if (rank > 0) {
                    sa0 = this->bwt.select(rank-1, c);
                }

                if (rank < number_of_runs_of_c) {
                    sa1 = this->bwt.select(rank, c);
                }

                run1 = this->bwt.run_of_position(sa1);

                const Triplet jump = [&] () -> Triplet {
                    if (rank == 0) {
                        // Check only succeeding -> we ignore thresholds in this case
                        return delay_succeeding_lce(sa1, run1, i);
                    } else if(rank >= number_of_runs_of_c) {
                        // Check only preceding -> we ignore thresholds in this case
                        return delay_preceding_lce(run1, rank, c, i);
                    }
                    // Check thresholds and boundary LCEs first
                    const size_t thr = thresholds[run1];
                    if (pos < thr) {
                        if (!lce_is_paused && thr_lce.skip_preceding_lce(run1, last_len)) {
                            const ri::ulint run0 = this->bwt.run_of_position(sa0);
                            const size_t ref0 = this->samples_last[run0];
                            const size_t len0 = last_len;
                            //const size_t ref1 = this->samples_start[run1];
                            //verbose("i = ", i, "last_ref = ", last_ref, "STORED UP for = ", ref0, " pos= ", pos, " thr= ", thr);
                            //verbose("i = ", i, "last_ref = ", last_ref, "STORED DOWN for = ", ref1);
                            return {sa0, ref0, len0};
                        } else {
                            return delay_preceding_lce(run1, rank, c, i);
                        }
                    } else {
                        if (!lce_is_paused && thr_lce.skip_succeeding_lce(run1, last_len)) {
                            const size_t ref1 = this->samples_start[run1];
                            const size_t len1 = last_len;
                            //verbose("i = ", i, "last_ref = ", last_ref, "STORED DOWN for = ", ref1);
                            return {sa1, ref1, len1};
                        } else {
                            return delay_succeeding_lce(sa1, run1, i);
                        }
                    }
                }();

                if (lce_cnt == lce_freq) {
                    const int last_delay = lce_freq - 1;
                    const size_t skipped_steps = (stored_it[last_delay] - stored_it[0]);
                    const int new_r_bound = last_len + skipped_steps;

                    const size_t lce = compute_lce(stored_sample_pos[last_delay], m-stored_it[last_delay], new_r_bound);

                    lce_is_paused = false;
                    if (lce == new_r_bound) {
                        //verbose("first delayed unbounded lce: ", compute_lce(stored_sample_pos[0], m-stored_it[0], 1000));
                        //verbose("last  delayed unbounded lce: ", compute_lce(stored_sample_pos[last_delay], m-stored_it[last_delay], 1000));
                        //verbose("lce: ", lce, " last_len+skipped_steps: ", last_len + skipped_steps, " last_len", last_len, " skipped_steps: ", skipped_steps, "i = ", i);
                        //verbose("stored_sample_pos:", stored_sample_pos[last_delay]);
                        // if we still process the same MEM, we know the lens right ahead
                        write_len_segment(skipped_steps+1);
                    } else  {
                        empty_stack();
                    }
                    lce_cnt = 0;
                } else {
                    write_len(1 + jump.len, lce_is_paused);
                }
                write_ref(jump.ref);
                pos = jump.sa;
            }
            pos = LF(pos, c); //! Perform one backward step
        }
        //verbose("lce_cnt: ", lce_cnt);
        //verbose("last len pushed: ", lens.back());

        if (lce_cnt > 0) {
            //verbose("emptying stack of size ", lce_cnt);
            // do the remaining LCEs
            const size_t last_delay = lce_cnt - 1;
            lce_is_paused = false;
            empty_stack();
            write_len_segment(m - 1 - stored_it[last_delay]);
        }
        //verbose("lens_size = ", lens.size(), " refs.size = ", refs.size());
        return std::make_pair(lens, refs);
    }

    /*
     * \param i position in the BWT
     * \param c character
     * \return lexicographic rank of cw in bwt
     */
    ulint LF(ri::ulint i, ri::uchar c)
    {
        // //if character does not appear in the text, return empty pair
        // if ((c == 255 and this->F[c] == this->bwt_size()) || this->F[c] >= this->F[c + 1])
        //     return {1, 0};
        //number of c before the interval
        ri::ulint c_before = this->bwt.rank(i, c);
        // number of c inside the interval rn
        ri::ulint l = this->F[c] + c_before;
        return l;
    }

    /* serialize the structure to the ostream
     * \param out     the ostream
     */
    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") // const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;

        out.write((char *)&this->terminator_position, sizeof(this->terminator_position));
        written_bytes += sizeof(this->terminator_position);
        written_bytes += my_serialize(this->F, out, child, "F");
        written_bytes += this->bwt.serialize(out);
        written_bytes += this->samples_last.serialize(out);

        written_bytes += samples_start.serialize(out, child, "samples_start");

        written_bytes += thresholds.serialize(out, child, "thresholds");

        written_bytes += thr_lce.serialize(out, child, "thr_lce");

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    /* load the structure from the istream
     * \param in the istream
     */
    void load(std::istream &in, const std::string& filename)
    {
        in.read((char *)&this->terminator_position, sizeof(this->terminator_position));
        my_load(this->F, in);
        this->bwt.load(in);
        this->r = this->bwt.number_of_runs();
        this->samples_last.load(in);
        this->samples_start.load(in);
        verbose("runof ", this->bwt.run_of_position(9734));

        thresholds.load(in, &this->bwt);
        thr_lce.load(in, &this->bwt);

        load_grammar(filename);
    }

    // // From r-index
    // ulint get_last_run_sample()
    // {
    //     return (samples_last[r - 1] + 1) % bwt.size();
    // }

protected :

    vector<ulint> build_F_(std::ifstream &heads, std::ifstream &lengths)
    {
        heads.clear();
        heads.seekg(0);
        lengths.clear();
        lengths.seekg(0);

        this->F = vector<ulint>(256, 0);
        int c;
        {
            ulint i = 0;
            while ((c = heads.get()) != EOF) {
                size_t length;
                lengths.read((char *)&length, 5);
                if (c > TERMINATOR)
                    this->F[c] += length;
                else {
                    this->F[TERMINATOR] += length;
                    this->terminator_position = i;
                }
                i++;
            }
        }
        for (ulint i = 255; i > 0; --i)
            this->F[i] = this->F[i - 1];
        this->F[0] = 0;
        for (ulint i = 1; i < 256; ++i)
            this->F[i] += this->F[i - 1];
        return this->F;
    }
};

#endif /* end of include guard: _MS_POINTERS_HH */
