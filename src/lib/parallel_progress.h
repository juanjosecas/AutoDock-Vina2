/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

#ifndef VINA_PARALLEL_PROGRESS_H
#define VINA_PARALLEL_PROGRESS_H

#include <boost/thread/mutex.hpp>
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include "console.h"
#include "incrementable.h"

struct parallel_progress : public incrementable {
    parallel_progress(std::function<void(double)>* c = NULL)
        : callback(c), count(0), value(0), active(false), terminal(false),
          estimate(true), next_log_percent(10), last_width(0) {}

    void init(std::size_t n, bool estimate_remaining = true) {
        count = n;
        value = 0;
        active = count > 0;
        terminal = vina_console::interactive();
        estimate = estimate_remaining;
        started = last_render = clock::now();
        if (active) render(false);
    }

    void operator++() {
        if (!active) return;
        boost::mutex::scoped_lock lock(self);
        ++value;
        // Keep callback frequency and serialization unchanged.
        if (callback) (*callback)(fraction());
        const clock::time_point now = clock::now();
        const unsigned percent = static_cast<unsigned>(100 * fraction());
        if ((terminal && now - last_render >= std::chrono::milliseconds(200)) ||
            (!terminal && percent >= next_log_percent)) {
            render(false);
            next_log_percent = (percent / 10 + 1) * 10;
            last_render = now;
        }
    }

    void finish() {
        if (!active) return;
        boost::mutex::scoped_lock lock(self);
        render(true);
        active = false;
    }

private:
    typedef std::chrono::steady_clock clock;
    double fraction() const { return std::min(1.0, static_cast<double>(value) / count); }
    void render(bool complete) {
        const double elapsed = std::chrono::duration<double>(clock::now() - started).count();
        std::ostringstream line;
        line << (complete ? "Search complete " : "Search ");
        if (terminal) {
            const unsigned filled = static_cast<unsigned>(24 * fraction());
            line << '[' << std::string(filled, '=') << std::string(24 - filled, ' ') << "] ";
        }
        line << std::fixed << std::setprecision(1) << 100 * fraction()
             << "% of step budget | " << elapsed << " s";
        if (!complete && estimate && value > 0 && value < count)
            line << " | ETA ~" << elapsed * (1 - fraction()) / fraction() << " s";
        std::string text = line.str();
        if (terminal) {
            std::cout << '\r' << (complete ? vina_console::success(text) : text);
            if (text.size() < last_width) std::cout << std::string(last_width - text.size(), ' ');
            last_width = text.size();
            if (complete) std::cout << '\n';
        } else {
            std::cout << text << '\n';
        }
        std::cout.flush();
    }
    boost::mutex self;
    std::function<void(double)>* callback;
    std::size_t count, value;
    bool active, terminal, estimate;
    unsigned next_log_percent;
    std::size_t last_width;
    clock::time_point started, last_render;
};

#endif
