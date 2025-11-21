//
// Created on 4/24/25.
//
#pragma once
#ifndef GEL_HMESH_RSRTIMER_H
#define GEL_HMESH_RSRTIMER_H

#include <chrono>
#include <iostream>
#include <string>
#include <unordered_map>

namespace Util {
    /// Timer class
    /// FIXME: The other ancient timer should probably be removed from the library and replaced with this
    struct RSRTimer {
        RSRTimer() = default;

        //RSRTimer(const RSRTimer&) = delete;
        //RSRTimer& operator=(const RSRTimer&) = delete;

    private:
        using time_point = std::chrono::high_resolution_clock::time_point;

        struct Time {
            time_point start;
            time_point end;
        };

        std::unordered_map<std::string, Time> m_times;

    public:
        RSRTimer& create(const std::string& name) {
            m_times.insert_or_assign(
                name,
                Time{std::chrono::high_resolution_clock::now(), std::chrono::high_resolution_clock::now()}
                );

            return *this;
        }

        RSRTimer& start(const std::string& name) {
            if (!m_times.contains(name)) {
                this->create(name);
            }
            m_times[name].start = std::chrono::high_resolution_clock::now();
            return *this;
        }

        void end(const std::string& name) {
            //int idx = idx_map[name];
            if (m_times.contains(name)) {
                m_times[name].end = std::chrono::high_resolution_clock::now();
            }
        }

        void show() const
        {
            for (const auto& [name, time] : m_times) {
                const auto count = std::chrono::duration_cast<std::chrono::milliseconds>(time.end - time.start).count();
                std::cout << "Spent "<< count << " milliseconds on " << name << "\n";
            }
        }
    };

    struct TimerGuard {
        TimerGuard() = delete;

        TimerGuard(RSRTimer& timer, const std::string& name) : timer(timer), name(name)
        {
            timer.start(name);
        }

        TimerGuard(const TimerGuard& other) = delete;

        ~TimerGuard()
        {
            timer.end(name);
        }
    private:
        RSRTimer& timer;
        std::string name;
    };
}

#endif //GEL_HMESH_RSRTIMER_H
