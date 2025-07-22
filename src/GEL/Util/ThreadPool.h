//
// Created by Cem Akarsubasi on 4/25/25.
//

#ifndef GEL_UTIL_THREADPOOL_H
#define GEL_UTIL_THREADPOOL_H

#include <thread>
#include <condition_variable>
#include <semaphore>
#include <functional>
#include <queue>
#include <vector>

namespace Util
{

class IExecutor {
public:
    virtual ~IExecutor() = default;

    [[nodiscard]]
    virtual size_t size() const = 0;

    virtual void addTask(std::function<void()>&& task) = 0;

    virtual void waitAll() = 0;
};

class ImmediatePool final : public IExecutor {
    using thread_t = std::thread;
    size_t m_size;
    std::vector<thread_t> m_threads;
public:
    explicit ImmediatePool(const size_t size = std::thread::hardware_concurrency()) : m_size{size} {}
    ~ImmediatePool() override
    {
        waitAll();
    }

    [[nodiscard]] size_t size() const override
    {
        return m_size;
    }

    void addTask(std::function<void()>&& task) override
    {
        m_threads.emplace_back(task);
    }

    void waitAll() override
    {
        for (auto& thread : m_threads) {
            if (thread.joinable()) {
                thread.join();
            }
        }
    }
};

/// TODO: fix the remaining race condition that seem to only trigger on ARM
/// TODO: maybe just use std::thread since we don't even rely on jthread

/// Apple Clang does not support jthread without an additional argument
/// We don't rely meaningfully on jthread, so a C+11 thread fallback is included

/// @brief a non-generic threadpool implementation
class ThreadPool final : public IExecutor {
#if defined(__APPLE__)
    using thread_t = std::thread;
#else
    using thread_t = std::jthread;
#endif

    std::mutex m_queue_mutex;
    std::counting_semaphore<> m_queue_semaphore{0};

    std::atomic_int m_number_working = 0;
    std::binary_semaphore m_number_working_condition{0};

    std::queue<std::function<void()>> m_function_queue;
    std::vector<thread_t> m_threads;

#if defined(__APPLE__)
    std::atomic_bool m_should_stop{false};
#endif

public:
    ThreadPool() = delete;

    /// A Threadpool should never be copied
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    /// @brief Construct a threadpool with the given thread count
    ///
    /// @throws std::invalid_argument if thread_count is 0
    explicit ThreadPool(const uint32_t thread_count) : m_threads{
        std::vector<thread_t>(thread_count)}
    {
        if (thread_count == 0) {
            throw std::invalid_argument("thread_count must be greater than 0");
        }
        for (auto& thread : m_threads) {
#if !defined(__APPLE__)
            thread = thread_t([this](const std::stop_token& stop_token) {
                while (!stop_token.stop_requested()) {
#else
            thread = thread_t([this]() {
                while (!m_should_stop.load()) {
#endif
                    std::function<void()> task;
                    {
                        while (m_function_queue.empty() &&
#if !defined(__APPLE__)
                            !stop_token.stop_requested()
#else
                            !m_should_stop.load()
#endif
                            )
                            m_queue_semaphore.acquire();
                        std::lock_guard lock(m_queue_mutex);

                        if (m_function_queue.empty()) {
                            continue;
                        }
                        task = std::move(m_function_queue.front());
                        m_function_queue.pop();
                    }
                    // TODO: exception handling to prevent a thread from going down
                    std::invoke(task);

                    {
                        // Only one worker thread should wake up the main thread. Only one worker thread will observe
                        // this atomic variable as 1.
                        if (m_number_working.fetch_sub(1) == 1) {
                            m_number_working_condition.release();
                        }
                    }
                }
            });
        }
    }

    ~ThreadPool() override
    {
        this->cancelAll();
        for (auto& thread : m_threads) {
            if (thread.joinable()) {
                thread.join();
            }
        }
    }

    /// @brief number of threads
    /// @return number of threads
    [[nodiscard]] size_t size() const override
    {
        return m_threads.size();
    }

    /// @brief Adds a task to the queue
    /// @param task task to be executed
    void addTask(std::function<void()>&& task) override
    {
        m_number_working.fetch_add(1);
        {
            std::lock_guard lock(m_queue_mutex);
            m_function_queue.push(task);
        }
        m_queue_semaphore.release();
    }

    /// @brief Waits until all threads have finished
    void waitAll() override
    {
        while (m_number_working.load() != 0)
            m_number_working_condition.acquire();
    }


    /// @brief Cancels every thread
    ///
    /// After this function, no other tasks should be added
    void cancelAll()
    {
#if !defined(__APPLE__)
        for (auto& t : m_threads) {
            t.request_stop();
        }
#else
        m_should_stop.store(true);
#endif
        m_queue_semaphore.release(static_cast<long>(this->size()));
    }
};

/// A facade of ThreadPool for trivially running parallel algorithms in one thread
class DummyPool final : public IExecutor
{
public:
    DummyPool() = default;

    [[nodiscard]]
    size_t size() const override
    {
        return 1;
    }

    void addTask(std::function<void()>&& task) override
    {
        task();
    }

    void waitAll() override
    {

    }
};

} // namespace GEL::Util

#endif // GEL_UTIL_THREADPOOL_H
