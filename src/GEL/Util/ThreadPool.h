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

// Might change in the future
using ThreadPool = ImmediatePool;

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
