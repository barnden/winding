#pragma once

#include <print>
#include <queue>
#include <thread>
#include <vector>
#include <mutex>
#include <condition_variable>

template <typename Job>
class ThreadPool {
    std::vector<std::thread> m_threads;

    std::queue<Job> m_queue;
    std::mutex m_queue_mutex;
    std::condition_variable m_condition;

    bool m_terminated = false;

    void worker()
    {
        while (true) {
            Job job;

            {
                auto lock = std::unique_lock(m_queue_mutex);

                m_condition.wait(lock,
                                 [&] { return !m_queue.empty() || m_terminated; });

                if (m_terminated)
                    return;

                job = m_queue.front();
                m_queue.pop();
            }

            job();
        }
    }

public:
    ThreadPool()
    {
        auto const num_threads = std::thread::hardware_concurrency();
        m_threads.reserve(num_threads);

        for (auto i = 0uz; i < num_threads; i++) {
            m_threads.push_back(std::thread(&ThreadPool::worker, this));
        }
    }

    ~ThreadPool()
    {
        stop();
    }

    void queue(Job&& job)
    {
        {
            auto lock = std::unique_lock(m_queue_mutex);
            m_queue.push(job);
        }

        m_condition.notify_one();
    }

    auto busy() -> bool
    {
        bool is_busy;
        {
            auto lock = std::unique_lock(m_queue_mutex);
            is_busy = !m_queue.empty();
        }

        return is_busy;
    }

    void stop()
    {
        {
            auto lock = std::unique_lock(m_queue_mutex);
            m_terminated = true;
        }

        m_condition.notify_all();

        for (auto&& thread : m_threads)
            thread.join();

        m_threads.clear();
    }
};
