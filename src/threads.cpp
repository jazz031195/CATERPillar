#include "threads.h"
#include <utility>

// Constructor
ThreadPool::ThreadPool(int numThreads) : m_threads(numThreads), stop(false) {
    for (int i = 0; i < m_threads; i++) {
        threads.emplace_back([this] {
            std::function<void()> task;
            while (true) {
                
                std::unique_lock<std::mutex> lock(mtx);
                cv.wait(lock, [this] {
                    return !tasks.empty() || stop;
                });
                
                if (stop) {
                    return;
                }

                task = std::move(tasks.front());
                tasks.pop();
                lock.unlock();

                task();
            }
        });
    }
}

// Destructor
ThreadPool::~ThreadPool() {
    {
        std::unique_lock<std::mutex> lock(mtx);
        stop = true;
    }
    cv.notify_all();

    for (auto& th : threads) {
        th.join();
    }
}
