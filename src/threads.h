#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <thread>
#include <vector>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <future>

class ThreadPool {
private:
    int m_threads;                                  // Number of threads
    std::vector<std::thread> threads;              // Vector to hold threads
    std::queue<std::function<void()>> tasks;       // Task queue
    std::mutex mtx;                                // Mutex for synchronization
    std::condition_variable cv;                    // Condition variable
    bool stop;                                     // Stop flag

public:
    explicit ThreadPool(int numThreads);           // Constructor
    ~ThreadPool();                                 // Destructor

    // Template function to execute a task
    template<class F, class... Args>
    auto enqueueTask(F&& f, Args&&... args) -> std::future<decltype(f(args...))> {
        using return_type = decltype(f(args...));

        auto task = std::make_shared<std::packaged_task<return_type()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );

        std::future<return_type> res = task->get_future();

        {
            std::unique_lock<std::mutex> lock(mtx);
            tasks.emplace([task]() { (*task)(); });
        }

        cv.notify_one();
        return res;
    }
};

#endif
