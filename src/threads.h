#ifndef THREADS_H
#define THREADS_H

#include "Eigen/Core"
#include <iostream>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <vector>

// ThreadPool class declaration
class ThreadPool {
public:
    // Constructor and destructor
    ThreadPool(size_t num_threads = std::thread::hardware_concurrency());
    ~ThreadPool();

    // Function to enqueue a task
    void enqueueTask(std::function<void()> task);

private:
    // Worker threads
    std::vector<std::thread> threads_;

    // Task queue
    std::queue<std::function<void()>> tasks_;

    // Synchronization primitives
    std::mutex queue_mutex_;
    std::condition_variable cv_;

    // Stop flag
    bool stop_ = false;
};

#endif // THREADS_H
