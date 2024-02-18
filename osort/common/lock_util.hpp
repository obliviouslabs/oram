#pragma once
#ifdef ENCLAVE_MODE
#ifndef TCS_NUM
#define TCS_NUM 1
#endif
#if TCS_NUM > 1
#include "sgx_spinlock.h"
// #include "sgx_thread.h"
struct Lock {
  sgx_spinlock_t _lock = SGX_SPINLOCK_INITIALIZER;
  void lock() { sgx_spin_lock(&_lock); }
  void unlock() { sgx_spin_unlock(&_lock); }
  // sgx_thread_mutex_t _lock = SGX_THREAD_MUTEX_INITIALIZER;
  // void lock() { sgx_thread_mutex_lock(&_lock); }
  // void unlock() { sgx_thread_mutex_unlock(&_lock); }
};
#else
struct Lock {
  void lock() {}
  void unlock() {}
};
#endif

#else
// #include <mutex>
// struct Lock {
//   std::mutex _lock;
//   void lock() { _lock.lock(); }
//   void unlock() { _lock.unlock(); }
//   Lock() {}
//   Lock(const Lock&) {}
// };
#include <atomic>
// TTASSpinlock
class Lock {
 public:
  Lock() : flag(false) {}
  Lock(const Lock&) {}

  void lock() {
    // First, test in a non-atomic way
    while (flag.load(std::memory_order_relaxed)) {
      // Spin without attempting to set the flag
      // Optionally, use std::this_thread::yield() to reduce contention
    }
    // Now attempt to set the flag atomically
    while (flag.exchange(true, std::memory_order_acquire)) {
      // If exchange returns true, the flag was already set, so keep spinning
    }
  }

  void unlock() { flag.store(false, std::memory_order_release); }

 private:
  std::atomic<bool> flag;
};
#endif
struct Critical {
  Lock& _lock;
  Critical(Lock& _lock) : _lock(_lock) { _lock.lock(); }
  ~Critical() { _lock.unlock(); }
};