#pragma once
#ifdef ENCLAVE_MODE
#ifndef TCS_NUM
#define TCS_NUM 1
#endif
#if TCS_NUM > 1
#include "sgx_spinlock.h"
struct Lock {
  sgx_spinlock_t _lock = SGX_SPINLOCK_INITIALIZER;
  void lock() { sgx_spin_lock(&_lock); }
  void unlock() { sgx_spin_unlock(&_lock); }
};
#else
struct Lock {
  void lock() {}
  void unlock() {}
};
#endif

#else
#include <mutex>
typedef std::mutex Lock;
#endif
struct Critical {
  Lock& _lock;
  Critical(Lock& _lock) : _lock(_lock) { _lock.lock(); }
  ~Critical() { _lock.unlock(); }
};