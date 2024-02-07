#pragma once
#include <fcntl.h>
#include <unistd.h>

#include <cstring>
#include <iostream>

bool lockFile(int fd) {
  struct flock fl;
  fl.l_type = F_WRLCK;  // Write lock
  fl.l_whence = SEEK_SET;
  fl.l_start = 0;
  fl.l_len = 0;  // Lock the whole file

  if (fcntl(fd, F_SETLKW, &fl) == -1) {
    std::cerr << "Error locking file: " << strerror(errno) << std::endl;
    return false;
  }

  return true;
}

bool unlockFile(int fd) {
  struct flock fl;
  fl.l_type = F_UNLCK;
  fl.l_whence = SEEK_SET;
  fl.l_start = 0;
  fl.l_len = 0;  // Unlock the whole file

  if (fcntl(fd, F_SETLK, &fl) == -1) {
    std::cerr << "Error unlocking file: " << strerror(errno) << std::endl;
    return false;
  }

  return true;
}

// RAII wrapper for file locking
class FileLocker {
 public:
  FileLocker(const std::string& path) : path_(path) {
    fd_ = open(path_.c_str(), O_RDWR);
    if (fd_ == -1) {
      std::cerr << "Error opening file: " << path << " " << strerror(errno)
                << std::endl;
      throw std::runtime_error("Error opening file");
    }
    if (!lockFile(fd_)) {
      throw std::runtime_error("Error locking file");
    }
  }

  int getFd() const { return fd_; }

  ~FileLocker() {
    if (!unlockFile(fd_)) {
      std::cerr << "Error unlocking file" << std::endl;
    }
    close(fd_);
  }

 private:
  std::string path_;
  int fd_;
};

bool clearFile(int fileDescriptor) {
  if (ftruncate(fileDescriptor, 0) == -1) {
    std::cerr << "Failed to clear the file." << std::endl;
    return false;
  }
  return true;
}

#include <semaphore>
class SemaphoreLock {
 public:
  explicit SemaphoreLock(std::counting_semaphore<>& semaphore)
      : sem(semaphore) {
    // Acquire a semaphore slot
    sem.acquire();
  }

  ~SemaphoreLock() {
    // Release the semaphore slot
    sem.release();
  }

  // Delete copy constructor and copy assignment operator
  SemaphoreLock(const SemaphoreLock&) = delete;
  SemaphoreLock& operator=(const SemaphoreLock&) = delete;

  // Define move constructor and move assignment operator, if needed

 private:
  std::counting_semaphore<>& sem;
};