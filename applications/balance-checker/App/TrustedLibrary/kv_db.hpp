#pragma once
#include <rocksdb/db.h>

#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include "rocksdb/write_batch.h"

template <typename K, typename V>
class KV_DB {
 public:
  KV_DB(const std::string& db_path) {
    rocksdb::Options options;
    options.create_if_missing = true;
    rocksdb::DB* rawDbPointer;
    rocksdb::Status status = rocksdb::DB::Open(options, db_path, &rawDbPointer);
    db_.reset(rawDbPointer);
    if (!status.ok()) {
      throw std::runtime_error("Unable to open/create database: " +
                               status.ToString());
    }

    // add a json meta data file under db_path to record the number of keys
  }

  void put(const K& key, const V& value) {
    rocksdb::Status status =
        db_->Put(rocksdb::WriteOptions(), serialize(key), serialize(value));
    if (!status.ok()) {
      throw std::runtime_error("Write failed: " + status.ToString());
    }
  }

  bool get(const K& key, V& value) {
    std::string valueStr;
    rocksdb::Status status =
        db_->Get(rocksdb::ReadOptions(), serialize(key), &valueStr);
    if (status.ok()) {
      value = deserialize<V>(valueStr);
      return true;
    } else if (status.IsNotFound()) {
      return false;
    } else {
      throw std::runtime_error("Read failed: " + status.ToString());
    }
  }

  bool get(const K& key, V& value, const V& defaultValue) {
    bool found = get(key, value);
    if (!found) {
      value = defaultValue;
    }
    return found;
  }

  bool del(const K& key) {
    rocksdb::Status status = db_->Delete(rocksdb::WriteOptions(), key);
    if (status.ok()) {
      return true;
    } else if (status.IsNotFound()) {
      return false;
    } else {
      throw std::runtime_error("Delete failed: " + status.ToString());
    }
  }

  rocksdb::WriteBatch newWriteBatch() { return rocksdb::WriteBatch(); }

  void writeBatch(rocksdb::WriteBatch& batch) {
    rocksdb::Status status = db_->Write(rocksdb::WriteOptions(), &batch);
    if (!status.ok()) {
      throw std::runtime_error("Write failed: " + status.ToString());
    }
  }

  class Iterator {
   public:
    Iterator(){};

    Iterator(rocksdb::Iterator* it) : it_(it) {}

    bool isValid() const { return it_.get() && it_->Valid(); }

    void seekToFirst() { it_->SeekToFirst(); }

    void next() { it_->Next(); }

    K key() const { return deserialize<K>(it_->key().ToString()); }

    V value() const { return deserialize<V>(it_->value().ToString()); }

    void reset() { it_.reset(); }

   private:
    std::unique_ptr<rocksdb::Iterator> it_;
  };

  Iterator getIterator() {
    return Iterator(db_->NewIterator(rocksdb::ReadOptions()));
  }

 private:
  std::unique_ptr<rocksdb::DB> db_;

  template <typename T>
  std::string serialize(const T& value) {
    if constexpr (std::is_same_v<T, std::string>) {
      return value;
    }
    std::ostringstream oss;
    oss << value;
    return oss.str();
  }

  template <typename T>
  static T deserialize(const std::string& str) {
    if constexpr (std::is_same_v<T, std::string>) {
      return str;
    }
    std::istringstream iss(str);
    T value;
    iss >> value;
    return value;
  }
};