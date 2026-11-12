#ifndef EXTERNAL_SORT_H
#define EXTERNAL_SORT_H

// Generic external (disk-backed) sort: chunked in-RAM sort + k-way merge via
// temporary files. Lets a caller sort/dedupe a sequence far larger than the
// memory budget it's given, at the cost of sequential disk I/O for the
// overflow. Record must be a POD/trivially-copyable type with no
// uninitialised padding bytes (records are
// written/read as raw bytes) so this stays format-agnostic -- it knows
// nothing about OSM, node ids, or ways; it is purely "sort records that
// don't fit in RAM".
//
// Usage:
//   ExternalSorter<Rec, Cmp> sorter(chunk_capacity, temp_dir, temp_prefix);
//   for (...) sorter.add(rec);
//   sorter.finish();
//   for (auto it = sorter.iterator(); it.has_next(); ) {
//     Rec r = it.next();
//     ...
//   }
// The iterator performs the k-way merge lazily (one record at a time), so
// the fully-sorted sequence is never itself materialized as one block. The
// iterator borrows the sorter's comparator, so the sorter must outlive it.

#include <vector>
#include <string>
#include <fstream>
#include <queue>
#include <cstdio>
#include <algorithm>
#include <stdexcept>

template <typename Record, typename Compare = std::less<Record> >
class ExternalSorter {
public:
  class Iterator {
  public:
    bool has_next() const { return !heap_.empty(); }

    Record next(){
      HeapEntry top = heap_.top();
      heap_.pop();
      Record result = top.rec;
      Record refill;
      if (read_one(top.run_idx, refill)) heap_.push(HeapEntry{refill, top.run_idx, cmp_});
      return result;
    }

  private:
    friend class ExternalSorter;

    struct HeapEntry {
      Record rec;
      size_t run_idx;
      const Compare *cmp;
      // priority_queue is a max-heap; negate the comparison to get the
      // smallest record out first (a min-heap over Record via Compare).
      bool operator<(const HeapEntry &other) const { return (*cmp)(other.rec, rec); }
    };

    explicit Iterator(const Compare *cmp) : cmp_(cmp) {}

    bool read_one(size_t run_idx, Record &out){
      std::ifstream &f = streams_[run_idx];
      if (!f.read(reinterpret_cast<char*>(&out), sizeof(Record))) return false;
      return true;
    }

    std::vector<std::ifstream> streams_;
    std::priority_queue<HeapEntry> heap_;
    const Compare *cmp_;
  };

  ExternalSorter(size_t chunk_capacity, std::string temp_dir, std::string temp_prefix, Compare cmp = Compare())
    : chunk_capacity_(chunk_capacity), temp_dir_(std::move(temp_dir)), temp_prefix_(std::move(temp_prefix)), cmp_(cmp) {
    buffer_.reserve(chunk_capacity_);
  }

  ~ExternalSorter(){
    for (const auto &path : run_paths_) std::remove(path.c_str());
  }

  ExternalSorter(const ExternalSorter&) = delete;
  ExternalSorter& operator=(const ExternalSorter&) = delete;

  void add(const Record &r){
    buffer_.push_back(r);
    if (buffer_.size() >= chunk_capacity_) flush_run();
  }

  // Call once after all add() calls, before iterating. Releases the
  // buffer's reserved capacity (not just its size) since no more add()
  // calls are expected afterward -- otherwise a finished-but-still-in-scope
  // sorter keeps holding a full chunk's worth of reserved (unused) memory
  // for as long as its iterator is being consumed, which adds up if another
  // sorter is actively buffering at the same time (e.g. a merge-join
  // feeding a second external sort).
  void finish(){
    if (!buffer_.empty()) flush_run();
    std::vector<Record>().swap(buffer_);
  }

  size_t n_runs() const { return run_paths_.size(); }

  Iterator iterator(){
    Iterator it(&cmp_);
    it.streams_.resize(run_paths_.size());
    for (size_t i = 0; i < run_paths_.size(); i++){
      it.streams_[i].open(run_paths_[i].c_str(), std::ios::binary);
      if (!it.streams_[i])
        throw std::runtime_error("Cannot reopen external sort run " + std::to_string(i) + " of " +
                                  std::to_string(run_paths_.size()) + ": " + run_paths_[i] +
                                  " (on Windows the per-process limit on concurrently open files is ~512, "
                                  "so use a larger chunk capacity / fewer runs)");
      Record first;
      if (it.streams_[i].read(reinterpret_cast<char*>(&first), sizeof(Record)))
        it.heap_.push(typename Iterator::HeapEntry{first, i, &cmp_});
    }
    return it;
  }

private:
  void flush_run(){
    std::sort(buffer_.begin(), buffer_.end(), cmp_);
    std::string path = temp_dir_ + "/" + temp_prefix_ + "_run" + std::to_string(run_paths_.size()) + ".tmp";
    std::ofstream f(path.c_str(), std::ios::binary);
    if (!f) throw std::runtime_error("Cannot create external sort temp file: " + path);
    f.write(reinterpret_cast<const char*>(buffer_.data()), buffer_.size() * sizeof(Record));
    run_paths_.push_back(path);
    buffer_.clear();
  }

  size_t chunk_capacity_;
  std::string temp_dir_;
  std::string temp_prefix_;
  Compare cmp_;
  std::vector<Record> buffer_;
  std::vector<std::string> run_paths_;
};

#endif
