#pragma once
#include "common/tracing/tracer.hpp"
#include <fstream>
#include <unordered_map>

extern bool g_disableProfiling;

enum class ProfilerEvent : uint32_t {
  Start,
  End,
  Count,
};

INLINE std::string toString(const ProfilerEvent& e) {
  switch(e) {
    case ProfilerEvent::Start: return "B";
    case ProfilerEvent::End: return "E";
    case ProfilerEvent::Count: return "C";
  }
  return "INVALID";
}


struct Empty { };
struct ProfilerLog {
  EventId eventId;
  uint64_t time;
  ProfilerEvent profType;
  uint64_t v;
};

// This profiles generates flamegraphs visible at https://github.com/jlfwong/speedscope
//
struct Profiler {
  using BlockProfiler_t = BlockTracer<Profiler, Empty>;
  using eint_t = std::underlying_type_t<EventId>;
  static Profiler g_Tracker;

  uint64_t blockRecCount;
  BlockProfiler_t* stack[MAX_TRACKER_REC];
  std::vector<ProfilerLog> logs;
  std::chrono::_V2::system_clock::time_point start;
  
  explicit Profiler() : stack{0}, logs{}, start{std::chrono::system_clock::now()} {
    std::cout << "Building profiler with " << TOTAL_TRACKERS << " total tracers." << std::endl;
    blockRecCount = 0;
  }

  INLINE void BeginBlock(BlockProfiler_t& that) {
    Assert(blockRecCount < MAX_TRACKER_REC);
    TrackEvent(that.eventId, ProfilerEvent::Start);
    stack[blockRecCount] = &that;
    blockRecCount += 1;
  }

  INLINE void EndBlock() {
    blockRecCount -= 1;
    stack[blockRecCount]->Finish();
  }

  INLINE void TrackEvent(const EventId eventId, ProfilerEvent eventType) {
    if (g_disableProfiling) return;
    auto currtime = std::chrono::system_clock::now();
    std::chrono::nanoseconds diff = currtime - start;
    logs.push_back({eventId, static_cast<uint64_t>(diff.count()), eventType});
  }

  INLINE void TrackValue(const EventId eventId, uint64_t val) {
    if (g_disableProfiling) return;
    auto currtime = std::chrono::system_clock::now();
    std::chrono::nanoseconds diff = currtime - start;
    logs.push_back({eventId, static_cast<uint64_t>(diff.count()), ProfilerEvent::Count, val});
  }
  
  INLINE void Measure(const EventId eventId, Empty) {
    if (g_disableProfiling) return;
    TrackEvent(eventId, ProfilerEvent::End);
  }

  ~Profiler() {
    Reset();
  }

  void Log(const std::string& filename_base= FLAMEGRAPHS_BASE_FOLDER "flameprofile") {
    if (g_disableProfiling) return;
    auto currtime = std::chrono::system_clock::now();
    std::chrono::nanoseconds endValue_ = currtime - start;
    uint64_t endValue = endValue_.count();

    // Log_SpeedScope(filename_base + "-speedscope.json", endValue);
    Log_Chrome(filename_base + "-chrome.json", endValue);
  }

  void Log_SpeedScope(const std::string& filename, uint64_t endValue) {
    std::cout << "Writing profile flamegraph to file \"" << filename << "\" (" << logs.size() << " events)" << std::endl;
    std::ofstream fout;
    std::vector<char> buf;
    int sz {1<<20};
    buf.resize(sz);
    memset(&buf[0],0,sz);
    fout.rdbuf()->pubsetbuf(&buf[0], sz);
    fout.open(filename);
    Assert(fout.is_open(), filename);
    fout << "\n\
      {\n\
        \"$schema\": \"https://www.speedscope.app/file-format-schema.json\",\n\
        \"version\": \"0.0.1\",\n\
        \"shared\": {\n\
          \"frames\": [\n";
    #define F(name) { \
        eint_t i = static_cast<eint_t>(EventId::name); \
        fout << "{\"name\": \"" << #name << "\"}, " << "\n"; \
      }
    #include "../tracer_events.hxx"
    fout << "{\"name\": \"invalid\"}";
    fout << "\n]\n},\n\
        \"profiles\": [\n\
          {\n\
            \"name\": \"profile.txt\",\n\
            \"type\": \"evented\",\n\
            \"unit\": \"none\",\n\
            \"startValue\": 0,\n\
            \"endValue\": ";
    fout << endValue << ", " << "\"events\": [\n";
    for (auto& log : logs) {
      bool isStart = log.profType == ProfilerEvent::Start;
      // This format does not support counters.
      if (log.profType != ProfilerEvent::Start && log.profType != ProfilerEvent::End) continue; 
      fout << "{\"at\": " << log.time << ","
           <<   "\"frame\": " << static_cast<eint_t>(log.eventId) << ","
           <<   "\"type\":" << (isStart ? "\"C\"" : "\"O\"")
           <<  "}," << "\n";
    }
    fout << "{\"at\": " << endValue << ","
           <<   "\"frame\": " << TOTAL_TRACKERS << ","
           <<   "\"type\":" << "\"O\""
           <<  "}," << std::endl;
    fout << "{\"at\": " << endValue << ","
           <<   "\"frame\": " << TOTAL_TRACKERS << ","
           <<   "\"type\":" << "\"C\""
           <<  "}" << std::endl;
    fout << "\n]\n}\n]\n}\n";
    fout << std::endl;
    fout.close();
  }

  void Log_Chrome(const std::string& filename, uint64_t endValue) {
    std::cout << "Writing profile flamegraph to file \"" << filename << "\" (" << logs.size() << " events)" << std::endl;
    std::ofstream fout;
    std::vector<char> buf{};
    int sz {1<<20};
    buf.resize(sz);
    fout.rdbuf()->pubsetbuf(&buf[0], sz);
    fout.open(filename);
    Assert(fout.is_open(), filename);
    std::unordered_map<EventId, std::string> m;
    fout << "[\n";

    #define F(name) { \
        m[EventId::name] = (#name); \
      }
    #include "../tracer_events.hxx"
    for (auto& log : logs) {
      if (log.profType == ProfilerEvent::Count) {
        fout << "{\"pid\": 1, "
             << "\"tid\": 1,";
      } else {
        fout << "{\"pid\": 0, "
             << "\"tid\": 0,";
      }
      fout << "\"ts\": " << log.time/1000.0  << ","
           <<   "\"name\": \"" << m[log.eventId] << "\","
           <<   "\"ph\":\"" << (toString(log.profType)) << "\"";
      if (log.profType == ProfilerEvent::Count) {
        fout << ", \"args\": {\"v\": " << log.v << "}";
      }
      fout <<  "}," << "\n";
    }
    fout << "{\"pid\": 0, \"tid\": 0, \"ts\": " << endValue/1000.0 << ","
           <<   "\"name\": \"EOL\","
           <<   "\"ph\":" << "\"i\""
           <<  "}" << std::endl;
    fout << "\n]\n";
    fout << std::endl;
    fout.close();
  }

  void Reset(bool log=true, const std::string& filename_base=FLAMEGRAPHS_BASE_FOLDER "flameProfile") {
    static int logCount = 0;
    if (log) {
      Log(filename_base + std::to_string(logCount));
      logCount++;
    }
    logs.clear();
    start = std::chrono::system_clock::now();
  }
};
