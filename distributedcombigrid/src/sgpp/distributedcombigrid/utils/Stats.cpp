#include "Stats.hpp"

namespace combigrid {

bool Stats::initialized_ = false;
bool Stats::finalized_ = false;
Stats::time_point Stats::init_time_;
Stats::time_point Stats::end_time_;
std::unordered_map<std::string, std::vector<Stats::Event>> Stats::event_;


void Stats::initialize() {
  assert(!initialized_);

  initialized_ = true;
  MPI_Barrier(MPI_COMM_WORLD);
  init_time_ = std::chrono::high_resolution_clock::now();
}

void Stats::finalize() {
  assert(initialized_);

  MPI_Barrier(MPI_COMM_WORLD);
  end_time_ = std::chrono::high_resolution_clock::now();

  finalized_ = true;
}

void Stats::startEvent(const std::string& name) {
  assert(initialized_);

  event_[name].emplace_back();
}

void Stats::stopEvent(const std::string& name) {
  assert(initialized_);
  // check if event was started
  assert(event_[name].back().end.time_since_epoch().count() == 0);

  event_[name].back().end = std::chrono::high_resolution_clock::now();
}

void Stats::write(const std::string& path) {
  assert(initialized_ && finalized_);

  using namespace std::chrono;

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::stringstream buffer;

  if (rank == 0) {
    buffer << "{" << std::endl;
  }

  buffer << "\"rank" << rank << "\": {" << std::endl;
  std::size_t event_count = 0;
  for (auto&& t = event_.begin(); t != event_.end(); ++t) {
    buffer << "\"" << t->first << "\"" << ": [" << std::endl;

    std::size_t element_count = 0;
    for (auto&& e = t->second.begin(); e != t->second.end(); ++e) {
      buffer << "[ "
          << duration_cast<milliseconds>(e->start - init_time_).count()
          << ", "
          << duration_cast<milliseconds>(e->end - init_time_).count()
          << " ]";
      if (++element_count != t->second.size()) {
        buffer << "," << std::endl;
      } else {
        buffer << std::endl;
      }
    }

    buffer << "]";
    if (++event_count != event_.size()) {
      buffer << "," << std::endl;
    } else {
      buffer << std::endl;
    }
  }
  buffer << "}";

  if (rank != size-1) {
    buffer << "," << std::endl;
  } else {
    buffer << std::endl << "}" << std::endl;
  }

  // get offset in file
  MPI_Offset len = buffer.str().size();
  MPI_Offset pos;
  MPI_Scan(&len, &pos, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  pos -= len;

  // get total file length
  MPI_Offset file_len;
  MPI_Allreduce(&len, &file_len, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

  // see: https://wickie.hlrs.de/platforms/index.php/MPI-IO
  MPI_Info info = MPI_INFO_NULL;
  if (file_len > 4*1024*1024 || 256 < size) {
      MPI_Info_create (&info);
      MPI_Info_set (info, "cb_align", "2");
      MPI_Info_set (info, "cb_nodes_list", "*:*");
      MPI_Info_set (info, "direct_io", "false");
      MPI_Info_set (info, "romio_ds_read", "disable");
      MPI_Info_set (info, "romio_ds_write", "disable");
      MPI_Info_set (info, "cb_nodes", "8");
  }

  // open file
  MPI_File fh;
  int err = MPI_File_open(MPI_COMM_WORLD, path.c_str(),
                          MPI_MODE_CREATE|MPI_MODE_WRONLY|MPI_MODE_EXCL,
                          info, &fh);
  if (err != MPI_SUCCESS)  {
    // file already existed, delete it and create new file
    if (rank == 0) {
      MPI_File_delete(path.c_str(), MPI_INFO_NULL);
    }
    MPI_File_open(MPI_COMM_WORLD, path.c_str(),
                  MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY,
                  info, &fh);
  }

  // write to single file with MPI-IO
  MPI_File_write_at_all(fh, pos, buffer.str().c_str(), (int)len, MPI_CHAR,
                        MPI_STATUS_IGNORE);
  MPI_File_close(&fh);
}

}
//end namespace combigrid
