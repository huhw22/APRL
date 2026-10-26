#define _GNU_SOURCE
#include <dlfcn.h>
#include <errno.h>
#include <fcntl.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/syscall.h>
#include <time.h>
#include <unistd.h>

#include <hdf5.h>

/* Test-only interposition. Nothing activates without the exact opt-in flag,
 * a recognized MPI rank environment, and an existing marker directory. No
 * MPI calls are made here. Production code and its watchdog are unchanged. */
static int active = 0;
static int rank_id = -1;
static int target_rank = -1;
static int write_phase = 0;
static int checkpoint_opened = 0;
static int injected = 0;
static int exit_fd = -1;
static const char *marker_directory = NULL;

static uint64_t monotonic_ns(void)
{
  struct timespec now;
  if (clock_gettime(CLOCK_MONOTONIC, &now) != 0) return 0;
  return (uint64_t)now.tv_sec * UINT64_C(1000000000) + (uint64_t)now.tv_nsec;
}

static int parse_rank(const char *value)
{
  char *end = NULL;
  long number;
  if (!value || !*value) return -1;
  errno = 0;
  number = strtol(value, &end, 10);
  if (errno || !end || *end || number < 0 || number > 1000000) return -1;
  return (int)number;
}

__attribute__((constructor)) static void configure_fault(void)
{
  const char *enabled = getenv("APRL_CHECKPOINT_FAULT_ENABLE");
  const char *phase = getenv("APRL_CHECKPOINT_FAULT_PHASE");
  const char *binary = getenv("APRL_CHECKPOINT_FAULT_BINARY");
  const char *rank_names[] = {"OMPI_COMM_WORLD_RANK", "PMI_RANK",
                             "PMIX_RANK", "SLURM_PROCID"};
  char path[4096];
  char executable[4096];
  ssize_t executable_length;
  unsigned int index;
  if (!enabled || strcmp(enabled, "1") != 0) return;
  /* The preload also reaches mpirun/srun and their helper processes. Only
   * the exact simulation executable may reserve per-rank evidence files. */
  if (!binary) return;
  executable_length = readlink("/proc/self/exe", executable, sizeof(executable)-1);
  if (executable_length < 0 || executable_length >= (ssize_t)sizeof(executable)-1)
    return;
  executable[executable_length] = '\0';
  if (strcmp(executable, binary) != 0) return;
  for (index = 0; index < sizeof(rank_names)/sizeof(rank_names[0]); ++index)
    {
      rank_id = parse_rank(getenv(rank_names[index]));
      if (rank_id >= 0) break;
    }
  target_rank = parse_rank(getenv("APRL_CHECKPOINT_FAULT_RANK"));
  marker_directory = getenv("APRL_CHECKPOINT_FAULT_MARKERS");
  if (rank_id < 0 || target_rank < 0 || !marker_directory ||
      !phase || (strcmp(phase, "create") && strcmp(phase, "write"))) return;
  write_phase = strcmp(phase, "write") == 0;
  if (snprintf(path, sizeof(path), "%s/rank-%06d.exit", marker_directory,
               rank_id) >= (int)sizeof(path)) return;
  exit_fd = open(path, O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC, 0600);
  if (exit_fd < 0) return;
  active = 1;
}

static int checkpoint_filename(const char *name)
{
  const char *base;
  size_t length;
  if (!name) return 0;
  base = strrchr(name, '/');
  base = base ? base + 1 : name;
  length = strlen(base);
  return strncmp(base, "step-", 5) == 0 && length >= 11 &&
    strcmp(base + length - 11, ".h5.partial") == 0;
}

static void mark_entry(void)
{
  char path[4096];
  uint64_t when = monotonic_ns();
  int fd;
  if (snprintf(path, sizeof(path), "%s/rank-%06d.enter", marker_directory,
               rank_id) >= (int)sizeof(path)) return;
  fd = open(path, O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC, 0600);
  if (fd >= 0) { (void)write(fd, &when, sizeof(when)); close(fd); }
}

static void block_selected_rank(void)
{
  char path[4096];
  char record[256];
  int fd, length;
  struct timespec remaining = {20, 0};
  if (!active || rank_id != target_rank || injected) return;
  injected = 1;
  if (snprintf(path, sizeof(path), "%s/injected.json", marker_directory) <
      (int)sizeof(path))
    {
      fd = open(path, O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC, 0600);
      if (fd >= 0)
        {
          length = snprintf(record, sizeof(record),
            "{\"rank\":%d,\"pid\":%ld,\"phase\":\"%s\",\"block_seconds\":20,"
            "\"monotonic_ns\":%llu}\n", rank_id, (long)getpid(),
            write_phase ? "write" : "create",
            (unsigned long long)monotonic_ns());
          if (length > 0 && length < (int)sizeof(record))
            (void)write(fd, record, (size_t)length);
          close(fd);
        }
    }
  while (nanosleep(&remaining, &remaining) != 0 && errno == EINTR) {}
}

hid_t H5Fcreate(const char *name, unsigned int flags, hid_t creation, hid_t access)
{
  typedef hid_t (*Create)(const char *, unsigned int, hid_t, hid_t);
  static Create real_create = NULL;
  hid_t file;
  int checkpoint = active && checkpoint_filename(name);
  if (!real_create) *(void **)(&real_create) = dlsym(RTLD_NEXT, "H5Fcreate");
  if (!real_create) return (hid_t)-1;
  if (checkpoint)
    {
      mark_entry();
      if (!write_phase) block_selected_rank();
    }
  file = real_create(name, flags, creation, access);
  if (checkpoint && file >= 0) checkpoint_opened = 1;
  return file;
}

herr_t H5Dwrite(hid_t dataset, hid_t memory_type, hid_t memory_space,
                hid_t file_space, hid_t transfer, const void *buffer)
{
  typedef herr_t (*Write)(hid_t, hid_t, hid_t, hid_t, hid_t, const void *);
  static Write real_write = NULL;
  if (!real_write) *(void **)(&real_write) = dlsym(RTLD_NEXT, "H5Dwrite");
  if (!real_write) return (herr_t)-1;
  if (active && write_phase && checkpoint_opened) block_selected_rank();
  return real_write(dataset, memory_type, memory_space, file_space, transfer, buffer);
}

/* The actual watchdog calls _exit(124). Record that fact with async-safe
 * clock_gettime/write, then preserve immediate Linux process-wide exit. No
 * artificial delay is inserted to make ranks appear to pass together. */
void _exit(int status)
{
  if (active && exit_fd >= 0)
    {
      struct ExitRecord { uint32_t magic; int32_t status; uint64_t time_ns; };
      const struct ExitRecord record = {UINT32_C(0x4150524c), status, monotonic_ns()};
      (void)write(exit_fd, &record, sizeof(record));
    }
  (void)syscall(SYS_exit_group, status);
  __builtin_unreachable();
}
