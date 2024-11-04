#pragma once

#include <array>
#include <chrono>
#include <iostream>
#include <numeric>
#include <span>
#include <string>
#include <vector>

#define FMT_STRING_ALIAS 1
#include <fmt/format.h>
#if FMT_VERSION >= 60000
#include <fmt/chrono.h>
#else
#include <fmt/time.h>
#endif

// Disable compiler warnings we would like to retain for our own code.
#ifdef __GNUC__
#pragma GCC diagnostic push
#ifndef __clang__
#pragma GCC diagnostic ignored "-Wcast-function-type"
#pragma GCC diagnostic ignored "-Wuseless-cast"
#endif
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#endif

#include <mpi.h>

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

// create a new alias to the display !method with some template
// specialization of the prefix and test functions
#define ALIAS_DISPLAY(name, prefix, test)                              \
    template <typename... Args>                                        \
    inline void name(const char* format, const Args&... args) {        \
        display<prefix, test>(format, fmt::make_format_args(args...)); \
    }

/**
 * Smaller helper methods for printing and error checking in MPI programs
 */
namespace mpi {

/**
 * Return true if MPI has been initialized
 */
inline bool initialized() {
    static int flag = 0;
    if (flag == 0) {
        if (MPI_Initialized(&flag) != MPI_SUCCESS) {
            throw std::runtime_error("failed to read MPI status");
        }
    }
    return flag;
}

/**
 * Return the current MPI rank of the specified communicator
 */
inline int rank(MPI_Comm comm) {
    int r = -1;
    if (initialized()) {
        if (MPI_Comm_rank(comm, &r) != MPI_SUCCESS) {
            throw std::runtime_error("failed to read MPI rank");
        }
    }
    return r;
}

/**
 * Return the current MPI rank
 */
inline int rank() {
    static int r = -1;
    if (r == -1) {
        r = rank(MPI_COMM_WORLD);
    }
    return r;
}

/**
 * Return the size of the specified communicator
 */
inline int size(MPI_Comm comm) {
    int s = -1;
    if (initialized()) {
        if (MPI_Comm_size(comm, &s) != MPI_SUCCESS) {
            throw std::runtime_error("failed to read MPI size");
        }
    }
    return s;
}

/**
 * Return the MPI system size
 */
inline int size() {
    static int s = -1;
    if (s == -1) {
        s = size(MPI_COMM_WORLD);
    }
    return s;
}

/**
 * Return true if either MPI has not been initialized yet, or we are
 * rank 0.
 */
inline bool first() {
    return rank() <= 0;
}

/**
 * Return true if either MPI has not been initialized yet, or we are
 * rank 0.
 */
inline bool first(MPI_Comm comm) {
    return rank(comm) <= 0;
}

/**
 * Always return true
 */
constexpr bool any() {
    return true;
}

/**
 * Display some text using the fmt library
 *
 * \arg format the message format string
 * \arg args arguments to format the message with
 * \tparam prefix a function returning a string to prefix the message with
 * \tparam test a function that has to return true for the message to be printed
 */
template <std::string prefix(void), bool test(void)>
void display(const char* format, fmt::format_args args) {
    if (test()) {
        std::clog << prefix() << fmt::vformat(format, args) << std::endl;
    }
}

/**
 * Return some time information
 * \return a string containing the hour/minute/seconds of the time at invocation
 */
inline std::string info_prefix() {
    auto now = std::chrono::system_clock::now();
    auto now_c = std::chrono::system_clock::to_time_t(now);
    return fmt::format("{:%T}  ", *std::localtime(&now_c));
}

/**
 * Return an "empty" string
 * \return a string full of spaces
 */
inline std::string show_prefix() {
    return "          ";
}

/**
 * Return an string signifying a warning
 * \return a prefix string containing "Warning:"
 */
inline std::string warn_prefix() {
    return "Warning:  ";
}

/**
 * Return an string signifying an error
 * \return a prefix string containing "ERROR:"
 */
inline std::string fail_prefix() {
    return "ERROR:    ";
}

/**
 * Display some information with a timestamp
 *
 * \sa display
 * \arg format the message format string
 * \arg args arguments to format the message with
 */
ALIAS_DISPLAY(info, info_prefix, any)

/**
 * Display some information
 *
 * \sa display
 * \arg format the message format string
 * \arg args arguments to format the message with
 */
ALIAS_DISPLAY(show, show_prefix, any)

/**
 * Display a warning
 *
 * \sa display
 * \arg format the message format string
 * \arg args arguments to format the message with
 */
ALIAS_DISPLAY(warn, warn_prefix, any)

/**
 * Display an error
 *
 * \sa display
 * \arg format the message format string
 * \arg args arguments to format the message with
 */
ALIAS_DISPLAY(fail, fail_prefix, any)

/**
 * \brief Abort the current execution process
 *
 * If MPI is initialized, use MPI_Abort, otherwise throw a
 * std::runtime_error.
 *
 * \arg errorcode to use with MPI_Abort
 * \arg format the message format string
 * \arg args arguments to format the message with
 */
template <typename... Args>
void abort(int errorcode, const char* format, const Args&... args) {
    std::string msg = fmt::format(fmt::runtime(format), args...);
    fail(msg.c_str());
    if (initialized()) {
        int size = msg.size();
        MPI_Error_string(errorcode, const_cast<char*>(msg.c_str()), &size);
        MPI_Abort(MPI_COMM_WORLD, errorcode);
        exit(errorcode);
    } else {
        throw std::runtime_error(msg);
    }
}

/**
 * \brief Abort the current execution process
 *
 * If MPI is initialized, use MPI_Abort with error code -1, otherwise
 * throw a std::runtime_error.
 *
 * \arg format the message format string
 * \arg args arguments to format the message with
 */
template <typename... Args>
void abort(const char* format, const Args&... args) {
    abort(-1, format, args...);
}

/**
 * \brief Terminate the current execution process
 *
 * If MPI is initialized, use MPI_Finalize, otherwise throw a
 * std::runtime_error.
 *
 * \arg errorcode to exit with
 * \arg format the message format string
 * \arg args arguments to format the message with
 */
template <typename... Args>
void terminate(int errorcode, const char* format, const Args&... args) {
    std::string msg = fmt::format(fmt::runtime(format), args...);
    fail(msg.c_str());
    if (initialized()) {
        MPI_Finalize();
        exit(errorcode);
    } else {
        throw std::runtime_error(msg);
    }
}

/**
 * \brief Terminate the current execution process
 *
 * If MPI is initialized, use MPI_Finalize and exit with -1, otherwise
 * throw a std::runtime_error.
 *
 * \arg format the message format string
 * \arg args arguments to format the message with
 */
template <typename... Args>
void terminate(const char* format, const Args&... args) {
    terminate(-1, format, args...);
}

namespace {

template <
    typename T,
    std::enable_if_t<!(std::is_fundamental<T>::value || std::is_pointer<T>::value)>* = nullptr>
inline void dump(const T&) {
    warn("Argument supplied: <unprintable>");
}

template <typename T, std::enable_if_t<std::is_pointer<T>::value>* = nullptr>
inline void dump(const T& arg) {
    warn("Argument supplied: {}", fmt::ptr(arg));
}

template <typename T, std::enable_if_t<std::is_fundamental<T>::value>* = nullptr>
inline void dump(const T& arg) {
    warn("Argument supplied: {}", arg);
}

// template<typename T,
//          typename = std::enable_if_t<std::is_pointer<T>::value>>
// inline void dump_one(const T& arg) {
//     warn("Argument supplied: {}", fmt::pointer(arg));
// }

template <typename Arg, typename... Args>
inline void dump(const Arg& arg, const Args&... args) {
    dump(arg);
    dump(args...);
}

}  // namespace

/**
 * \brief Test a return value and abort if needed
 *
 * By default, test for MPI_SUCCESS.
 *
 * \sa abort
 * \tparam success the code to use for succes
 * \arg returncode the value to test against success
 * \arg format the message format string used when aborting
 * \arg args arguments to format the message with when aborting
 */
template <typename Function, typename... Args>
inline void ensure_(const char* file,
                    const char* name,
                    const int line,
                    Function& f,
                    const Args&... args) {
    int error_code = f(args...);
    if (error_code != MPI_SUCCESS) {
        char buf[BUFSIZ];
        int len, error_class;
        MPI_Error_class(error_code, &error_class);
        MPI_Error_string(error_class, buf, &len);
        warn("Error class:   {}", buf);
        MPI_Error_string(error_code, buf, &len);
        warn("Error message: {}", buf);
        warn("Function called:   {}:{}:{}", file, name, line);
        dump(args...);
        abort(error_code, "See previous warnings");
    }
}

// This is just sick, but oh well…
#define ensure(...) ensure_(__FILE__, __func__, __LINE__, __VA_ARGS__)

namespace {

template <typename T>
MPI_Datatype register_and_commit() {
    static_assert(std::is_trivially_copyable<T>::value,
                  "Only trivially copyable types are supported.");
    MPI_Datatype kind;
    ensure(MPI_Type_contiguous, sizeof(T), MPI_BYTE, &kind);
    ensure(MPI_Type_commit, &kind);
    return kind;
}

}  // namespace

template <typename T>
inline MPI_Datatype datatype() {
    static MPI_Datatype kind = register_and_commit<T>();
    return kind;
}

template <>
inline MPI_Datatype datatype<size_t>() {
    static_assert(sizeof(size_t) == sizeof(unsigned long int),
                  "Need to redefine MPI type for size_t");
    return MPI_UNSIGNED_LONG;
}

template <>
inline MPI_Datatype datatype<float>() {
    return MPI_FLOAT;
}

template <>
inline MPI_Datatype datatype<double>() {
    return MPI_DOUBLE;
}

/**
 * Methods that can be called everywhere, but are only executed on the
 * first rank
 */
namespace rank0 {
/**
 * Display some information with a timestamp if called from rank 0
 *
 * \sa display
 * \arg format the message format string
 * \arg args arguments to format the message with
 */
ALIAS_DISPLAY(info, info_prefix, first)

/**
 * Display some information if called from rank 0
 *
 * \sa display
 * \arg format the message format string
 * \arg args arguments to format the message with
 */
ALIAS_DISPLAY(show, show_prefix, first)

/**
 * Display a warning if called from rank 0
 *
 * \sa display
 * \arg format the message format string
 * \arg args arguments to format the message with
 */
ALIAS_DISPLAY(warn, warn_prefix, first)

/**
 * Display an error if called from rank 0
 *
 * \sa display
 * \arg format the message format string
 * \arg args arguments to format the message with
 */
ALIAS_DISPLAY(fail, fail_prefix, first)
}  // namespace rank0

/**
 * \brief Async counter that uses one-sided MPI communications
 */
template <size_t size, typename T>
class Counter {
  public:
    Counter(std::initializer_list<T> values = {}, MPI_Comm comm = MPI_COMM_WORLD) {
        MPI_Win_allocate(mpi::first(comm) ? size * sizeof(T) : 0,
                         sizeof(T),
                         MPI_INFO_NULL,
                         comm,
                         &buffer_,
                         &window_);
        if (mpi::first(comm)) {
            int i = 0;
            for (const auto& v: values) {
                buffer_[i] = v;
                ++i;
            }
        }
        MPI_Barrier(comm);
    }

    ~Counter() {
        MPI_Win_free(&window_);
    }

    T add(T value, int index) {
        T ret;
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, window_);
        MPI_Fetch_and_op(&value, &ret, datatype<T>(), 0, index, MPI_SUM, window_);
        MPI_Win_unlock(0, window_);
        return ret;
    }

  private:
    T* buffer_;
    MPI_Win window_;
};

/**
 * \brief Counter that uses whole-comm MPI communications
 */
template <size_t size, typename T>
class SyncCounter {
  public:
    SyncCounter(std::array<T, size> values = {}, MPI_Comm comm = MPI_COMM_WORLD)
        : comm_(comm)
        , values_(values) {}

    T add(T value, int index) {
        std::vector<T> vs(::mpi::size(comm_));
        ensure(MPI_Allgather, &value, 1, datatype<T>(), vs.data(), 1, datatype<T>(), comm_);

        auto local_offset = std::accumulate(vs.begin(), vs.begin() + rank(comm_), T{});
        auto total_count = std::accumulate(vs.begin(), vs.end(), T{});
        auto current_offset = values_[index] + local_offset;
        values_[index] += total_count;
        return current_offset;
    }

  private:
    MPI_Comm comm_;
    std::array<T, size> values_;
};

struct MPI {
    MPI(int* argc, char*** argv) {
        if (!mpi::initialized()) {
            int provided;
            MPI_Init_thread(argc, argv, MPI_THREAD_FUNNELED, &provided);
        }
    }

    ~MPI() {
        MPI_Finalize();
    }
};

/// RAII style wrapper for MPI handles.
/** Idiomatic C tends to expose handles to resources which need to be
 *  acquired (allocated, opened) and released (free, closed) manually. Idiomatic
 *  C++ on the other hand prefers the RAII approach to dealing with resources.
 *  RAII simply means that the life time of the resource is tied to the life
 *  time of the C++ object representing the resource.
 *
 *  The term *resource* shall refer to an object which needs to be
 *  acquired before use; must be released exactly once; and must never
 *  be used after it's been released.
 *
 *  The term *handle* shall refer to an object with shallow copy semantics that
 *  identifies the resource it's a handle for. Often either a pointer to an
 *  opaque object or an integer identifier.
 *
 *  Observations:
 *    - the details of acquiring the resource can often be
 *      configured by the user and is a very intentional act.
 *
 *    - there's only one way of releasing the resource, and it's easily
 *      forgotten.
 *
 *  Therefore, `Resource` doesn't provide means of acquiring the resource,
 *  only for taking over ownership of a resource from its handle.
 *
 *  This class offers a CRTP base class for wrapping MPI handles as
 *  resources. It implements the common API of resources:
 *
 *    - Take ownership of the resource through the constructor.
 *    - Forbid copying of the resource, but allow moving.
 *    - Provide an interface for dropping ownership.
 *    - Access to the raw handle, e.g., to be passed to other functions.
 *
 *  Note, `Resource` implements unique ownership. If shared ownership is needed
 *  one can maybe use a `std::shared_ptr<Derived>`.
 *
 *  Finally, sometimes one needs to release a resource before its scope ends. One
 *  can do so as follows:
 *
 *     SomeResource::free(resource.drop_ownership());
 *
 * \tparam Derived  The type of the C++ resource. It must support two static
 *                  methods:
 *                    - `invalid_handle()` which returns a handle representing
 *                      the null handle.
 *
 *                    - `free(Handle &)` which releases the resource.
 *
 * \tparam Handle   The type of the C handle to be wrapped as a C++ resource.
 */
template <class Derived, class Handle>
class Resource {
  public:
    /// Acquires ownership of the resource.
    explicit Resource(Handle handle)
        : handle_(handle) {}


    Resource(Resource&& other) noexcept
        : handle_(other.drop_ownership()) {}


    Resource(const Resource&) = delete;

    ~Resource() {
        if (handle_ != Derived::invalid_handle()) {
            Derived::free(handle_);
        }
    }

    Resource& operator=(Resource&& other) noexcept {
        handle_ = other.drop_ownership();
        return (*this);
    }

    Resource& operator=(const Resource&) = delete;

    /// Return the C handle to the resource.
    /** This does not drop ownership of the resource. Therefore, the handle
     *  must not be used to release the resource. Use `drop_ownership` to
     *  gain ownership of the resource.
     */
    Handle operator*() const {
        return handle_;
    }

    /// Returns the handle to and drops ownership of the resource.
    /** Note that after dropping ownership, this object represents
     *  the null resource.
     */
    Handle drop_ownership() {
        auto tmp = handle_;
        handle_ = Derived::invalid_handle();

        return tmp;
    }

  private:
    Handle handle_;
};


class File: public Resource<File, MPI_File> {
  private:
    using super = Resource<File, MPI_File>;

  public:
    // This is a bitfield, with names for mixed states.
    enum Mode : unsigned int { READ = 0x01, WRITE = 0x02, READ_WRITE = 0x03 };

    using super::Resource;
    using super::operator=;

    static File open(const std::string& filename, Mode mode, MPI_Comm comm);
    static int mpi_mode(Mode mode);

    static void free(MPI_File file);
    static MPI_File invalid_handle() noexcept;

    template <class T>
    void read_scalar_at(T& data, size_t offset) const {
        read_array_at(std::span{&data, 1ul}, offset);
    }

    template <class T>
    void write_scalar_at(const T& data, size_t offset) {
        write_array_at(std::span{&data, 1ul}, offset);
    }

    template <class T, std::size_t N>
    void read_array_at(const std::span<T, N>& data, size_t offset) const {
        static_assert(std::is_trivially_copyable<T>::value,
                      "Only trivially copyable types can be copied byte-by-byte.");

        auto count = data.size();

        if (count > 0) {
            auto type = mpi::datatype<T>();
            auto begin = data.data();
            mpi::ensure(MPI_File_read_at, **this, offset, begin, count, type, MPI_STATUS_IGNORE);
        }
    }

    template <class T, class A>
    void read_array_at(std::vector<T, A>& data, size_t offset) const {
        read_array_at(std::span<T>{data}, offset);
    }

    template <class T, std::size_t N>
    void write_array_at(const std::span<T, N>& data, size_t offset) {
        static_assert(std::is_trivially_copyable<T>::value,
                      "Only trivially copyable types can be copied byte-by-byte.");

        auto count = data.size();
        if (count > 0) {
            auto type = mpi::datatype<T>();
            auto begin = data.data();
            mpi::ensure(MPI_File_write_at, **this, offset, begin, count, type, MPI_STATUS_IGNORE);
        }
    }

    template <class T, class A>
    void write_array_at(const std::vector<T, A>& data, size_t offset) {
        write_array_at(std::span<const T>{data}, offset);
    }
};

}  // namespace mpi
