//
// This file contains helper routines and "small" classes which don't have their own source file.
//

#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <dirent.h>
#include <cstring>

#include "ch_vdif_assembler_internals.hpp"

using namespace std;

namespace ch_vdif_assembler {
#if 0
};  // pacify emacs c-mode!
#endif


// -------------------------------------------------------------------------------------------------
//
// Thread utils


thread_base::thread_base(const string &name_)
    : name(name_)
{ }


// static member function
void thread_base::_spawn(thread_base *thread)
{
    xassert(thread);

#if THREAD_DEBUG >= 1
    cout << (string("spawning ") + thread->name) << endl;
#endif

    int err = pthread_create(&thread->pthread, NULL, thread_base::_pthread_main, reinterpret_cast<void *> (thread));

    if (err != 0) {
	cout << thread->name << ": pthread_create() failed!" << endl;
	delete thread;
	throw runtime_error("pthread_create() failed!\n");
    }

}

//
// static member function
//
// FIXME think carefuly about corner cases, e.g. exception thrown in stringstream::str()
//
void *thread_base::_pthread_main(void *arg)
{
    const char *status = "exited normally";
    thread_base *p = reinterpret_cast<thread_base *> (arg);
    xassert(p);

#if THREAD_DEBUG >= 1
    cout << (p->name + string(": thread starting")) << endl;
#endif
    
    p->timer.start_running();

    try {
	p->thread_body();
    } catch (...) {
	status = "threw exception";
    }

    p->timer.stop_running();

#if THREAD_DEBUG >= 1
    // FIXME report busyfrac by default in assembler/processing threads
    stringstream ss;
    ss << p->name << ": " << status
       << ", busyfrac=" << p->timer.busyfrac() 
       << "\n";

    cout << ss.str() << flush;
#endif

    delete p;
    return NULL;
}


// -------------------------------------------------------------------------------------------------
//
// Chunk pools


chunk_pool::chunk_pool(int nbytes_per_chunk_, bool set_zero_)
    : nbytes_per_chunk(nbytes_per_chunk_), set_zero(set_zero_)
{
    xassert(nbytes_per_chunk > 0);
    pthread_mutex_init(&mutex, NULL);
}


chunk_pool::~chunk_pool()
{
    clear();
    pthread_mutex_destroy(&mutex);
}


uint8_t *chunk_pool::get_chunk()
{
    uint8_t *ret = NULL;
    pthread_mutex_lock(&mutex);
    
    if (pointer_pool.size() > 0) {
	ret = pointer_pool[pointer_pool.size()-1];
	pointer_pool.pop_back();
	pthread_mutex_unlock(&mutex);
    }
    else {
	pthread_mutex_unlock(&mutex);
    
	if (posix_memalign(reinterpret_cast<void **> (&ret), constants::cache_line_size, nbytes_per_chunk) != 0)
	    throw std::runtime_error("couldn't allocate memory");
    }

    xassert(ret != NULL);

    if (set_zero)
	memset(ret, 0, nbytes_per_chunk);

    return ret;
}


void chunk_pool::put_chunk(uint8_t *buf)
{
    xassert(buf != NULL);

    pthread_mutex_lock(&mutex);
    pointer_pool.push_back(buf);
    pthread_mutex_unlock(&mutex);
}


void chunk_pool::clear()
{
    vector<uint8_t *> v;

    pthread_mutex_lock(&mutex);
    pointer_pool.swap(v);
    pthread_mutex_unlock(&mutex);

    for (unsigned int i = 0; i < v.size(); i++)
	free(v[i]);
}


vdif_chunk_pool::vdif_chunk_pool(int packet_count_, bool set_zero_)
    : chunk_pool(packet_count_ * constants::packet_nbytes + 2 * vdif_chunk::pad, set_zero_),
      packet_count(packet_count_)
{
    xassert(packet_count > 0);
}


assembled_chunk_pool::assembled_chunk_pool(int assembler_nt_)
    : chunk_pool(assembler_nt_ * 2 * constants::chime_nfreq, true),
      assembler_nt(assembler_nt_)
{
    xassert(assembler_nt > 0);
    xassert(assembler_nt % constants::cache_line_size == 0);
}


// -------------------------------------------------------------------------------------------------
//
// vdif_chunk


vdif_chunk::vdif_chunk(const shared_ptr<vdif_chunk_pool> &pool_, int seq_id_)
    : pool(pool_), seq_id(seq_id_)
{
    this->buf0 = pool->get_chunk();
    this->buf = this->buf0 + pad;
    this->capacity = pool->packet_count;
    this->size = 0;
    this->is_on_disk = false;

    // placeholder; will eventually get set in assembler_nerve_center::stream_put_chunk() or assembler_nerve_center::trigger()
    this->want_on_disk = false;
}


vdif_chunk::vdif_chunk(const string &filename, int seq_id)
{
    struct stat s;
    int err = stat(filename.c_str(), &s);

    if (err < 0) {
	cout << filename << ": stat() failed: " << strerror(errno) << endl;
	throw runtime_error(strerror(errno));
    }
    
    if (s.st_size % constants::packet_nbytes) {
	cout << filename << ": warning: file size (=" << s.st_size << ")"
	     << " is not divisible by packet_nbytes (=" << constants::packet_nbytes << "),"
	     << " truncating file" << endl;
	// we now fall through here, instead of throwing an exception...
    }

    int npackets = s.st_size / constants::packet_nbytes;
    int nbytes_unpadded = npackets * constants::packet_nbytes;

    int fd = open(filename.c_str(), O_RDONLY);

    if (fd < 0) {
	cout << filename << ": open() failed: " << strerror(errno) << endl;
	throw runtime_error(strerror(errno));
    }

    this->capacity = npackets;
    this->buf0 = reinterpret_cast<uint8_t *> (malloc(nbytes_unpadded + 2*pad));
    this->buf = buf0 + pad;
    this->size = capacity;   // in anticipation of reading the whole buffer
    this->seq_id = seq_id;
    this->is_on_disk = true;
    this->want_on_disk = false;

    int pos = 0;

    while (pos < nbytes_unpadded) {
	int ret = read(fd, &buf[pos], nbytes_unpadded - pos);

	if (ret <= 0) {
	    close(fd);
	    free(this->buf0);
	    this->buf0 = this->buf = NULL;
	    this->capacity = this->size = 0;

	    cout << filename << ": error in read(), or unexpected end-of-file: " << strerror(errno) << endl;
	    throw runtime_error(strerror(errno));
	}
	
	pos += ret;
    }

    close(fd);
}


void vdif_chunk::write(const string &filename)
{
    int fd = open(filename.c_str(), O_WRONLY | O_CREAT, 0666);

    if (fd < 0) {
	cout << filename << ": " << strerror(errno) << endl;
	throw runtime_error(strerror(errno));
    }

    ssize_t pos = 0;
    ssize_t sz = capacity * constants::packet_nbytes;

    while (pos < sz) {
	int ret = ::write(fd, &buf[pos], sz);

	if (ret <= 0) {
	    close(fd);
	    cout << filename << ": error on write(): " << strerror(errno) << endl;
	    throw runtime_error(strerror(errno));
	}

	pos += ret;
    }

    close(fd);
}


vdif_chunk::~vdif_chunk()
{
    if (pool)
	pool->put_chunk(buf0);
    else
	free(buf0);
	
    this->buf0 = this->buf = NULL;
    this->capacity = this->size = 0;
}


// -------------------------------------------------------------------------------------------------
//
// assembled_chunk


assembled_chunk::assembled_chunk(const shared_ptr<assembled_chunk_pool> &pool_, int64_t t0_)
    : pool(pool_), 
      buf(pool_->get_chunk()),
      t0(t0_), nt(pool_->assembler_nt), pcount(0)
{ }


assembled_chunk::~assembled_chunk()
{
    pool->put_chunk(const_cast<uint8_t *> (buf));
    const_cast<const uint8_t *&> (this->buf) = NULL;
}


bool assembled_chunk::is_zero() const
{
    // FIXME could make an assembly language kernel for this, but it's 
    // currently only used for non speed critical unit testing

    const int64_t *p64 = reinterpret_cast<const int64_t *> (buf);
    int n64 = (constants::chime_nfreq * 2 * nt) / 8;
    bool ret = true;

    // optimize for case where 'true' is returned
    for (int i = 0; i < n64; i++)
	ret = p64[i] ? false : ret;  // compiler emits conditional move, not branch

    return ret;
}


bool assembled_chunk::is_equal(const assembled_chunk &a) const
{
    xassert(this->t0 == a.t0);
    xassert(this->nt == a.nt);

    int nbytes = constants::chime_nfreq * 2 * nt;
    return memcmp(this->buf, a.buf, nbytes) == 0;
}


// -------------------------------------------------------------------------------------------------
//
// vdif_assembler


vdif_assembler::vdif_assembler(bool write_to_disk, int rbuf_size, int abuf_size, int assembler_nt)
    : nc(make_shared<assembler_nerve_center> (write_to_disk, rbuf_size, abuf_size, assembler_nt)),
      killer(make_shared<assembler_killer> (nc, "last reference to assembler dropped"))
{ 
    string dataset_name = make_dataset_name();

    spawn_assembler_thread(nc);
    
    for (int ithread = 0; ithread < constants::num_disks; ithread++) {
	string data_dir = make_data_dir(dataset_name, ithread);
	spawn_disk_writer_thread(nc, data_dir, ithread);
    }
}


vdif_assembler::~vdif_assembler()
{
    int ndrops_assembler, ndrops_disk_writer, ntot;
    nc->get_drop_stats(ndrops_assembler, ndrops_disk_writer, ntot);

    double afrac = (ntot > 0) ? ((double)ndrops_assembler / (double)ntot) : 0.0;
    double dfrac = (ntot > 0) ? ((double)ndrops_disk_writer / (double)ntot) : 0.0;

    stringstream ss;
    ss << "assembler: " << ntot << " buffers processed, "
       << ndrops_assembler << " assembler drops (frac=" << afrac << "), "
       << ndrops_disk_writer << " disk writer drops (frac=" << dfrac << ")\n";

    string s = ss.str();
    cout << s.c_str() << flush;
}


void vdif_assembler::register_processor(const shared_ptr<vdif_processor> &p)
{
    spawn_processing_thread(nc, p);
}


void vdif_assembler::run(const shared_ptr<vdif_stream> &s)
{
    this->start_async(s);
    this->wait_until_end();
}


void vdif_assembler::start_async(const shared_ptr<vdif_stream> &s)
{
    cout << "assembler: start\n" << flush;
    nc->stream_start(s->is_realtime);

    try {
	s->spawn_threads(this->nc);
    }
    catch (...) {
	nc->kill_assembler("error when spawning stream thread(s)");
	throw runtime_error("error when spawning stream thread(s)");
    }
}


void vdif_assembler::wait_until_end()
{
    nc->wait_until_end();
    cout << "assembler: normal exit\n" << flush;
}


// -------------------------------------------------------------------------------------------------
//
// assembler_killer helper class


assembler_killer::assembler_killer()
    : nc(), killmsg("[ should never see this ]")
{ }


assembler_killer::assembler_killer(const shared_ptr<assembler_nerve_center> &nc_, const char *killmsg_)
    : nc(nc_), killmsg(killmsg_)
{ }


assembler_killer::~assembler_killer()
{
    if (nc)
	nc->kill_assembler(killmsg);
}

void assembler_killer::set_victim(const shared_ptr<assembler_nerve_center> &nc_, const char *killmsg_)
{
    nc = nc_;
    killmsg = killmsg_;
}

void assembler_killer::let_live()
{
    nc = shared_ptr<assembler_nerve_center> ();
}


// -------------------------------------------------------------------------------------------------
//
// processor_handle


processor_handle::processor_handle(const string &name_, const shared_ptr<assembler_nerve_center> &nc_)
    : name(name_), nc(nc_),
      ichunk(-1),   // if ichunk is negative, then the first call to processor_get_chunk() will initialize
      ndrops(0), nprocessed(0)
{
    xassert(nc);
    nc->processor_start();
}


processor_handle::~processor_handle()
{
    nc->processor_end(ichunk);
    nc = shared_ptr<assembler_nerve_center> ();

    double dropfrac = 0.0;
    if ((ndrops > 0) || (nprocessed > 0))
	dropfrac = (double)ndrops / (double)(ndrops + nprocessed);
    
    stringstream ss;
    ss << name << ": " << nprocessed << " chunks processed, " 
       << ndrops << " dropped (dropfrac=" << dropfrac << ")\n";

    string s = ss.str();
    cout << s.c_str() << flush;
}


shared_ptr<assembled_chunk> processor_handle::get_next_chunk(thread_timer &timer)
{
    int nd = 0;
    shared_ptr<assembled_chunk> chunk = nc->processor_get_chunk(ichunk, nd, timer);

    if (chunk)
	nprocessed++;

    if (nd > 0) {
	cout << (string("  !!!! ") + name + " is running slow, can't keep up with assembler\n") << flush;
	ndrops += nd;
    }

    return chunk;
}


// -------------------------------------------------------------------------------------------------


void xmkdir(const string &dirname)
{
    int err = mkdir(dirname.c_str(), 0777);

    if (!err)
	return;

    if (errno != EEXIST) {
	stringstream ss;
	ss << "couldn't create directory " << dirname << ": " << strerror(errno);
	
	string err_msg = ss.str();
	cerr << err_msg << "\n";
	throw runtime_error(err_msg);
    }
    
    struct stat s;
    err = stat(dirname.c_str(), &s);
    
    if (err < 0) {
	stringstream ss;
	ss << "couldn't stat file " << dirname << ": " << strerror(errno);
	
	string err_msg = ss.str();
	cerr << err_msg << "\n";
	throw runtime_error(err_msg);
    }

    if (!S_ISDIR(s.st_mode)) {
	stringstream ss;
	ss << "couldn't create directory " << dirname << ": file already exists and is not a directory";
	
	string err_msg = ss.str();
	cerr << err_msg << "\n";
	throw runtime_error(err_msg);
    }
}


bool is_empty_dir(const string &dirname)
{
    DIR *dir = opendir(dirname.c_str());
    if (!dir)
	throw runtime_error(dirname + ": opendir() failed: " + strerror(errno));

    ssize_t name_max = pathconf(dirname.c_str(), _PC_NAME_MAX);
    name_max = min(name_max, (ssize_t)4096);

    vector<char> buf(sizeof(struct dirent) + name_max + 1);
    struct dirent *entry = reinterpret_cast<struct dirent *> (&buf[0]);
    
    struct dirent *result = nullptr;

    for (;;) {
	int err = readdir_r(dir, entry, &result);	
	if (err)
	    throw runtime_error(dirname + ": readdir_r() failed");
	if (!result)
	    return true;
	if (strcmp(entry->d_name, ".") && strcmp(entry->d_name, ".."))
	    return false;
    }
}


string make_dataset_name()
{
    char data_time[64];
    time_t rawtime;
    struct tm *timeinfo;

    time(&rawtime);
    timeinfo = gmtime(&rawtime);
    strftime(data_time, sizeof(data_time), "%Y%m%dT%H%M%SZ", timeinfo);

    // string dataset_name = string(data_time) + "_chime_beamformed";
    string dataset_name = string(data_time) + "_vdif_assembler";
    return dataset_name;
}


string make_data_dir(const string &dataset_name, int disk_id)
{
    xassert(disk_id >= 0);
    xassert(disk_id < constants::num_disks);

    string outdir = string("/drives/E/") + to_string(disk_id) + string("/") + dataset_name;
    // xmkdir(outdir);
    return outdir;
}


}   // namespace ch_vdif_assembler
