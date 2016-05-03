#include <sys/socket.h>
#include <sys/epoll.h>
#include <arpa/inet.h>
#include <errno.h>
#include <cstring>

#include "ch_vdif_assembler_internals.hpp"

using namespace std;

namespace ch_vdif_assembler {
#if 0
};  // pacify emacs c-mode!
#endif


struct network_thread : public thread_base 
{
    shared_ptr<assembler_nerve_center> nc;

    network_thread(const shared_ptr<assembler_nerve_center> &nc_)
	: thread_base("network thread"),
	  nc(nc_)
    {
	xassert(nc);
    }

    virtual ~network_thread() { }

    // stream_start() gets called prior to spawning thread, but thread is responsible for calling stream_end()
    virtual void thread_body()
    {
	static const int packets_per_chunk = 50000;
	static const int packet_size = constants::packet_nbytes;
	static const int max_events = 100;
	static const int udp_port_number = 10251;
	
	// kill assembler if we throw an exception somewhere
	assembler_killer killer(nc, "network thread threw exception");
	
	// the "false" in the constructor is "set_zero=false"
	shared_ptr<vdif_chunk_pool> vpool = make_shared<vdif_chunk_pool> (packets_per_chunk, false);

	struct epoll_event events[max_events];
	struct epoll_event ev;

	int epoll_fd = epoll_create(10);

	if (epoll_fd < 0) {
	    cout << "network thread: epoll_create() failed: " << strerror(errno) << endl;
	    throw runtime_error(strerror(errno));
	}

	int sock_fd = socket(AF_INET, SOCK_DGRAM | SOCK_NONBLOCK, IPPROTO_UDP);
	
	if (sock_fd < 0) {
	    cout << "network thread: socket() failed: " << strerror(errno) << endl;
	    throw runtime_error(strerror(errno));
	}
	
	struct sockaddr_in server_address;
	memset(&server_address, 0, sizeof(server_address));
	
	server_address.sin_family = AF_INET;
	inet_pton(AF_INET, "0.0.0.0", &server_address.sin_addr);
	server_address.sin_port = htons(udp_port_number);
	
	if (::bind(sock_fd, (struct sockaddr *) &server_address, sizeof(server_address)) < 0) {
	    cout << "network thread: bind() failed: " << strerror(errno) << endl;
	    throw runtime_error(strerror(errno));
	}
	
	int n = 256 * 1024 * 1024;
	if (setsockopt(sock_fd, SOL_SOCKET, SO_RCVBUF, (void *) &n, sizeof(n)) < 0) {
	    cout << "network thread: setsockopt() failed: " << strerror(errno) << endl;
	    throw runtime_error(strerror(errno));
	}
	
	ev.events = EPOLLIN;
	ev.data.fd = sock_fd;
	if (epoll_ctl(epoll_fd, EPOLL_CTL_ADD, sock_fd, &ev) < 0) {
	    cout << "network thread: epoll_ctl() failed: " << strerror(errno) << endl;
	    throw runtime_error(strerror(errno));
	}

	int seq_id = 0;
	struct timeval tv0 = get_time();

	shared_ptr<vdif_chunk> chunk = make_shared<vdif_chunk> (vpool, seq_id);
	cout << "network thread: start\n";
	
	for (;;) {
	    int num_events = epoll_wait(epoll_fd, events, max_events, -1);
	    
	    if (num_events < 0) {
		cout << "network thread: epoll_wait() failed: " << strerror(errno) << endl;
		throw runtime_error(strerror(errno));
	    }
	    
	    for (int i = 0; i < num_events; i++) {
		if (events[i].data.fd != sock_fd)
		    continue;
		
		uint8_t *buf = chunk->buf + chunk->size * packet_size;
		ssize_t bytes_read = read(sock_fd, buf, packet_size);
		
		// FIXME silently drop?
		if (bytes_read != packet_size)
		    continue;
		
		chunk->size++;
		if (chunk->size < chunk->capacity)
		    continue;
		
		// Last packet in chunk received
		struct timeval tv1 = get_time();
		double gbps = 8.0e-9 * packet_size * chunk->capacity / time_diff(tv0,tv1);
		cout << (to_string(gbps) + " Gbps\n") << flush;

		nc->stream_put_chunk(chunk, timer);
		tv0 = tv1;
		seq_id++;

		chunk = make_shared<vdif_chunk> (vpool, seq_id);		
	    }
	}

	nc->stream_end();
	killer.let_live();
    }
};


struct network_stream : public vdif_stream
{
    network_stream()
	: vdif_stream(true)   // is_realtime = true
    { }

    virtual ~network_stream() { }

    virtual void spawn_threads(const shared_ptr<assembler_nerve_center> &nc)
    {
	xassert(nc);
	nc->check_alive();
	spawn_thread<network_thread> (nc);
    }
};


shared_ptr<vdif_stream> make_network_stream()
{
    return make_shared<network_stream> ();
}


}   // namespace ch_vdif_assembler
