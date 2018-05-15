#ifndef NETWORK_UTILS_H
#define NETWORK_UTILS_H

#include <string>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <netinet/in.h>
#include <unistd.h>
#include <netdb.h>
#include <strings.h>
#include <csignal>
#include <errno.h>
#include <unistd.h>
#include <assert.h>
#include <arpa/inet.h>
#include <iostream>

/*
 * C++ abstraction layer to simplify the usage of
 * TCP Internet Sockets.
 * Additionally some higher level functionality
 * is given in the NetworkUtils class.
 * Attention: there is no crypto stuff implemented
 * here. Man in the middle and other attacks are
 * possible.
 */
class Socket {
  public:
    Socket();
    virtual bool init() = 0;
    bool isInitialized() const;
    int getFileDescriptor() const;
  protected:
    int sockfd = -1;
    bool initialized = false;
};

class ClientSocket : public Socket {
  friend class ServerSocket; // allows the server to accept a client

  public:
    ClientSocket(const std::string& host,
          const int port);
    bool init();
    bool sendall(const std::string& mesg) const;
    bool sendallPrefixed(const std::string& mesg) const;
    bool recvall(std::string& buff, size_t len, int flags = 0) const;
    bool recvallPrefixed(std::string& buff, int flags = 0) const;
    std::string getRemoteHost() const;
    unsigned short getRemotePort() const;
    bool isReadable(int timeout = -1) const;
    bool isWriteable(int timeout = -1) const;

  private:
    int remotePort;
    std::string remoteHost;
};

class ServerSocket : public Socket {
  unsigned short port;

  public:
    ServerSocket(const unsigned short port);
    bool init();
    ClientSocket acceptClient() const;
};

class NetworkUtils {
  public:
    static bool tunnel(ClientSocket sender, ClientSocket receiver, unsigned int buffsize = 16384);
    static bool isInteger(const std::string& s);
};

#endif
