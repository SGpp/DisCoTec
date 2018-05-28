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
#include <vector>
#include <fstream>

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
    virtual ~Socket();
    virtual bool init() = 0;
    bool isInitialized() const;
    int getFileDescriptor() const;
  protected:
    int sockfd_ = -1;
    bool initialized_ = false;
};

class ClientSocket : public Socket {
  friend class ServerSocket; // allows the server to accept a client

  public:
    ClientSocket(const std::string& host,
          const int port);
    ~ClientSocket();
    bool init();
    bool sendall(const std::string& mesg) const;
    bool sendall(const char* buf, size_t len) const;
    bool sendallPrefixed(const std::string& mesg) const;
    bool recvall(std::string& buff, size_t len, int flags = 0) const;
    bool recvall(char* &buf, size_t len, int flags = 0) const;
    bool recvallBinaryToFile(const std::string& filename, size_t len, 
        size_t chunksize = 2048) const;
    bool recvallPrefixed(std::string& buff, int flags = 0) const;
    std::string getRemoteHost() const;
    unsigned short getRemotePort() const;
    bool isReadable(int timeout = -1) const;
    bool isWriteable(int timeout = -1) const;

  private:
    ClientSocket(); // for ServerSocket only
    int remotePort_;
    std::string remoteHost_;
};

class ServerSocket : public Socket {
  public:
    ServerSocket();
    ServerSocket(const unsigned short port);
    ~ServerSocket();
    bool init();
    ClientSocket* acceptClient() const;
    int getPort();
  private:
    int port_;

};

class NetworkUtils {
  public:
    static bool tunnel(const ClientSocket* sender, const ClientSocket* receiver, unsigned int buffsize = 16384);

    static bool isInteger(const std::string& s);

    static bool isLittleEndian();

    template<typename T>
      static bool sendBinary(const std::vector<T>& buf, const ClientSocket& socket);

    template<typename T>
      static bool recvBinary(std::vector<T>& buf, size_t len, const ClientSocket& socket);
};

#endif
