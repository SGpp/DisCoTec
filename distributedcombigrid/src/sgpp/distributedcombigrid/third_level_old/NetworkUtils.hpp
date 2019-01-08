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
    virtual ~ClientSocket();
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
    int getRemotePort() const;
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
    virtual ~ServerSocket();
    ServerSocket(const unsigned short port);
    bool init();
    ClientSocket* acceptClient() const;
    int getPort();
  private:
    int port_;
};

class NetworkUtils {
  public:
    static bool forward(const ClientSocket* sender, const ClientSocket* receiver,
        size_t size = 0, size_t  chunksize = 2048);

    static bool isInteger(const std::string& s);

    static bool isLittleEndian();

    static void split(const std::string& s, const char c,
        std::vector<std::string>& tokens);

    template<typename T>
      static bool sendBinary(const std::vector<T>& buf, const ClientSocket& socket)
      {
        if (!socket.isInitialized())
          return false;
        bool success = false;
        size_t dataSize = buf.size() * sizeof(T);
        char endianness = NetworkUtils::isLittleEndian();
        // send endianness header
        success = socket.sendall(&endianness, 1);
        if (!success)
          return false;
        success = socket.sendall( (char*)(buf.data()), dataSize );
        if (!success)
          return false;
        std::cout << "sendBinary finished!" << std::endl;
        return true;
      }

    // TODO Optimize memory usage (hopefully this wont cause thrashing and os recognizes the linear access)
    // Optimize same endianess: no need for copying values to buf.
    template<typename T>
      static bool recvBinary(std::vector<T>& buf, size_t len, const ClientSocket& socket)
      {
        if (!socket.isInitialized())
          return false;
        bool success = false;
        //recv endianness first
        char* endianness = nullptr;
        success = socket.recvall(endianness, 1);
        if (!success)
          return false;
        //recv data
        char* raw = nullptr;
        size_t dataSize = len*sizeof(T);
        success = socket.recvall(raw, dataSize);
        if (!success)
          return false;
        buf.clear();
        buf.reserve(len);
        if ((bool) *endianness == NetworkUtils::isLittleEndian()) { //same endianness
          T* values = reinterpret_cast<T*>(raw);
          // copy values to buf
          for (size_t i = 0; i < len; i++) { // linear access to poss. huge array
            buf.push_back(values[i]);
          }
        } else { // different endianness
          char* rightOrder = new char[sizeof(T)];
          // extract values from raw buffer in right order
          // reverts byte order for each sizeof(T) big block
          for (size_t i = 0; i < len; i+=sizeof(T)) { // somehow linear access
            for (size_t j = sizeof(T)-1; j >= 0; j--) {
              rightOrder[sizeof(T)-1 - j] = raw[i+j];
              buf.push_back( reinterpret_cast<T &>(rightOrder) );
            }
          }
          delete[] rightOrder;
        }
        delete[] endianness;
        delete[] raw;

        std::cout << "recvBinary finished!" << std::endl;
        return true;
      }
};

#endif
