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
#include <memory>

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

    bool sendallPrefixed(const char* buf, size_t len) const;

    template<typename FG_ELEMENT>
    bool sendallBinary(const std::vector<FG_ELEMENT>& buf, int flags = 0) const;

    template<typename FG_ELEMENT>
    bool sendallBinaryPrefixed(const std::vector<FG_ELEMENT>& buf, int flags = 0) const;

    bool recvall(std::string& buff, size_t len, int flags = 0) const;

    bool recvall(char* buf, size_t len, int flags = 0)  const;

    bool recvallPrefixed(std::string& buff, int flags = 0) const;

    bool recvallPrefixed(char* buf, int flags = 0) const;

    template<typename FG_ELEMENT>
    bool recvallBinary(std::vector<FG_ELEMENT>& buff, bool& endianness, int flags = 0) const;

    template<typename FG_ELEMENT>
    bool recvallBinaryPrefixed(std::vector<FG_ELEMENT>& buff, bool& endiannes,
        int flags = 0) const;

    template<typename FG_ELEMENT>
    bool recvallBinaryPrefixedCorrectInPlace(std::vector<FG_ELEMENT>& data,
        size_t chunksize = 2048, int flags = 0) const;

    bool recvallBinaryToFile(const std::string& filename, size_t len,
        size_t chunksize = 2048, int flags = 0) const;

    std::string getRemoteHost() const;

    int getRemotePort() const;

    bool isReadable(int timeout = -1) const;

    bool isWriteable(int timeout = -1) const;

  private:
    ClientSocket(); // for ServerSocket only

    bool sendSize();

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

    // TODO Optimize same endianess: no need for copying values to buf.
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
          for (size_t i = 0; i < len; i++) {
            buf.push_back(values[i]);
          }
        } else { // different endianness
          char* rightOrder = new char[sizeof(T)];
          // extract values from raw buffer in right order
          // reverts byte order for each sizeof(T) big block
          for (size_t i = 0; i < len; i+=sizeof(T)) {
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

    template<typename T>
    static void reverseEndianness(std::vector<T>& buf)
    {
      char* rawBuf = reinterpret_cast<char*>(buf.data());
      size_t rawSize = buf.size() * sizeof(T);
      for (size_t i = 0; i < rawSize; i+=sizeof(T))
      {
        // swap entrys of block
        for (size_t j = 0; j < sizeof(T)/2; j++) {
          rawBuf[i+j] ^= rawBuf[i+(sizeof(T))-j];
          rawBuf[i+(sizeof(T))-j] ^= rawBuf[i+j];
          rawBuf[i+j] ^= rawBuf[i+(sizeof(T))-j];
        }
      }
    }
};

template<typename FG_ELEMENT>
bool ClientSocket::recvallBinaryPrefixedCorrectInPlace(std::vector<FG_ELEMENT>& buff, size_t chunksize, int flags) const {
  assert(isInitialized() && "Client Socket not initialized");
  // receive endianness
  char temp = ' ';
  assert(recvall(&temp, 1)&& "Receiving Endianess failed");
  bool endianness = temp;

  // receive length
  int err;
  long n = -1;
  std::string lenstr = "";
  do {
    n = recv(sockfd_, &temp, 1, flags);
    err = errno;
    if (n <= 0)
      break;
    lenstr += temp;
  } while (temp != '#');
  if (n == -1) {
    std::cerr << "Receive of length failed: " << std::to_string(err)
      << std::endl;
    return false;
  } else if (n == 0) {
    std::cerr << "Receive of length failed, sender terminated too early: " <<
      std::to_string(err) << std::endl;
    return false;
  }
  lenstr.pop_back();
  assert(NetworkUtils::isInteger(lenstr) && "Received length is not a number");
  size_t len = (size_t) std::stoi(lenstr);

  buff.resize(len);
  char* rawBuffer = reinterpret_cast<char*>(buff.data());

  // receive binary data and correct endianness in place
  std::vector<char> tempBuff(chunksize);
  std::vector<char> overhead(sizeof(FG_ELEMENT));
  size_t total = 0;
  size_t oldSize = 0;
  size_t newSize = 0;
  while(total < len) {
    n = recv(sockfd_, tempBuff.data(), len-total, flags);

    err = errno;
    if (n <= 0)
      break;

    // tempBuff example for sizeof(FG_ELEMENT)=4:
    //
    //        recv k-1                     recv k
    //                     old     tail               new
    //  [ oo|xxxx|...|xxxx|o ] -> [ ooo|xxxx|...|xxxx|oo ]
    //                              0                  n-1
    if (total == 0) //first run
      oldSize = n % sizeof(FG_ELEMENT);
    else
      oldSize = newSize;

    if (endianness != NetworkUtils::isLittleEndian()) {
      size_t tailSize = sizeof(FG_ELEMENT) - oldSize;
      size_t rawSize = total + tailSize;
      newSize = (tailSize - n) % sizeof(FG_ELEMENT);
      size_t blocksSize = (n - tailSize - newSize);
      if (tailSize > 0) {
        // fill remaining block
        for (size_t i = 0; i < tailSize; i++)
          rawBuffer[total-oldSize-1 + i] = tempBuff[i];
      }
      if (blocksSize > 0) {
        // append correct sized blocks
        for (size_t i = rawSize; i < rawSize+n-newSize; i+=sizeof(FG_ELEMENT)) {
          for (size_t j = sizeof(FG_ELEMENT)-1; j >= 0; j--) {
            rawBuffer[i+(sizeof(FG_ELEMENT)-1-j)] = tempBuff[i-rawSize+j];
          }
        }
      }
      if (newSize > 0) {
        rawSize += n-newSize;
        // append new overhead
        for (size_t i = rawSize; i < rawSize+newSize; i++){
          size_t j = n - (i-rawSize) -1;
          rawBuffer[i] = tempBuff[j];
        }
      }
    }

    total += n;
  }
  if ( n == -1) {
    std::cerr << "Receive failed: " << std::to_string(err) << std::endl;
    return false;
  } else if (n == 0) {
    std::cerr << "Receive failed, sender terminated too early: " <<
      std::to_string(err) << std::endl;
    return false;
  }
  return true;
}

template<typename FG_ELEMENT>
bool ClientSocket::recvallBinaryPrefixed(std::vector<FG_ELEMENT>& buff, bool& endianness, int flags) const
{
  assert(isInitialized() && "Client Socket not initialized");

  // receive endianness
  char temp = ' ';
  assert(recvall(&temp, 1)&& "Receiving Endianess failed");
  endianness = temp;

  // receive binary data
  char* rawBuf = reinterpret_cast<char*>(buff.data());
  return recvallPrefixed(rawBuf, flags);

}

template<typename FG_ELEMENT>
bool ClientSocket::recvallBinary(std::vector<FG_ELEMENT>& buff,
    bool& endianness, int flags) const {
  // receive endianness
  char temp = ' ';
  assert(recvall(&temp, 1)&& "Receiving Endianess failed");
  endianness = temp;

  // receive binary data
  char* rawBuf = reinterpret_cast<char*>(buff.data());
  size_t rawSize = buff.size()*sizeof(FG_ELEMENT);
  std::cout << "receiving " << rawSize << "Bytes" << std::endl;
  return recvall(rawBuf, rawSize, flags);
}


template<typename FG_ELEMENT>
bool ClientSocket::sendallBinary(const std::vector<FG_ELEMENT>& buf, int flags) const {
  const char* rawBuf = reinterpret_cast<const char*>(buf.data());
  size_t rawSize = buf.size() * sizeof(FG_ELEMENT);
  std::cout << "sending " << rawSize << "Bytes" << std::endl;

  char isLittleEndian = static_cast<char>(NetworkUtils::isLittleEndian());
  bool success = sendall(&isLittleEndian, 1);
  if (!success)
    return false;
  return sendall(rawBuf, rawSize);
}

// Sends the given buffer in binary form. The first byte sent specifies
// if the system is little-endian. The next bytes specify the length of
// the following buffer and are sperated by an # from the following raw data. 
template<typename FG_ELEMENT>
bool ClientSocket::sendallBinaryPrefixed(const std::vector<FG_ELEMENT>& buf, int flags) const {
  const char* rawBuf = reinterpret_cast<const char*>(buf.data());
  size_t rawSize = buf.size() * sizeof(FG_ELEMENT);

  char isLittleEndian = static_cast<char>(NetworkUtils::isLittleEndian());
  bool success = sendall(&isLittleEndian, 1);
  if (!success)
    return false;
  return sendallPrefixed(rawBuf, rawSize);
}


#endif
