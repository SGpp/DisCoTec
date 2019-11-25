#ifndef NETWORK_UTILS_H
#define NETWORK_UTILS_H

#include <stdio.h>
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

    bool recvall(std::string& buff, size_t len, int flags = 0) const;

    bool recvall(char* buf, size_t len, int flags = 0)  const;

    bool recvallPrefixed(std::string& buff, int flags = 0) const;

    bool recvallPrefixed(char* buf, int flags = 0) const;

    template <typename FG_ELEMENT, typename FUNC>
    bool recvallBinaryAndReduceInPlace(FG_ELEMENT* const buff, size_t buffSize, FUNC reduceOp, size_t chunksize, int flags) const;

    template<typename FG_ELEMENT>
    bool recvallBinaryAndCorrectInPlace(FG_ELEMENT* const buff, size_t buffSize,
                                        size_t chunksize, int flags) const;

    template<typename FG_ELEMENT>
    bool recvallBinaryPrefixedCorrectInPlace(std::vector<FG_ELEMENT>& data,
        size_t chunksize = 2048, int flags = 0) const;

    bool recvallBinaryToFile(const std::string& filename, size_t len,
        size_t chunksize = 2048, int flags = 0) const;

    std::string getRemoteHost() const;

    int getRemotePort() const;

    bool isReadable(int timeoutSec = -1) const;

    bool isWriteable(int timeoutSec = -1) const;

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
    static const int noTimeout = -1;

    static bool forward(const ClientSocket& sender, const ClientSocket& receiver,
        size_t chunksize = 2048, size_t size = 0);

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
        // receive endianness
        char temp = ' ';
        assert(socket.recvall(&temp, 1)&& "Receiving Endianess failed");
        bool endianness = temp;
        //recv data
        char* raw = nullptr;
        size_t dataSize = len*sizeof(T);
        success = socket.recvall(raw, dataSize);
        if (!success)
          return false;
        buf.clear();
        buf.reserve(len);
        if (endianness == NetworkUtils::isLittleEndian()) { //same endianness
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
        delete[] raw;

        std::cout << "recvBinary finished!" << std::endl;
        return true;
      }

    template<typename T>
    inline static T reverseEndianness(const T& value)
    {
      size_t retValue;
      char* rawValue = reinterpret_cast<char*>(&value);
      char* rawRet = reinterpret_cast<char*>(retValue);

      // reverse bytes
      for (size_t j = 0; j < sizeof(T); j++)
          rawRet[j] = rawValue[sizeof(T)-j];

      return retValue;
    }
};

template<typename FG_ELEMENT>
bool ClientSocket::recvallBinaryAndCorrectInPlace(FG_ELEMENT* const buff, size_t buffSize, size_t chunksize, int flags) const {
  assert(isInitialized() && "Client Socket not initialized");
  // receive endianness
  char temp = ' ';
  assert(recvall(&temp, 1)&& "Receiving Endianess failed");
  bool hasSameEndianness = bool(temp) && NetworkUtils::isLittleEndian();

  // for recv()
  int err;
  ssize_t recvd = -1;
  size_t totalRecvd = 0;
  std::vector<char> recvBuff(chunksize);
  size_t rawSize = buffSize * sizeof(FG_ELEMENT);

  // receive binary data and correct endianness in place
  std::vector<char> remaining;
  remaining.reserve(sizeof(FG_ELEMENT));
  ssize_t head_k1 = 0;
  ssize_t head_k = 0;
  ssize_t tail_k = 0;
  ssize_t numAdded = 0;
  while(totalRecvd < rawSize) {
    recvd = recv(sockfd_, recvBuff.data(), std::min(rawSize-totalRecvd, chunksize), flags);

    err = errno;
    if (recvd <= 0)
      break;

    // recvBuff example for sizeof(FG_ELEMENT)=4:
    //
    //        recv k-1                     recv k
    //
    //    tail_k1          head_k1  tail_k            head_k
    //    v                v        v                 v
    //  [ oo|xxxx|...|xxxx|o ] -> [ ooo|xxxx|...|xxxx|oo ]
    //                              0                  recvd-1

    tail_k = (sizeof(FG_ELEMENT) - head_k1) % sizeof(FG_ELEMENT);
    head_k = (recvd-tail_k) % sizeof(FG_ELEMENT);

    // push tail into remaining
    for (long i = 0; i < tail_k; ++i)
      remaining.push_back(recvBuff[size_t(i)]);

    // if remaining has sizeof(FG_ELEMENT) entries, append to buff and clear
    // remaining
    if (remaining.size() == sizeof(FG_ELEMENT)) {
      FG_ELEMENT& remainingVal = *reinterpret_cast<FG_ELEMENT*>(remaining.data());
      if (hasSameEndianness)
        buff[numAdded++] =  remainingVal;
      else
        buff[numAdded++] = reverseEndianness(remainingVal);
      remaining.clear();
    }

    // append correct sized blocks
    const ssize_t&& numCorrect = recvd-head_k;
    for (ssize_t i = tail_k; i < numCorrect; i+=sizeof(FG_ELEMENT)) {
      if (hasSameEndianness) {
        FG_ELEMENT& val = *reinterpret_cast<FG_ELEMENT*>(recvBuff[i]);
        buff[numAdded++] = val;
      } else {
        FG_ELEMENT& val = *reinterpret_cast<FG_ELEMENT*>(recvBuff[i]);
        buff[numAdded++] = reverseEndianness(val);
        numAdded++;
      }
    }

    // push head into remaining
    for (ssize_t i = recvd-head_k; i < recvd; ++i)
      remaining.push_back(recvBuff[size_t(i)]);

    // overwrite old head
    head_k1 = head_k;

    totalRecvd += (size_t)recvd;
  }
  if ( recvd == -1) {
    std::cerr << "Receive failed: " << std::to_string(err) << std::endl;
    return false;
  } else if (recvd == 0) {
    std::cerr << "Receive failed, sender terminated too early: " <<
      std::to_string(err) << std::endl;
    return false;
  }
  return true;
}

template <typename FG_ELEMENT, typename FUNC>
bool ClientSocket::recvallBinaryAndReduceInPlace(FG_ELEMENT* const buff, size_t buffSize, FUNC reduceOp, size_t chunksize, int flags) const {
  assert(isInitialized() && "Client Socket not initialized");
  // receive endianness
  char temp = ' ';
  assert(recvall(&temp, 1)&& "Receiving Endianess failed");
  bool hasSameEndianness = bool(temp) && NetworkUtils::isLittleEndian();

  // for recv()
  int err;
  ssize_t recvd = -1;
  std::vector<char> recvBuff(chunksize);
  size_t totalRecvd = 0;
  size_t rawSize = buffSize * sizeof(FG_ELEMENT);

  // for oddly received data
  std::vector<char> remaining;
  remaining.reserve(sizeof(FG_ELEMENT));
  ssize_t head_k1 = 0;
  ssize_t head_k = 0;
  ssize_t tail_k = 0;
  ssize_t numAdded = 0;

  while(totalRecvd < rawSize) {
    recvd = recv(sockfd_, recvBuff.data(), std::min(rawSize-totalRecvd, chunksize), flags);
    err = errno;
    if (recvd <= 0)
      break;

    // recvBuff example for sizeof(FG_ELEMENT)=4:
    //
    //        recv k-1                     recv k
    //
    //    tail_k1          head_k1  tail_k            head_k
    //    v                v        v                 v
    //  [ oo|xxxx|...|xxxx|o ] -> [ ooo|xxxx|...|xxxx|oo ]
    //                              0                  recvd-1

    tail_k = (sizeof(FG_ELEMENT) - head_k1) % sizeof(FG_ELEMENT);
    head_k = (recvd-tail_k) % sizeof(FG_ELEMENT);

    // push tail into remaining
    for (long i = 0; i < tail_k; ++i)
      remaining.push_back(recvBuff[size_t(i)]);

    // if remaining has sizeof(FG_ELEMENT) entries, add to buff and clear
    // remaining
    if (remaining.size() == sizeof(FG_ELEMENT)) {
      FG_ELEMENT remainingVal = *reinterpret_cast<FG_ELEMENT*>(remaining.data());
      if (hasSameEndianness)
        buff[numAdded++] = reduceOp(buff[numAdded], remainingVal);
      else
        buff[numAdded++] = reduceOp(buff[numAdded], reverseEndianness(remainingVal));
      remaining.clear();
    }

    // add correctly received data
    const ssize_t&& numCorrect = recvd-head_k;
    for (ssize_t i = tail_k; i < numCorrect; i+=sizeof(FG_ELEMENT)) {
      FG_ELEMENT val = *reinterpret_cast<FG_ELEMENT*>(recvBuff[i]);
      if (hasSameEndianness)
        buff[numAdded++] = reduceOp(buff[numAdded], val);
      else
        buff[numAdded++] = reduceOp(buff[numAdded], reverseEndianness(val));
    }

    // push head into remaining
    for (ssize_t i = recvd-head_k; i < recvd; ++i)
      remaining.push_back(recvBuff[size_t(i)]);

    // overwrite old head
    head_k1 = head_k;
  }

  if ( recvd == -1) {
    std::cerr << "Receive failed: " << std::to_string(err) << std::endl;
    return false;
  } else if (recvd == 0) {
    std::cerr << "Receive failed, sender terminated too early: " <<
      std::to_string(err) << std::endl;
    return false;
  }
  return true;
}

/** Sends the given buffer in binary form.
 * The first byte sent specifies the endianess of the system and is followed by
 * the raw binary data.
 */
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

#endif
