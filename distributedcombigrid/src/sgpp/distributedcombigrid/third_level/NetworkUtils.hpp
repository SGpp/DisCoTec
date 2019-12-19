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

template <typename FG_ELEMENT>
using ReduceFcn = FG_ELEMENT(*)(const FG_ELEMENT&, const FG_ELEMENT&);

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

    bool sendall(const char* buff, size_t len) const;

    bool sendallPrefixed(const std::string& mesg) const;

    bool sendallPrefixed(const char* buff, size_t len) const;

    template<typename FG_ELEMENT>
    bool sendallBinary(const FG_ELEMENT* buff, size_t buffSize, int flags = 0) const;

    bool recvall(std::string& buff, size_t len, int flags = 0) const;

    bool recvall(char* buff, size_t len, int flags = 0)  const;

    bool recvallPrefixed(std::string& buff, int flags = 0) const;

    bool recvallPrefixed(std::unique_ptr<char[]>& buff, size_t& len,
                         int flags = 0) const;

    template <typename FG_ELEMENT>
    bool recvallBinaryAndReduceInPlace(FG_ELEMENT* buff,
                 size_t buffSize,
                 ReduceFcn<FG_ELEMENT> reduceOp,
                 size_t chunksize = 2048, int flags = 0) const;

    template<typename FG_ELEMENT>
    bool recvallBinaryAndCorrectInPlace(FG_ELEMENT* buff, size_t buffSize,
                                        size_t chunksize = 2048,
                                        int flags = 0) const;

    // TODO
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
    inline static T reverseEndianness(const T& value)
    {
      T retValue;
      const char* rawValue = reinterpret_cast<const char*>(&value);
      char* rawRet = reinterpret_cast<char*>(&retValue);

      // reverse bytes
      for (size_t j = 0; j < sizeof(T); j++)
          rawRet[j] = rawValue[sizeof(T)-j];

      return retValue;
    }
};

template<typename FG_ELEMENT>
bool ClientSocket::recvallBinaryAndCorrectInPlace(FG_ELEMENT* buff,
                           size_t buffSize, size_t chunksize, int flags) const {
  assert(isInitialized() && "Client Socket not initialized");
  // receive endianness
  char temp = ' ';
  assert(recvall(&temp, 1) && "Receiving Endianess failed");
  bool hasSameEndianness = bool(temp) == NetworkUtils::isLittleEndian();

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
  ssize_t numAppended = 0;
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
    head_k = (recvd - tail_k) % sizeof(FG_ELEMENT);

    // push tail into remaining
    for (long i = 0; i < tail_k; ++i)
      remaining.push_back(recvBuff[size_t(i)]);

    // if remaining has sizeof(FG_ELEMENT) entries, append to buff and clear
    // remaining
    if (remaining.size() == sizeof(FG_ELEMENT)) {
      FG_ELEMENT& remainingVal = *reinterpret_cast<FG_ELEMENT*>(remaining.data());
      if (hasSameEndianness)
        buff[numAppended++] =  remainingVal;
      else
        buff[numAppended++] = NetworkUtils::reverseEndianness(remainingVal);
      remaining.clear();
    }

    // append correct sized blocks
    for (ssize_t i = tail_k; i < recvd - head_k; i+=sizeof(FG_ELEMENT)) {
      if (hasSameEndianness) {
        FG_ELEMENT& val = *reinterpret_cast<FG_ELEMENT*>(&recvBuff[i]);
        buff[numAppended++] = val;
      } else {
        FG_ELEMENT& val = *reinterpret_cast<FG_ELEMENT*>(&recvBuff[i]);
        buff[numAppended++] = NetworkUtils::reverseEndianness(val);
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

template <typename FG_ELEMENT>
bool ClientSocket::recvallBinaryAndReduceInPlace(FG_ELEMENT* buff,
                size_t buffSize,
                FG_ELEMENT (*reduceOp) (const FG_ELEMENT&, const FG_ELEMENT&),
                size_t chunksize, int flags) const {
  assert(isInitialized() && "Client Socket not initialized");
  // receive endianness
  char temp = ' ';
  assert(recvall(&temp, 1)&& "Receiving Endianess failed");
  bool hasSameEndianness = bool(temp) == NetworkUtils::isLittleEndian();

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
  ssize_t numReduced = 0;

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
        buff[numReduced++] = reduceOp(buff[numReduced], remainingVal);
      else
        buff[numReduced++] = reduceOp(buff[numReduced], NetworkUtils::reverseEndianness(remainingVal));
      remaining.clear();
    }

    // add correctly received data
    const ssize_t&& numCorrect = recvd-head_k;
    for (ssize_t i = tail_k; i < numCorrect; i+=sizeof(FG_ELEMENT)) {
      FG_ELEMENT val = *reinterpret_cast<FG_ELEMENT*>(&recvBuff[i]);
      if (hasSameEndianness)
        buff[numReduced++] = reduceOp(buff[numReduced], val);
      else
        buff[numReduced++] = reduceOp(buff[numReduced], NetworkUtils::reverseEndianness(val));
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

/** Sends the given buffer in binary form.
 * The first byte sent specifies the endianess of the system and is followed by
 * the raw binary data.
 */
template<typename FG_ELEMENT>
bool ClientSocket::sendallBinary(const FG_ELEMENT* buff, size_t buffSize, int flags) const {
  const char* rawBuf = reinterpret_cast<const char*>(buff);
  size_t rawSize = buffSize * sizeof(FG_ELEMENT);

  char endianFlag = static_cast<char>(NetworkUtils::isLittleEndian());
  bool success = sendall(&endianFlag, 1);
  if (!success)
    return false;
  return sendall(rawBuf, rawSize);
}

#endif
